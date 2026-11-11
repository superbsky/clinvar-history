"""Command-line entry point for the ClinVar history pipeline."""

from datetime import date
from pathlib import Path

import click
import httpx

from .attribution import attribute_event_tables
from .manifest import (
    default_index_urls,
    discover_manifest,
    read_manifest,
    write_manifest,
)
from .pipeline import (
    http_vcf_source,
    ingest_submission_releases,
    ingest_vcf_releases,
    resolve_release_checksums,
)


@click.group()
def main() -> None:
    """Build compact longitudinal ClinVar datasets."""


@main.command("build-manifest")
@click.option("--start", type=click.DateTime(formats=["%Y-%m-%d"]), required=True)
@click.option("--end", type=click.DateTime(formats=["%Y-%m-%d"]), required=True)
@click.option(
    "--output",
    type=click.Path(path_type=Path),
    default=Path("data/manifest/releases.parquet"),
    show_default=True,
)
def build_manifest(start, end, output: Path) -> None:
    """Discover official releases and write a bounded Parquet manifest."""
    start_date: date = start.date()
    end_date: date = end.date()
    if end_date < start_date:
        raise click.BadParameter("--end must not precede --start", param_hint="--end")

    with httpx.Client(timeout=30, follow_redirects=True) as client:
        releases = discover_manifest(
            default_index_urls(start_date, end_date),
            start=start_date,
            end=end_date,
            client=client,
        )
        releases = resolve_release_checksums(releases, client=client)
    write_manifest(releases, output)
    from .checkpoints import stamp_checkpoint
    from .runner import manifest_request

    stamp_checkpoint(output, manifest_request(start_date, end_date))
    click.echo(f"Wrote {len(releases)} releases to {output}")


@main.command("ingest-vcf")
@click.option(
    "--manifest",
    "manifest_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "--root",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("data"),
    show_default=True,
)
@click.option(
    "--keep-raw",
    is_flag=True,
    default=False,
    help="Retain downloaded .vcf.gz files (disabled by default).",
)
def ingest_vcf(manifest_path: Path, root: Path, keep_raw: bool) -> None:
    """Stream VCF releases into compact state and event tables."""
    releases = read_manifest(manifest_path)
    with httpx.Client(timeout=None, follow_redirects=True) as client:
        releases = resolve_release_checksums(releases, client=client)

        def source(release):
            return http_vcf_source(
                release,
                client=client,
                keep_raw=keep_raw,
                raw_dir=root / "raw" / "vcf",
                expected_md5=release.md5,
            )

        summary = ingest_vcf_releases(releases, root=root, source=source)
    click.echo(
        f"Processed {summary.processed} releases; skipped {summary.skipped} completed"
    )


@main.command("ingest-submissions")
@click.option(
    "--manifest",
    "manifest_path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "--root",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("data"),
    show_default=True,
)
@click.option(
    "--keep-raw",
    is_flag=True,
    default=False,
    help="Retain downloaded .txt.gz files (disabled by default).",
)
def ingest_submissions(manifest_path: Path, root: Path, keep_raw: bool) -> None:
    """Stream submission summaries into compact SCV state and events."""
    releases = read_manifest(manifest_path)
    with httpx.Client(timeout=None, follow_redirects=True) as client:
        releases = resolve_release_checksums(releases, client=client)

        def source(release):
            return http_vcf_source(
                release,
                client=client,
                keep_raw=keep_raw,
                raw_dir=root / "raw" / "submissions",
                expected_md5=release.md5,
            )

        summary = ingest_submission_releases(
            releases,
            root=root,
            source=source,
        )
    click.echo(
        f"Processed {summary.processed} releases; skipped {summary.skipped} completed"
    )


@main.command("attribute-events")
@click.option(
    "--root",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("data"),
    show_default=True,
)
@click.option(
    "--output",
    type=click.Path(path_type=Path, dir_okay=False),
    default=Path("results/attributions.parquet"),
    show_default=True,
)
@click.option(
    "--max-lag-days",
    type=click.IntRange(min=0),
    default=35,
    show_default=True,
    help="Maximum observation-date lag after interval-overlap filtering.",
)
def attribute_events(root: Path, output: Path, max_lag_days: int) -> None:
    """Build conservative candidate attributions for aggregate changes."""
    attributions = attribute_event_tables(
        root,
        output,
        max_lag_days=max_lag_days,
    )
    click.echo(f"Wrote {len(attributions)} attributions to {output}")


@main.command("run-study")
@click.option("--root", type=click.Path(path_type=Path), default=Path("data/full"))
@click.option(
    "--results", type=click.Path(path_type=Path), default=Path("results/full-v1")
)
@click.option("--start", type=click.DateTime(formats=["%Y-%m-%d"]), required=True)
@click.option("--end", type=click.DateTime(formats=["%Y-%m-%d"]), required=True)
@click.option(
    "--offline",
    is_flag=True,
    help="Use completed local data; never access the network.",
)
@click.option("--refresh-manifest", is_flag=True)
@click.option("--max-lag-days", type=click.IntRange(min=0), default=35)
@click.option("--attempts", type=click.IntRange(min=1), default=3)
@click.option("--retry-delay", type=click.FloatRange(min=0, max=60), default=10)
def run_study(root, results, start, end, **options):
    """Run ingestion and analyses with semantic checkpoint validation."""
    from .runner import run_study as run

    try:
        run(root.resolve(), results.resolve(), start.date(), end.date(), **options)
    except ValueError as error:
        raise click.ClickException(str(error)) from error


@main.command("adopt-legacy")
@click.option("--root", type=click.Path(path_type=Path, exists=True), required=True)
@click.option("--results", type=click.Path(path_type=Path, exists=True), required=True)
@click.option("--start", type=click.DateTime(formats=["%Y-%m-%d"]), required=True)
@click.option("--end", type=click.DateTime(formats=["%Y-%m-%d"]), required=True)
def adopt_legacy(root, results, start, end):
    """Audit and register retained v0.1 artifacts for managed updates."""
    from .runner import adopt_legacy as adopt

    try:
        adopt(root, results, start.date(), end.date())
    except ValueError as error:
        raise click.ClickException(str(error)) from error


@main.command("review-xml")
@click.option("--root", type=click.Path(path_type=Path, exists=True), required=True)
@click.option(
    "--attributions", type=click.Path(path_type=Path, exists=True), required=True
)
@click.option("--output", type=click.Path(path_type=Path), required=True)
@click.option("--limit", type=click.IntRange(1, 100), default=12)
@click.option("--seed", type=int, default=20261111)
@click.option("--max-record-bytes", type=click.IntRange(min=1), default=2_000_000)
@click.option("--max-total-bytes", type=click.IntRange(min=1), default=20_000_000)
@click.option("--offline", is_flag=True)
def review_xml(**options):
    """Review a bounded sample of unexplained transitions against current XML."""
    from .xml_review import review_sample

    try:
        result = review_sample(**options)
    except ValueError as error:
        raise click.ClickException(str(error)) from error
    click.echo(
        f"XML sample: {result['status_counts']}; {result['total_xml_bytes']} bytes"
    )


@main.command("compare-xml-versions")
@click.option("--variation-id", type=click.IntRange(min=1), required=True)
@click.option("--old-version", type=click.IntRange(min=1), required=True)
@click.option("--new-version", type=click.IntRange(min=1), required=True)
@click.option("--output", type=click.Path(path_type=Path), required=True)
@click.option("--max-record-bytes", type=click.IntRange(min=1), default=2_000_000)
@click.option("--offline", is_flag=True)
def compare_xml_versions(**options):
    """Compare SCV evidence in two explicitly versioned VCV XML records."""
    from .xml_review import compare_versions

    try:
        result = compare_versions(**options)
    except ValueError as error:
        raise click.ClickException(str(error)) from error
    click.echo(
        f"Compared versions: {len(result['scv_changes'])} SCV differences; historical causality unconfirmed"
    )


if __name__ == "__main__":
    main()
