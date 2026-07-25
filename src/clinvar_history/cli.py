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


if __name__ == "__main__":
    main()
