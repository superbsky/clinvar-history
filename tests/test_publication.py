import re
import subprocess
import tomllib
from pathlib import Path

import yaml
from markdown_it import MarkdownIt

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_URL = "https://github.com/superbsky/clinvar-history"
PINNED_ACTION = re.compile(r"[^@\s]+@[0-9a-f]{40}")
MARKDOWN = MarkdownIt("commonmark")


def _inline_markdown_targets(document: Path):
    for block_token in MARKDOWN.parse(document.read_text()):
        for token in (block_token, *(block_token.children or ())):
            if token.type == "link_open":
                yield token.attrGet("href")
            elif token.type == "image":
                yield token.attrGet("src")


def _unresolved_local_markdown_links(markdown_files, root: Path):
    resolved_root = root.resolve()
    unresolved = []

    for document in markdown_files:
        document_label = document.relative_to(resolved_root)
        for raw_target in _inline_markdown_targets(document):
            if raw_target.lower().startswith(("#", "http://", "https://", "mailto:")):
                continue
            target = raw_target.split("#", 1)[0]
            if not target:
                continue
            failure = f"{document_label} -> {target}"
            target_path = Path(target)
            if target_path.is_absolute():
                unresolved.append(failure)
                continue
            resolved_target = (document.parent / target_path).resolve()
            if not resolved_target.is_relative_to(resolved_root):
                unresolved.append(failure)
                continue
            if not resolved_target.exists():
                unresolved.append(failure)

    return unresolved


def _readme_section(readme: str, title: str) -> str:
    match = re.search(
        rf"^## {re.escape(title)}\s*$\n(?P<body>.*?)(?=^## |\Z)",
        readme,
        flags=re.MULTILINE | re.DOTALL,
    )
    assert match is not None, f"Missing README section: {title}"
    return re.sub(r"\s+", " ", match.group("body")).strip()


def test_required_publication_files_exist():
    required = ["LICENSE", "CITATION.cff", "SECURITY.md"]

    missing = [name for name in required if not (ROOT / name).is_file()]

    assert missing == []


def test_license_contains_canonical_mit_terms():
    license_text = (ROOT / "LICENSE").read_text()
    normalized = re.sub(r"\s+", " ", license_text).strip()

    assert license_text.startswith("MIT License\n\n")
    assert "Copyright (c) 2026 Andrew Ponomarev" in license_text.splitlines()
    required_terms = [
        (
            "Permission is hereby granted, free of charge, to any person obtaining "
            "a copy of this software and associated documentation files (the "
            '"Software"), to deal in the Software without restriction, including '
            "without limitation the rights to use, copy, modify, merge, publish, "
            "distribute, sublicense, and/or sell copies of the Software, and to "
            "permit persons to whom the Software is furnished to do so, subject "
            "to the following conditions:"
        ),
        (
            "The above copyright notice and this permission notice shall be "
            "included in all copies or substantial portions of the Software."
        ),
        (
            'THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, '
            "EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF "
            "MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. "
            "IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY "
            "CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, "
            "TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE "
            "SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."
        ),
    ]

    for expected in required_terms:
        assert expected in normalized


def test_project_metadata_is_ready_for_publication():
    with (ROOT / "pyproject.toml").open("rb") as file:
        project = tomllib.load(file)["project"]

    assert project["name"] == "clinvar-history"
    assert project["version"] == "0.1.0"
    assert project["readme"] == "README.md"
    assert project["license"] == "MIT"
    assert project["license-files"] == ["LICENSE"]
    assert project["authors"] == [
        {
            "name": "Andrew Ponomarev",
            "email": "19895289+superbsky@users.noreply.github.com",
        }
    ]
    assert {
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    } <= set(project["classifiers"])
    assert project["urls"]["Repository"] == REPOSITORY_URL
    assert project["urls"]["Report"].endswith(
        "/blob/main/docs/full-backfill-results.md"
    )


def test_build_backend_is_pinned_for_reproducible_builds():
    with (ROOT / "pyproject.toml").open("rb") as file:
        build_system = tomllib.load(file)["build-system"]

    assert build_system == {
        "requires": ["hatchling==1.31.0"],
        "build-backend": "hatchling.build",
    }


def test_citation_metadata_identifies_the_release():
    citation = (ROOT / "CITATION.cff").read_text()
    expected_fields = [
        r"^cff-version:\s*1\.2\.0\s*$",
        r'^version:\s*["\']0\.1\.0["\']\s*$',
        r'^date-released:\s*["\']2026-07-25["\']\s*$',
        r"^license:\s*MIT\s*$",
        r"^\s*-\s*family-names:\s*Ponomarev\s*$",
        r"^\s*given-names:\s*Andrew\s*$",
        rf"^repository-code:\s*{re.escape(REPOSITORY_URL)}\s*$",
    ]

    for pattern in expected_fields:
        assert re.search(pattern, citation, flags=re.MULTILINE), pattern


def test_security_policy_states_the_research_scope():
    security = (ROOT / "SECURITY.md").read_text()

    assert "# Security Policy" in security
    assert "not a clinical interpretation tool" in security.lower()


def test_generated_data_and_build_outputs_stay_untracked():
    tracked_pathspecs = [
        "data",
        "artifacts",
        "results",
        "logs",
        "dist",
        ":(glob)**/*.egg-info/**",
    ]
    tracked = subprocess.run(
        ["git", "ls-files", "--", *tracked_pathspecs],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert tracked.stdout.strip() == ""

    examples = [
        "data/example.parquet",
        "artifacts/example.json",
        "results/example.csv",
        "logs/example.log",
        "dist/example.whl",
        "src/clinvar_history.egg-info/PKG-INFO",
    ]
    not_ignored = [
        path
        for path in examples
        if subprocess.run(
            ["git", "check-ignore", "-q", path],
            cwd=ROOT,
            check=False,
        ).returncode
        != 0
    ]

    assert not_ignored == []


def test_ci_is_pinned_and_uses_only_offline_quality_gates():
    workflow_path = ROOT / ".github" / "workflows" / "ci.yml"
    workflow = yaml.load(workflow_path.read_text(), Loader=yaml.BaseLoader)

    assert workflow["on"] == {
        "push": {"branches": ["main"]},
        "pull_request": "",
    }
    assert workflow["permissions"] == {"contents": "read"}
    assert workflow["concurrency"] == {
        "group": "${{ github.workflow }}-${{ github.ref }}",
        "cancel-in-progress": "true",
    }
    assert set(workflow["jobs"]) == {"quality"}

    quality = workflow["jobs"]["quality"]
    assert quality["strategy"] == {
        "fail-fast": "false",
        "matrix": {"python-version": ["3.12", "3.13"]},
    }

    steps = quality["steps"]
    assert [step["name"] for step in steps] == [
        "Check out repository",
        "Install uv",
        "Sync dependencies",
        "Lint",
        "Check formatting",
        "Run tests",
        "Build distribution",
    ]
    assert {
        "runs-on": quality.get("runs-on"),
        "timeout-minutes": quality.get("timeout-minutes"),
        "persist-credentials": steps[0].get("with", {}).get("persist-credentials"),
    } == {
        "runs-on": "ubuntu-24.04",
        "timeout-minutes": "20",
        "persist-credentials": "false",
    }

    action_refs = [step["uses"] for step in steps if "uses" in step]
    assert all(PINNED_ACTION.fullmatch(ref) for ref in action_refs), action_refs
    assert action_refs == [
        "actions/checkout@de0fac2e4500dabe0009e67214ff5f5447ce83dd",
        "astral-sh/setup-uv@08807647e7069bb48b6ef5acd8ec9567f424441b",
    ]
    assert steps[0]["with"] == {"persist-credentials": "false"}
    assert steps[1]["with"] == {
        "version": "0.10.10",
        "python-version": "${{ matrix.python-version }}",
        "enable-cache": "true",
    }

    expected_commands = [
        "uv sync --locked --all-groups",
        "uv run ruff check .",
        "uv run ruff format --check .",
        "uv run pytest -q",
        "uv build",
    ]
    run_commands = [step["run"] for step in steps if "run" in step]
    forbidden = [
        "run_full_backfill.sh",
        "build-manifest",
        "ingest-vcf",
        "ingest-submissions",
    ]
    rejected = [
        command
        for command in forbidden
        if any(command in run_command for run_command in run_commands)
    ]

    assert rejected == [], f"networked CI commands are forbidden: {rejected}"
    assert run_commands == expected_commands
    assert all("if" not in step for step in steps[:-1])
    assert steps[-1]["if"] == "matrix.python-version == '3.13'"


def test_internal_working_documents_stay_unpublished():
    pathspecs = [
        "internal",
        "docs/plans",
        ":(glob)docs/*-prompt.md",
        ":(glob)docs/*-draft.md",
        "docs/contributor-pilot.md",
        "docs/research-tables/liang-crosscheck/chronology.json",
    ]
    tracked = subprocess.run(
        ["git", "ls-files", "--", *pathspecs],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert tracked.stdout.strip() == ""
    examples = [
        "internal/notes.md",
        "docs/plans/example.md",
        "docs/review-prompt.md",
        "docs/manuscript-draft.md",
        "docs/contributor-pilot.md",
        "docs/research-tables/liang-crosscheck/chronology.json",
    ]
    for path in examples:
        assert (
            subprocess.run(
                ["git", "check-ignore", "-q", path], cwd=ROOT, check=False
            ).returncode
            == 0
        ), path


def test_readme_has_publication_header_and_canonical_full_run():
    readme = (ROOT / "README.md").read_text()
    run = _readme_section(readme, "Run")

    required_header = (
        "actions/workflows/ci.yml/badge.svg",
        "https://img.shields.io/badge/Python-3.12%2B",
        "uv sync --locked",
    )
    for expected in required_header:
        assert expected in readme

    assert "scripts/run_full_backfill.sh" in run
    assert run.index("scripts/run_full_backfill.sh") < run.index(
        "analysis/07_build_final_report.py"
    )
    assert "Equivalent manual stages and customization" in run
    assert run.count("--manifest data/full/manifest/releases.parquet") == 2
    assert run.count("--root data/full") == 9
    assert run.count("--end 2026-07-24") == 2
    assert run.count("--study-end 2026-07-24") == 2
    assert "--data-root data/full" in run
    assert "--results-root results/full-v1" in run
    expected_outputs = (
        "attributions.parquet",
        "renovo_quarterly.parquet",
        "sampling_bias.parquet",
        "transition_matrix.parquet",
        "vus_time_to_resolution.parquet",
        "event_mechanisms.parquet",
        "scv_vus_time_to_resolution.parquet",
    )
    for filename in expected_outputs:
        assert f"--output results/full-v1/{filename}" in run
    assert "--attributions results/full-v1/attributions.parquet" in run
    assert "2026-07-31" not in readme
    assert "modern full releases have not yet been benchmarked" not in readme


def test_readme_defines_checkpoint_reuse_as_structural_not_semantic():
    readme = (ROOT / "README.md").read_text()
    run = _readme_section(readme, "Run")

    required = (
        "only structural/readability-based, not semantic validation",
        "Parquet is readable",
        "existing manifest",
        "changed date horizon or configuration",
        "fresh manifest/data/result paths",
        "remove the existing manifest and checkpoints",
    )
    for expected in required:
        assert expected in run


def test_readme_full_backfill_results_section_has_audited_claims():
    readme = (ROOT / "README.md").read_text()
    results = _readme_section(readme, "Full-backfill results")

    required = (
        "417 releases total (338 VCF + 79 SCV), all completed",
        "479,525 canonical VCV transitions",
        "217,177",
        "2,315,022",
        "10.66x",
        "99.91%",
        "80.58%",
        "79.35%",
        "[full English report](docs/full-backfill-results.md)",
    )
    for expected in required:
        assert expected in results


def test_readme_resource_profile_has_audited_storage_and_guidance():
    readme = (ROOT / "README.md").read_text()
    resources = _readme_section(readme, "Resource profile")

    required = (
        "44,250,913,390 bytes (41.21 GiB)",
        "829.5 MiB",
        "seven result Parquet files: 28.6 MiB",
        "raw source-release files retained: zero",
        "8–10 GiB",
        "32 GiB",
        "planning margins, not observed full-run peak-disk or peak-RSS measurements",
    )
    for expected in required:
        assert expected in resources


def test_readme_data_licensing_and_citation_sections_are_scoped():
    readme = (ROOT / "README.md").read_text()
    licensing = _readme_section(readme, "Data source and licensing")
    citation = _readme_section(readme, "Citation")

    licensing_claims = (
        "does not redistribute those releases",
        "https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/",
        "https://www.ncbi.nlm.nih.gov/home/about/policies/",
        "[MIT license](LICENSE)",
        "original code and documentation",
        "third-party submitters",
        "direct diagnosis",
        "medical decision-making",
    )
    for expected in licensing_claims:
        assert expected in licensing
    assert "[CITATION.cff](CITATION.cff)" in citation
    assert "Cite this repository" in citation


def test_readme_pilot_cadence_is_explicitly_superseded():
    readme = (ROOT / "README.md").read_text()
    pilot = _readme_section(readme, "Real-data pilot")

    assert "superseded by the full-backfill cadence analysis above" in pilot
    assert "schema and restart validation" in pilot


def test_markdown_parser_returns_only_rendered_targets(tmp_path):
    document = tmp_path / "links.md"
    document.write_text(
        """```text
[fenced example](fenced.md)
````

``[code example](code.md)``

[reference link][reference]

[reference]: reference.md

[![figure](figure.svg)](report.md)
[balanced path](a(b).md)
"""
    )

    targets = list(_inline_markdown_targets(document))

    assert len(targets) == 4
    assert set(targets) == {"reference.md", "figure.svg", "report.md", "a(b).md"}


def test_local_markdown_link_rejects_absolute_target(tmp_path):
    root = tmp_path / "repository"
    root.mkdir()
    outside = tmp_path / "outside.md"
    outside.write_text("outside")
    document = root / "README.md"
    document.write_text(f"[outside]({outside})\n")

    unresolved = _unresolved_local_markdown_links([document], root)

    assert unresolved == [f"README.md -> {outside}"]


def test_local_markdown_link_rejects_parent_escape(tmp_path):
    root = tmp_path / "repository"
    root.mkdir()
    outside = tmp_path / "outside.md"
    outside.write_text("outside")
    document = root / "README.md"
    document.write_text("[outside](../outside.md)\n")

    unresolved = _unresolved_local_markdown_links([document], root)

    assert unresolved == ["README.md -> ../outside.md"]


def test_local_markdown_link_rejects_symlink_escape(tmp_path):
    root = tmp_path / "repository"
    root.mkdir()
    outside = tmp_path / "outside.md"
    outside.write_text("outside")
    symlink = root / "linked.md"
    symlink.symlink_to(outside)
    document = root / "README.md"
    document.write_text("[outside](linked.md)\n")

    unresolved = _unresolved_local_markdown_links([document], root)

    assert unresolved == ["README.md -> linked.md"]


def test_local_markdown_links_resolve():
    markdown_files = [ROOT / "README.md", *sorted((ROOT / "docs").rglob("*.md"))]
    unresolved = _unresolved_local_markdown_links(markdown_files, ROOT)

    assert unresolved == [], "Unresolved local Markdown links:\n" + "\n".join(
        unresolved
    )
