import json
import os
import subprocess
import sys
from datetime import date
from pathlib import Path
from xml.etree import ElementTree

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from clinvar_history.reporting import (
    ReportValidationError,
    audit_full_backfill,
    build_full_report_summary,
    render_report_figures,
    summarize_mechanisms,
    summarize_outcomes,
    summarize_quarters,
    summarize_sampling,
    summarize_transitions,
)

EXPECTED_FIGURES = {
    "vus-quarterly-trend.svg",
    "temporal-resolution.svg",
    "transition-matrix.svg",
    "candidate-mechanisms.svg",
    "vus-episode-outcomes.svg",
}


def _quarterly_table() -> pa.Table:
    return pa.Table.from_pylist(
        [
            {
                "period_start": date(2020, 1, 1),
                "period_end": date(2020, 4, 1),
                "vus_previous": 100,
                "vus_current": 120,
                "reclassified_vus": 2,
                "net_vus_change": 20,
                "reclassification_rate": 2 / 120,
                "discovery_rate": 20 / 120,
                "discovery_reclassification_ratio": 10.0,
            },
            {
                "period_start": date(2020, 4, 1),
                "period_end": date(2020, 7, 1),
                "vus_previous": 120,
                "vus_current": 200,
                "reclassified_vus": 4,
                "net_vus_change": 80,
                "reclassification_rate": 4 / 200,
                "discovery_rate": 80 / 200,
                "discovery_reclassification_ratio": 20.0,
            },
            {
                "period_start": date(2020, 7, 1),
                "period_end": date(2020, 7, 20),
                "vus_previous": 200,
                "vus_current": 250,
                "reclassified_vus": 3,
                "net_vus_change": 50,
                "reclassification_rate": 3 / 250,
                "discovery_rate": 50 / 250,
                "discovery_reclassification_ratio": 50 / 3,
            },
        ]
    )


def test_summarize_quarters_separates_complete_and_partial_intervals():
    summary = summarize_quarters(_quarterly_table(), split_date=date(2020, 4, 1))

    assert summary.initial_vus == 100
    assert summary.final_vus == 250
    assert summary.net_vus_growth == 150
    assert summary.fold_change == pytest.approx(2.5)
    assert summary.total_reclassified == 9
    assert summary.complete_quarters == 2
    assert summary.partial_intervals == 1
    assert summary.complete.mean_reclassification == pytest.approx(
        (2 / 120 + 4 / 200) / 2
    )
    assert summary.before_split.quarters == 1
    assert summary.after_split.quarters == 1


def test_summarize_quarters_rejects_inconsistent_inventory_chain():
    table = _quarterly_table().set_column(
        2,
        "vus_previous",
        pa.array([100, 999, 200]),
    )

    with pytest.raises(ReportValidationError, match="inventory chain"):
        summarize_quarters(table, split_date=date(2020, 4, 1))


def test_summarize_sampling_recomputes_recall_and_missed_counts():
    table = pa.Table.from_pylist(
        [
            {
                "cadence": "weekly",
                "true_transitions": 100,
                "observed_transitions": 99,
                "event_recall": 0.99,
                "missed_recurrent_transitions": 1,
                "missed_conflict_episodes": 0,
                "mean_event_date_error_days": 3.5,
            },
            {
                "cadence": "monthly",
                "true_transitions": 100,
                "observed_transitions": 80,
                "event_recall": 0.8,
                "missed_recurrent_transitions": 20,
                "missed_conflict_episodes": 7,
                "mean_event_date_error_days": 14.0,
            },
        ]
    )

    rows = summarize_sampling(table)

    assert rows["monthly"].recall == pytest.approx(0.8)
    assert rows["monthly"].missed_transitions == 20
    assert rows["monthly"].missed_conflict_episodes == 7


def test_summarize_sampling_rejects_declared_recall_mismatch():
    table = pa.Table.from_pylist(
        [
            {
                "cadence": "monthly",
                "true_transitions": 100,
                "observed_transitions": 80,
                "event_recall": 0.9,
                "missed_recurrent_transitions": 20,
                "missed_conflict_episodes": 7,
                "mean_event_date_error_days": 14.0,
            }
        ]
    )

    with pytest.raises(ReportValidationError, match="recall"):
        summarize_sampling(table)


def test_summarize_mechanisms_sorts_and_recomputes_shares():
    table = pa.Table.from_pylist(
        [
            {"primary_mechanism": "unexplained", "count": 3},
            {"primary_mechanism": "scv_added", "count": 7},
        ]
    )

    summary = summarize_mechanisms(table)

    assert [row.mechanism for row in summary] == ["scv_added", "unexplained"]
    assert summary[0].share == pytest.approx(0.7)
    assert summary[1].share == pytest.approx(0.3)


def test_summarize_outcomes_reports_counts_shares_and_observed_medians():
    table = pa.Table.from_pylist(
        [
            {"duration_days": 100, "outcome": "censored"},
            {"duration_days": 200, "outcome": "censored"},
            {"duration_days": 10, "outcome": "benign"},
            {"duration_days": 20, "outcome": "benign"},
            {"duration_days": 30, "outcome": "pathogenic"},
        ]
    )

    summary = summarize_outcomes(table)

    assert summary.total == 5
    assert summary.counts["censored"] == 2
    assert summary.shares["benign"] == pytest.approx(0.4)
    assert summary.median_days["benign"] == 15
    assert summary.median_days["pathogenic"] == 30


def test_summarize_outcomes_rejects_unknown_outcome():
    table = pa.Table.from_pylist([{"duration_days": 10, "outcome": "unknown"}])

    with pytest.raises(ReportValidationError, match="outcome"):
        summarize_outcomes(table)


def test_summarize_outcomes_accepts_source_other_exit_value():
    table = pa.Table.from_pylist([{"duration_days": 10, "outcome": "other_exit"}])

    summary = summarize_outcomes(table)

    assert summary.counts["other_exit"] == 1


def test_summarize_transitions_preserves_direction_and_total():
    table = pa.Table.from_pylist(
        [
            {"old_classification": "vus", "new_classification": "benign", "count": 4},
            {"old_classification": "benign", "new_classification": "vus", "count": 1},
        ]
    )

    summary = summarize_transitions(table)

    assert summary.total == 5
    assert summary.counts[("vus", "benign")] == 4
    assert summary.counts[("benign", "vus")] == 1


def _write_table(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows), path)


def _write_full_report_fixture(tmp_path: Path) -> tuple[Path, Path]:
    data_root = tmp_path / "data"
    results_root = tmp_path / "results"
    manifest_rows = [
        {
            "release_date": date(2020, 1, 1),
            "source_month": None,
            "kind": "vcf",
            "url": "https://example.test/a.vcf.gz",
            "size_bytes": 100,
            "md5": "a" * 32,
        },
        {
            "release_date": date(2020, 2, 1),
            "source_month": "2020-02",
            "kind": "submission_summary",
            "url": "https://example.test/b.txt.gz",
            "size_bytes": 200,
            "md5": None,
        },
    ]
    _write_table(data_root / "manifest" / "releases.parquet", manifest_rows)

    connection = duckdb.connect(str(data_root / "clinvar.duckdb"))
    connection.execute(
        """
        CREATE TABLE pipeline_runs (
            kind VARCHAR,
            url VARCHAR,
            release_date DATE,
            status VARCHAR,
            source_md5 VARCHAR,
            official_md5 VARCHAR
        )
        """
    )
    connection.executemany(
        "INSERT INTO pipeline_runs VALUES (?, ?, ?, ?, ?, ?)",
        [
            (
                "vcf",
                "https://example.test/a.vcf.gz",
                date(2020, 1, 1),
                "completed",
                "a" * 32,
                "a" * 32,
            ),
            (
                "submission_summary",
                "https://example.test/b.txt.gz",
                date(2020, 2, 1),
                "completed",
                "b" * 32,
                None,
            ),
        ],
    )
    connection.close()

    _write_table(
        data_root / "state" / "vcv_latest.parquet",
        [{"variation_id": 1, "classification": "vus"}],
    )
    _write_table(
        data_root / "state" / "scv_latest.parquet",
        [{"scv_accession": "SCV1", "classification": "vus"}],
    )
    _write_table(
        data_root / "events" / "vcv" / "release_date=2020-01-01" / "events.parquet",
        [
            {
                "variation_id": variation_id,
                "event_type": "classification_changed",
                "observed_at": date(2020, 1, 1),
            }
            for variation_id in range(1, 4)
        ],
    )
    _write_table(
        data_root / "events" / "scv" / "release_date=2020-02-01" / "events.parquet",
        [
            {
                "scv_accession": "SCV1",
                "event_type": "scv_added",
                "observed_at": date(2020, 2, 1),
            }
        ],
    )

    _write_table(
        results_root / "attributions.parquet",
        [
            {
                "variation_id": variation_id,
                "observed_at": date(2020, 1, 1),
                "primary_mechanism": "scv_added",
                "candidate_mechanisms": ["scv_added"],
                "supporting_event_ids": [f"SCV{variation_id}"],
            }
            for variation_id in range(1, 4)
        ],
    )
    _write_table(
        results_root / "event_mechanisms.parquet",
        [{"primary_mechanism": "scv_added", "count": 3}],
    )
    pq.write_table(_quarterly_table(), results_root / "renovo_quarterly.parquet")
    _write_table(
        results_root / "sampling_bias.parquet",
        [
            {
                "cadence": cadence,
                "true_transitions": 3,
                "observed_transitions": observed,
                "event_recall": observed / 3,
                "missed_recurrent_transitions": 3 - observed,
                "missed_conflict_episodes": 0,
                "mean_event_date_error_days": error,
            }
            for cadence, observed, error in (
                ("weekly", 3, 3.0),
                ("monthly", 2, 14.0),
                ("quarterly", 1, 45.0),
            )
        ],
    )
    _write_table(
        results_root / "transition_matrix.parquet",
        [
            {"old_classification": "vus", "new_classification": "benign", "count": 2},
            {
                "old_classification": "vus",
                "new_classification": "pathogenic",
                "count": 1,
            },
        ],
    )
    _write_table(
        results_root / "vus_time_to_resolution.parquet",
        [
            {
                "variation_id": 1,
                "started_at": date(2020, 1, 1),
                "ended_at": date(2020, 1, 11),
                "duration_days": 10,
                "outcome": "benign",
                "review_status": "single submitter",
                "gene_info": "GENE:1",
                "conditions": "condition",
                "variant_type": "SNV",
            }
        ],
    )
    _write_table(
        results_root / "scv_vus_time_to_resolution.parquet",
        [
            {
                "scv_accession": "SCV1",
                "variation_id": 1,
                "started_at": date(2020, 1, 1),
                "ended_at": date(2020, 1, 11),
                "duration_days": 10,
                "outcome": "pathogenic",
                "submitter": "lab",
                "collection_method": "clinical testing",
                "review_status": "single submitter",
                "conditions": "condition",
            }
        ],
    )
    return data_root, results_root


def test_audit_accepts_consistent_artifacts(tmp_path):
    data_root, results_root = _write_full_report_fixture(tmp_path)

    audit = audit_full_backfill(
        data_root,
        results_root,
        expected_vcf_releases=1,
        expected_scv_releases=1,
    )

    assert audit.ok
    assert audit.vcf_releases == 1
    assert audit.scv_releases == 1
    assert audit.completed_ledger_rows == 2
    assert audit.available_official_checksums == 1
    assert audit.cross_table_transition_total == 3


def test_audit_rejects_manifest_ledger_mismatch(tmp_path):
    data_root, results_root = _write_full_report_fixture(tmp_path)
    connection = duckdb.connect(str(data_root / "clinvar.duckdb"))
    connection.execute("DELETE FROM pipeline_runs WHERE kind = 'submission_summary'")
    connection.close()

    with pytest.raises(ReportValidationError, match="manifest/ledger"):
        audit_full_backfill(data_root, results_root)


def test_audit_rejects_checksum_mismatch(tmp_path):
    data_root, results_root = _write_full_report_fixture(tmp_path)
    connection = duckdb.connect(str(data_root / "clinvar.duckdb"))
    connection.execute(
        "UPDATE pipeline_runs SET source_md5 = ? WHERE kind = 'vcf'",
        ["c" * 32],
    )
    connection.close()

    with pytest.raises(ReportValidationError, match="checksum"):
        audit_full_backfill(data_root, results_root)


def test_audit_rejects_cross_table_total_mismatch(tmp_path):
    data_root, results_root = _write_full_report_fixture(tmp_path)
    _write_table(
        results_root / "event_mechanisms.parquet",
        [{"primary_mechanism": "scv_added", "count": 2}],
    )

    with pytest.raises(ReportValidationError, match="transition totals"):
        audit_full_backfill(data_root, results_root)


def test_audit_rejects_invalid_episode_duration(tmp_path):
    data_root, results_root = _write_full_report_fixture(tmp_path)
    table = pq.read_table(results_root / "vus_time_to_resolution.parquet")
    duration_index = table.schema.get_field_index("duration_days")
    table = table.set_column(duration_index, "duration_days", pa.array([9]))
    pq.write_table(table, results_root / "vus_time_to_resolution.parquet")

    with pytest.raises(ReportValidationError, match="duration"):
        audit_full_backfill(data_root, results_root)


def test_audit_counts_hard_linked_files_once(tmp_path):
    data_root, results_root = _write_full_report_fixture(tmp_path)
    before = audit_full_backfill(data_root, results_root).retained_data_bytes
    source = data_root / "state" / "vcv_latest.parquet"
    duplicate = data_root / "state" / "vcv_checkpoint.parquet"
    os.link(source, duplicate)

    after = audit_full_backfill(data_root, results_root).retained_data_bytes

    assert after == before


def test_build_full_report_summary_reads_every_result_layer(tmp_path):
    data_root, results_root = _write_full_report_fixture(tmp_path)

    summary = build_full_report_summary(
        data_root,
        results_root,
        expected_vcf_releases=1,
        expected_scv_releases=1,
    )

    assert summary.audit.cross_table_transition_total == 3
    assert summary.quarterly.final_vus == 250
    assert summary.sampling["monthly"].observed_transitions == 2
    assert summary.mechanisms[0].mechanism == "scv_added"
    assert summary.transitions.total == 3
    assert summary.vcv_outcomes.counts["benign"] == 1
    assert summary.scv_outcomes.counts["pathogenic"] == 1


def test_render_report_figures_are_accessible_and_deterministic(tmp_path):
    data_root, results_root = _write_full_report_fixture(tmp_path)
    summary = build_full_report_summary(data_root, results_root)
    figures_dir = tmp_path / "figures"

    first_paths = render_report_figures(summary, figures_dir)
    first_bytes = {name: path.read_bytes() for name, path in first_paths.items()}
    second_paths = render_report_figures(summary, figures_dir)

    assert set(first_paths) == EXPECTED_FIGURES
    assert set(second_paths) == EXPECTED_FIGURES
    for name, path in second_paths.items():
        content = path.read_bytes()
        assert content == first_bytes[name]
        root = ElementTree.fromstring(content)
        namespace = "{http://www.w3.org/2000/svg}"
        assert root.tag == f"{namespace}svg"
        assert root.find(f"{namespace}title") is not None
        assert root.find(f"{namespace}desc") is not None
        ids = [
            element.attrib["id"] for element in root.iter() if "id" in element.attrib
        ]
        assert len(ids) == len(set(ids))
        labelled_by = root.attrib["aria-labelledby"].split()
        assert set(labelled_by).issubset(ids)


def test_svg_base_style_does_not_override_direct_text_colors(tmp_path):
    data_root, results_root = _write_full_report_fixture(tmp_path)
    summary = build_full_report_summary(data_root, results_root)
    path = render_report_figures(summary, tmp_path / "figures")[
        "temporal-resolution.svg"
    ]
    root = ElementTree.parse(path).getroot()
    style = root.find("{http://www.w3.org/2000/svg}style")

    assert style is not None
    base_text_rule = style.text.split("}", 1)[0]
    assert "fill:" not in base_text_rule


def test_transition_matrix_column_labels_are_painted_after_cells(tmp_path):
    data_root, results_root = _write_full_report_fixture(tmp_path)
    summary = build_full_report_summary(data_root, results_root)
    path = render_report_figures(summary, tmp_path / "figures")["transition-matrix.svg"]
    elements = list(ElementTree.parse(path).getroot())
    cell_indexes = [
        index
        for index, element in enumerate(elements)
        if element.attrib.get("class") == "matrix-cell"
    ]
    label_indexes = [
        index
        for index, element in enumerate(elements)
        if element.attrib.get("class") == "matrix-column-label"
    ]

    assert len(cell_indexes) == 64
    assert len(label_indexes) == 8
    assert min(label_indexes) > max(cell_indexes)


def test_final_report_cli_audits_and_writes_figures(tmp_path):
    data_root, results_root = _write_full_report_fixture(tmp_path)
    figures_dir = tmp_path / "cli-figures"
    project_root = Path(__file__).parents[1]

    completed = subprocess.run(
        [
            sys.executable,
            str(project_root / "analysis" / "07_build_final_report.py"),
            "--data-root",
            str(data_root),
            "--results-root",
            str(results_root),
            "--figures-dir",
            str(figures_dir),
            "--expected-vcf-releases",
            "1",
            "--expected-scv-releases",
            "1",
        ],
        cwd=project_root,
        env={**os.environ, "PYTHONPATH": str(project_root / "src")},
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    payload = json.loads(completed.stdout)
    assert payload["audit"]["completed_ledger_rows"] == 2
    assert payload["headline"]["cross_table_transition_total"] == 3
    assert set(payload["figures"]) == EXPECTED_FIGURES
    assert {path.name for path in figures_dir.glob("*.svg")} == EXPECTED_FIGURES
