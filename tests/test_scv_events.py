from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import replace
from datetime import date
from pathlib import Path

import pyarrow.parquet as pq
from click.testing import CliRunner

from clinvar_history.classification import ClassificationState
from clinvar_history.cli import main
from clinvar_history.manifest import ReleaseFile, ReleaseKind
from clinvar_history.pipeline import IngestionSummary, ingest_submission_releases
from clinvar_history.scv_diff import (
    ScvEventType,
    consolidate_scv_records,
    diff_scv_states,
)
from clinvar_history.storage import read_scv_state
from clinvar_history.submission_parser import parse_submissions

FIXTURES = Path(__file__).parent / "fixtures"


def _records(name):
    with (FIXTURES / name).open() as stream:
        return list(parse_submissions(stream))


def test_scv_diff_separates_reclassification_from_other_updates():
    events = list(
        diff_scv_states(
            _records("submission_a.txt"),
            _records("submission_b.txt"),
            observed_at=date(2024, 2, 1),
        )
    )
    by_accession = {}
    for event in events:
        by_accession.setdefault(event.scv_accession, []).append(event.event_type)

    assert by_accession["SCV000000001"] == [
        ScvEventType.SCV_VERSION_CHANGED,
        ScvEventType.SCV_RECLASSIFIED,
    ]
    assert by_accession["SCV000000002"] == [ScvEventType.SCV_VERSION_CHANGED]
    assert by_accession["SCV000000003"] == [
        ScvEventType.SCV_VERSION_CHANGED,
        ScvEventType.REVIEW_STATUS_CHANGED,
        ScvEventType.CONDITION_SET_CHANGED,
    ]
    assert by_accession["SCV000000004"] == [ScvEventType.SCV_REMOVED]
    assert by_accession["SCV000000005"] == [ScvEventType.SCV_ADDED]
    assert ScvEventType.SCV_RECLASSIFIED not in by_accession["SCV000000002"]


def test_scv_events_preserve_old_and_new_values():
    events = list(
        diff_scv_states(
            _records("submission_a.txt"),
            _records("submission_b.txt"),
            observed_at=date(2024, 2, 1),
        )
    )
    reclassified = next(
        event for event in events if event.event_type is ScvEventType.SCV_RECLASSIFIED
    )

    assert reclassified.variation_id == 123
    assert reclassified.old_version == 1
    assert reclassified.new_version == 2
    assert reclassified.old_raw_classification == "Uncertain significance"
    assert reclassified.new_raw_classification == "Likely pathogenic"
    assert reclassified.old_submitter == "Lab A"
    assert reclassified.new_submitter == "Lab A"
    assert reclassified.old_collection_method == "clinical testing"
    assert reclassified.new_collection_method == "clinical testing"


def test_duplicate_scv_prefers_classification_over_historical_missing_row():
    pathogenic = _records("submission_a.txt")[-1]
    missing = replace(
        pathogenic,
        raw_classification="not provided",
        classification=ClassificationState.MISSING,
        classification_ambiguous=False,
        last_evaluated=date(2023, 8, 1),
    )
    older_missing = replace(missing, last_evaluated=date(2023, 7, 1))

    events = list(
        diff_scv_states(
            [],
            [missing, older_missing, pathogenic],
            observed_at=date(2024, 1, 1),
        )
    )

    assert len(events) == 1
    assert events[0].event_type is ScvEventType.SCV_ADDED
    assert events[0].new_classification is ClassificationState.PATHOGENIC
    assert consolidate_scv_records([missing, older_missing, pathogenic]) == [pathogenic]


def test_duplicate_scv_with_conflicting_classes_prefers_latest_evaluation():
    pathogenic = _records("submission_a.txt")[-1]
    older_pathogenic = replace(
        pathogenic,
        last_evaluated=date(2014, 7, 1),
    )
    later_vus = replace(
        pathogenic,
        raw_classification="Uncertain significance",
        classification=ClassificationState.VUS,
        last_evaluated=date(2019, 3, 25),
    )

    assert consolidate_scv_records([older_pathogenic, later_vus]) == [later_vus]
    assert consolidate_scv_records([later_vus, older_pathogenic]) == [later_vus]


def test_duplicate_scv_with_tied_conflicting_classes_becomes_conflict():
    pathogenic = _records("submission_a.txt")[-1]
    vus = replace(
        pathogenic,
        raw_classification="Uncertain significance",
        classification=ClassificationState.VUS,
    )

    forward = consolidate_scv_records([pathogenic, vus])
    reverse = consolidate_scv_records([vus, pathogenic])

    assert forward == reverse
    assert forward[0].classification is ClassificationState.CONFLICT
    assert "Pathogenic" in forward[0].raw_classification
    assert "Uncertain significance" in forward[0].raw_classification


@contextmanager
def _fixture_source(release: ReleaseFile) -> Iterator[Iterator[bytes]]:
    fixture = release.url.removeprefix("fixture://")
    with (FIXTURES / fixture).open("rb") as stream:
        yield iter(stream)


def test_submission_ingestion_writes_compact_state_and_partitioned_events(tmp_path):
    releases = [
        ReleaseFile(
            release_date=date(2024, 1, 1),
            kind=ReleaseKind.SUBMISSION_SUMMARY,
            url="fixture://submission_a.txt",
            size_bytes=1,
        ),
        ReleaseFile(
            release_date=date(2024, 2, 1),
            kind=ReleaseKind.SUBMISSION_SUMMARY,
            url="fixture://submission_b.txt",
            size_bytes=1,
        ),
    ]

    first = ingest_submission_releases(
        releases,
        root=tmp_path,
        source=_fixture_source,
    )
    second = ingest_submission_releases(
        releases,
        root=tmp_path,
        source=_fixture_source,
    )

    assert first == IngestionSummary(processed=2, skipped=0)
    assert second == IngestionSummary(processed=0, skipped=2)
    assert len(read_scv_state(tmp_path / "state" / "scv_latest.parquet")) == 4
    event_files = sorted((tmp_path / "events" / "scv").rglob("*.parquet"))
    assert len(event_files) == 2
    assert sum(pq.read_table(path).num_rows for path in event_files) == 12


def test_ingest_submissions_cli_is_available():
    result = CliRunner().invoke(main, ["ingest-submissions", "--help"])

    assert result.exit_code == 0
    assert "--manifest" in result.output
    assert "--root" in result.output
    assert "--keep-raw" in result.output
