import gzip

import pytest

pytest.importorskip("lxml")
from pathlib import Path

from clinvar_history.checkpoints import atomic_json
from clinvar_history.research_archives import extract_archived_records


def test_old_archive_extraction_distinguishes_included_from_interpreted_and_retains_provenance(
    tmp_path,
):
    content = b"""<ClinVarVariationRelease ReleaseDate="2022-08-01"><VariationArchive VariationID="1" Version="2" RecordType="interpreted"><InterpretedRecord><ReviewStatus>criteria provided</ReviewStatus><Interpretations><Interpretation><Description>Uncertain significance</Description></Interpretation></Interpretations><ClinicalAssertionList><ClinicalAssertion SubmissionDate="2022-01-01"><ClinVarAccession Accession="SCV000000001" Version="3"/><ReviewStatus>criteria provided</ReviewStatus><Interpretation DateLastEvaluated="2021-12-01"><Description>Uncertain significance</Description></Interpretation></ClinicalAssertion></ClinicalAssertionList></InterpretedRecord></VariationArchive><VariationArchive VariationID="2" RecordType="included"/></ClinVarVariationRelease>"""
    archive = tmp_path / "source.xml.gz"
    with gzip.open(archive, "wb") as f:
        f.write(content)
    atomic_json({"md5": "fixture"}, archive.with_suffix(".gz.json"))
    result = extract_archived_records(archive, [1, 2, 3], tmp_path / "out")
    assert result["missing"] == [3]
    assert result["scanned_records"] == 2
    one, two = result["records"]
    assert one["aggregate_classifications"] == ["Uncertain significance"]
    assert one["assertions"][0]["accession"] == "SCV000000001"
    assert one["assertions"][0]["version"] == "3"
    assert two["record_type"] == "included"
    assert not two["assertions"]
    assert (Path(tmp_path) / "out/1.xml").exists()
