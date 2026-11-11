import httpx
import pytest

from clinvar_history.xml_review import fetch_xml, parse_vcv_xml

XML = b"""<ClinVarResult-Set><VariationArchive VariationID="42" Accession="VCV000000042" Version="2">
<ClassifiedRecord><Classifications><GermlineClassification>
<Description>Uncertain significance</Description><Explanation>Two disagreeing submissions</Explanation>
</GermlineClassification></Classifications><ClinicalAssertionList>
<ClinicalAssertion ContributesToAggregateClassification="false">
<ClinVarAccession Accession="SCV000000001" Version="3" SubmitterName="Lab"/>
<Classification DateLastEvaluated="2020-01-01"><ReviewStatus>criteria provided</ReviewStatus>
<GermlineClassification>Benign</GermlineClassification></Classification>
<Comment>Evidence details</Comment><AttributeSet><Attribute Type="AssertionMethod">ACMG</Attribute></AttributeSet>
</ClinicalAssertion><ClinicalAssertion><ClinVarAccession Accession="SCV000000002"/>
<Classification><OncogenicityClassification>Oncogenic</OncogenicityClassification></Classification>
</ClinicalAssertion></ClinicalAssertionList></ClassifiedRecord></VariationArchive></ClinVarResult-Set>"""


def test_xml_preserves_richer_evidence_and_separates_germline():
    result = parse_vcv_xml(XML, 42)
    assert result["current_germline_classification"] == "Uncertain significance"
    assert result["aggregate_explanations"] == ["Two disagreeing submissions"]
    assert result["assertions"][0]["contributes"] == "false"
    assert result["assertions"][0]["assertion_methods"] == ["ACMG"]
    assert result["assertions"][1]["classification"] is None
    assert result["assertions"][1]["classification_type"] == "other"


def test_xml_rejects_wrong_variation_and_entities():
    with pytest.raises(ValueError, match="VariationID"):
        parse_vcv_xml(XML, 43)
    with pytest.raises(ValueError, match="DTD"):
        parse_vcv_xml(b"<!DOCTYPE x><x/>", 42)


def test_fetch_rejects_response_beyond_record_budget():
    def handler(request):
        assert request.url.params["is_variationid"] == "true"
        return httpx.Response(200, content=XML)

    with (
        httpx.Client(transport=httpx.MockTransport(handler)) as client,
        pytest.raises(ValueError, match="budget"),
    ):
        fetch_xml(42, client=client, max_bytes=10)


def test_version_pair_requires_exact_versions_and_reuses_verified_cache(
    tmp_path, monkeypatch
):
    import xml.etree.ElementTree as ET

    import clinvar_history.xml_review as review

    old = ET.fromstring(XML)
    old.find("VariationArchive").set("Version", "1")
    assertions = old.find(".//ClinicalAssertionList")
    assertions.remove(assertions[-1])
    old_bytes = ET.tostring(old)
    calls = []

    def handler(request):
        calls.append(request.url.params["id"])
        assert "is_variationid" not in request.url.params
        return httpx.Response(
            200, content=old_bytes if request.url.params["id"].endswith(".1") else XML
        )

    real_client = httpx.Client
    monkeypatch.setattr(
        review.httpx,
        "Client",
        lambda **kwargs: real_client(transport=httpx.MockTransport(handler)),
    )
    monkeypatch.setattr(review.time, "sleep", lambda _: None)
    result = review.compare_versions(42, 1, 2, tmp_path)
    assert calls == ["VCV000000042.1", "VCV000000042.2"]
    assert any(
        c["change"] == "added" and c["accession"] == "SCV000000002"
        for c in result["scv_changes"]
    )
    assert review.compare_versions(42, 1, 2, tmp_path, offline=True) == result
    assert len(calls) == 2
    (tmp_path / "VCV000000042.1.xml").write_bytes(b"bad")
    with pytest.raises(ValueError, match="corrupt"):
        review.compare_versions(42, 1, 2, tmp_path, offline=True)


def test_version_pair_rejects_api_fallback_to_current_version(tmp_path, monkeypatch):
    import clinvar_history.xml_review as review

    real_client = httpx.Client
    monkeypatch.setattr(
        review.httpx,
        "Client",
        lambda **kwargs: real_client(
            transport=httpx.MockTransport(lambda _: httpx.Response(200, content=XML))
        ),
    )
    with pytest.raises(ValueError, match="different VCV version"):
        review.compare_versions(42, 1, 2, tmp_path)
