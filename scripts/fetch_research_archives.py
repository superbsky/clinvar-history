"""Fetch the two dated XML archives used by the bounded historical investigation."""

import shutil
from pathlib import Path

from clinvar_history.research_archives import download_archive


def main():
    destination = Path("artifacts/research-xml")
    destination.mkdir(parents=True, exist_ok=True)
    files = [("2022-08", 2034256048), ("2022-09", 1999297913)]
    remaining = sum(
        size
        for month, size in files
        if not (destination / f"ClinVarVariationRelease_{month}.xml.gz").exists()
    )
    if shutil.disk_usage(destination).free < remaining + 512 * 1024**2:
        raise ValueError(
            "Need remaining compressed bytes plus 512 MiB free before acquisition"
        )
    for month, size in files:
        name = f"ClinVarVariationRelease_{month}.xml.gz"
        print(
            download_archive(
                "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/VCV_xml_old_format/archive/2022/"
                + name,
                destination / name,
                expected_bytes=size,
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
