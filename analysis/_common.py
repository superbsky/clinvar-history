"""Shared offline helpers for command-line analysis scripts."""

import argparse
from dataclasses import asdict, is_dataclass
from enum import Enum
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


def base_parser(description: str, default_output: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--root", type=Path, default=Path("data"))
    parser.add_argument("--output", type=Path, default=Path(default_output))
    return parser


def _plain(value):
    if isinstance(value, Enum):
        return value.value
    if is_dataclass(value):
        return {key: _plain(item) for key, item in asdict(value).items()}
    if isinstance(value, dict):
        return {str(_plain(key)): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    return value


def write_rows(rows, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pylist([_plain(row) for row in rows])
    pq.write_table(table, output, compression="zstd")
