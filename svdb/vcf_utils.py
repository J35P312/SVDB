"""Shared VCF file utilities used across SVDB modules."""
import gzip
from typing import Dict, IO, Tuple


def open_vcf(path: str) -> IO[str]:
    """Open a plain or gzip-compressed VCF file for reading.

    Returns a text-mode file handle suitable for use as a context manager::

        with open_vcf(path) as lines:
            for line in lines:
                ...
    """
    opener = gzip.open if path.endswith(".gz") else open
    return opener(path, "rt")


def normalize_chrom(name: str) -> str:
    """Strip leading 'chr'/'Chr'/'CHR' prefix from a chromosome name."""
    return name.replace("chr", "").replace("Chr", "").replace("CHR", "")


def parse_info_field(info_str: str) -> Dict[str, str]:
    """Parse a VCF INFO field string into a key→value dict.

    FLAG entries (no '=') are silently skipped — they carry no value.
    Uses maxsplit=1 so values that contain '=' are preserved intact.
    """
    result: Dict[str, str] = {}
    for tag in info_str.split(";"):
        parts = tag.split("=", 1)
        if len(parts) > 1:
            result[parts[0]] = parts[1]
    return result


def parse_ci(info_value: str) -> Tuple[int, int]:
    """Parse a CIPOS or CIEND INFO value into an (lower, upper) interval pair.

    Handles both single-value ('500') and two-value ('-100,200') forms,
    with or without surrounding parentheses. Returns absolute integer values.

    Examples::

        parse_ci("100,200")   -> (100, 200)
        parse_ci("(-50,50)")  -> (50, 50)
        parse_ci("300")       -> (300, 300)
    """
    parts = info_value.replace("(", "").replace(")", "").split(",")
    lower = abs(int(parts[0]))
    upper = abs(int(parts[1])) if len(parts) > 1 else lower
    return lower, upper
