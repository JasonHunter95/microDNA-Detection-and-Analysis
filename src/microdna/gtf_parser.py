#!/usr/bin/env python3
"""
gtf_parser.py - Shared utilities for parsing GTF attribute strings.

Provides functions to extract gene metadata from GTF-style attribute strings,
commonly found in BED files derived from GENCODE annotations via gtf2bed.
"""

import re
from typing import Any


def extract_gene_info(attribute_field: str | None) -> dict[str, Any]:
    """
    Parse a GTF-style attribute string to extract gene metadata.

    GTF attribute format: key "value"; key "value"; ...
    Example: gene_id "ENSG00000223972.5"; gene_name "DDX11L1"; gene_type "transcribed_unprocessed_pseudogene";

    Args:
        attribute_field: The raw attribute string from GTF/BED file.
            Can be None or NaN from pandas.

    Returns:
        Dictionary containing:
        - gene_id: ENSEMBL gene ID (e.g., "ENSG00000223972.5")
        - gene_name: Gene symbol (e.g., "DDX11L1")
        - gene_type: Gene biotype (e.g., "protein_coding", "lncRNA")
        - transcript_id: ENSEMBL transcript ID
        - transcript_name: Transcript name
        - level: Annotation confidence level (1=verified, 2=manually annotated, 3=automated)
        - gene_status: Gene status (e.g., "KNOWN", "NOVEL")

        Missing fields are set to None.
    """
    # Handle None or NaN values
    if attribute_field is None:
        attribute_field = ""
    else:
        attribute_field = str(attribute_field)

    # Define patterns for each field
    patterns = {
        "gene_id": r'gene_id "([^"]+)"',
        "gene_name": r'gene_name "([^"]+)"',
        "gene_type": r'gene_type "([^"]+)"',
        "transcript_id": r'transcript_id "([^"]+)"',
        "transcript_name": r'transcript_name "([^"]+)"',
        "level": r"level (\d+)",
        "gene_status": r'gene_status "([^"]+)"',
    }

    result: dict[str, Any] = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, attribute_field)
        if match:
            value = match.group(1)
            # Convert level to int
            result[key] = int(value) if key == "level" else value
        else:
            result[key] = None

    return result


def extract_basic_gene_info(attribute_field: str | None) -> dict[str, str]:
    """
    Parse a GTF attribute string to extract basic gene identifiers.

    Simplified version that returns only gene_id, gene_name, and gene_type,
    with "." as default for missing values.

    Args:
        attribute_field: The raw attribute string from GTF/BED file.

    Returns:
        Dictionary with gene_id, gene_name, gene_type (defaults to ".").
    """
    if not isinstance(attribute_field, str):
        return {"gene_id": ".", "gene_name": ".", "gene_type": "."}

    gene_id = re.search(r'gene_id "([^"]+)"', attribute_field)
    gene_name = re.search(r'gene_name "([^"]+)"', attribute_field)
    gene_type = re.search(r'gene_type "([^"]+)"', attribute_field)

    return {
        "gene_id": gene_id.group(1) if gene_id else ".",
        "gene_name": gene_name.group(1) if gene_name else ".",
        "gene_type": gene_type.group(1) if gene_type else ".",
    }
