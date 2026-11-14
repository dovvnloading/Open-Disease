"""Load release-98 human Reactome sets onto a fixed, shared Entrez-gene universe.

Source: Reactome, CC BY 4.0 (https://reactome.org/license). All source files are
pinned release exports. The public products are the consistently defined sets
and an aggregate, outcome-blind count of changes during projection.
"""
from __future__ import annotations

from collections import defaultdict
import csv
from pathlib import Path
import re

from open_disease.geo_microarray import sha256_file

SOURCE_BASE = "https://download.reactome.org/98/"
FILES = {
    "NCBI2Reactome_All_Levels.v98.txt":
        "272ca73c16eb285be1698e07c81912335a1d3804e125f789f5160dff150be243",
    "ReactomePathwaysRelation.v98.txt":
        "aa3609ee852442899b400bf4817b36311d820883ecd75259ed77531d93d8f1dc",
    "ReactomePathways.v98.txt":
        "b8103a50d625ffb6acd18cddffcce0ceb0fdf36b19392d104c32bb1fcc7e2f9f",
}
DISEASE = "R-HSA-1643685"
PRIOR_TARGET = 2354
HUMAN_ID = re.compile(r"R-HSA-([1-9][0-9]*)")
NCBI_ID = re.compile(r"[1-9][0-9]*")


def _rows(path, width):
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE):
            if len(row) != width:
                raise ValueError("Unexpected number of columns in Reactome annotation")
            yield row


def _disease_descendants(names, relations):
    if names.get(DISEASE) != "Disease":
        raise ValueError("Expected human Disease chapter")
    children = defaultdict(set)
    for parent, child in relations:
        if parent in names or child in names:
            if parent not in names or child not in names:
                raise ValueError("Hierarchy edge crosses the human pathway catalogue")
            children[parent].add(child)
    excluded = set()
    pending = [DISEASE]
    while pending:
        parent = pending.pop()
        if parent not in excluded:
            excluded.add(parent)
            pending.extend(children[parent])
    return excluded


def project_memberships(names, members, excluded, shared_genes):
    """Require shared coverage before deterministic deduplication by membership."""
    allowed = set(shared_genes) - {PRIOR_TARGET}
    seen = set()
    sets = {}
    qc = {"excluded_disease_sets": 0, "excluded_size_sets": 0, "excluded_coverage_sets": 0,
          "excluded_duplicate_sets": 0, "tested_pathways": 0}
    for pathway_id in sorted(members, key=lambda key: int(key.rsplit("-", 1)[-1])):
        original = members[pathway_id]
        if pathway_id in excluded:
            qc["excluded_disease_sets"] += 1
            continue
        covered = tuple(sorted(original & allowed))
        if not 15 <= len(covered) <= 200:
            qc["excluded_size_sets"] += 1
            continue
        if 2 * len(covered) < len(original):
            qc["excluded_coverage_sets"] += 1
            continue
        if covered in seen:
            qc["excluded_duplicate_sets"] += 1
            continue
        seen.add(covered)
        sets[pathway_id] = (names[pathway_id], covered)
    qc["tested_pathways"] = len(sets)
    return sets, qc


def load_pathways(directory, shared_genes):
    """Return {stable_id: (name, sorted Entrez IDs)} and aggregate QC for v98."""
    directory = Path(directory)
    for filename, digest in FILES.items():
        if sha256_file(directory / filename) != digest:
            raise ValueError(f"Incorrect Reactome release checksum for {filename}")
    names = {}
    for pathway, name, species in _rows(directory / "ReactomePathways.v98.txt", 3):
        if species != "Homo sapiens":
            continue
        if not HUMAN_ID.fullmatch(pathway) or not name or pathway in names:
            raise ValueError("Invalid or repeated human ID in Reactome catalogue")
        names[pathway] = name
    excluded = _disease_descendants(names, _rows(directory / "ReactomePathwaysRelation.v98.txt", 2))
    members = defaultdict(set)
    human_rows = invalid = 0
    for gene, pathway, url, name, evidence, species in _rows(directory / "NCBI2Reactome_All_Levels.v98.txt", 6):
        if species != "Homo sapiens":
            continue
        human_rows += 1
        if (names.get(pathway) != name or evidence not in {"IEA", "TAS"}
                or url != "https://reactome.org/PathwayBrowser/#/" + pathway):
            raise ValueError("Human gene-pathway mapping differs from the catalogue")
        if not NCBI_ID.fullmatch(gene):
            invalid += 1
            continue
        members[pathway].add(int(gene))
    sets, projection = project_memberships(names, members, excluded, shared_genes)
    return sets, {
        "release": 98, "human_catalogue_pathways": len(names), "human_mapping_rows": human_rows,
        "nonnumeric_mapping_rows": invalid, "human_gene_pathway_pairs": sum(map(len, members.values())),
        "human_gene_bearing_pathways": len(members), "distinct_source_genes": len(set().union(*members.values())),
        "disease_subtree_pathways": len(excluded), **projection,
    }
