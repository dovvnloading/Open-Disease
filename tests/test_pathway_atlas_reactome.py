"""Outcome-independent Reactome projections, including overlap and cutoff edges."""
import pytest

from scripts.pathway_atlas_reactome import DISEASE, PRIOR_TARGET, _disease_descendants, project_memberships


def test_disease_subtree_includes_descendants_with_a_second_parent():
    names = {DISEASE: "Disease", "R-HSA-2": "Other", "R-HSA-3": "D1", "R-HSA-4": "D2"}
    edges = [(DISEASE, "R-HSA-3"), ("R-HSA-3", "R-HSA-4"), ("R-HSA-2", "R-HSA-4")]
    assert _disease_descendants(names, edges) == {DISEASE, "R-HSA-3", "R-HSA-4"}
    with pytest.raises(ValueError, match="crosses"):
        _disease_descendants(names, [(DISEASE, "R-MMU-1")])


def test_projection_uses_full_denominator_and_smallest_numeric_id():
    universe = set(range(1, 202)) | {PRIOR_TARGET}
    fifteen = set(range(1, 16))
    outside = set(range(301, 316))
    members = {
        "R-HSA-10": fifteen, "R-HSA-2": fifteen | outside,
        "R-HSA-4": fifteen | outside | {PRIOR_TARGET}, "R-HSA-5": set(range(1, 15)) | {PRIOR_TARGET},
        "R-HSA-6": set(range(1, 202)), "R-HSA-7": set(range(1, 201)),
        DISEASE: fifteen,
    }
    names = {k: k for k in members}
    sets, qc = project_memberships(names, members, {DISEASE}, universe)
    assert list(sets) == ["R-HSA-2", "R-HSA-7"]
    assert sets["R-HSA-2"][1] == tuple(range(1, 16))
    assert qc == {"excluded_disease_sets": 1, "excluded_size_sets": 2, "excluded_coverage_sets": 1,
                  "excluded_duplicate_sets": 1, "tested_pathways": 2}


def test_pinned_parser_discards_nonnumeric_ids_and_collapses_evidence(tmp_path, monkeypatch):
    from hashlib import sha256
    from scripts import pathway_atlas_reactome as reactome

    stable = "R-HSA-25"
    entries = []
    for gene in [*range(1, 16), 2, "NM_1"]:
        evidence = "TAS" if len(entries) < 15 else "IEA"
        entries.append(f"{gene}\t{stable}\thttps://reactome.org/PathwayBrowser/#/{stable}"
                       f"\tLeaf\t{evidence}\tHomo sapiens")
    sources = {
        "ReactomePathways.v98.txt": f"{DISEASE}\tDisease\tHomo sapiens\n{stable}\tLeaf\tHomo sapiens\n",
        "ReactomePathwaysRelation.v98.txt": "",
        "NCBI2Reactome_All_Levels.v98.txt": "\n".join(entries) + "\n",
    }
    hashes = {}
    for filename, contents in sources.items():
        path = tmp_path / filename
        path.write_bytes(contents.encode())
        hashes[filename] = sha256(path.read_bytes()).hexdigest()
    monkeypatch.setattr(reactome, "FILES", hashes)
    sets, qc = reactome.load_pathways(tmp_path, range(1, 16))
    assert sets == {stable: ("Leaf", tuple(range(1, 16)))}
    assert qc["human_mapping_rows"] == 17
    assert qc["nonnumeric_mapping_rows"] == 1
    assert qc["human_gene_pathway_pairs"] == 15
    (tmp_path / "NCBI2Reactome_All_Levels.v98.txt").write_bytes(b"corrupted")
    with pytest.raises(ValueError, match="checksum"):
        reactome.load_pathways(tmp_path, range(1, 16))
