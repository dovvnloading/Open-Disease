"""Regression checks for allele orientation, ambiguous strands and sparse BGZF ranges."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from india_gwas_panukb import read_site
from india_gwas_portability import orient


def test_opposite_allele_and_coordinates():
    panel = dict(chromosome="2", pos="100", risk="C", alleles={"C", "T"})
    external = dict(chromosome="2", base_pair_location="100", effect_allele="T", other_allele="C")
    assert orient(panel, external) == ("aligned", -1)
    assert orient(panel, dict(external, base_pair_location="101")) == ("coordinate_mismatch", None)


def test_strand_flip_is_ambiguous_near_half():
    panel = dict(chromosome="1", pos="90", risk="A", alleles={"A", "T"},
                 india_freq=.408, europe_risk_freq=.498)
    external = dict(chromosome="1", base_pair_location="90", effect_allele="A", other_allele="T",
                    effect_allele_frequency=".495")
    assert orient(panel, external) == ("strand_ambiguous", None)
    assert orient(dict(panel, india_freq=.12, europe_risk_freq=.1),
                  dict(external, effect_allele_frequency=".11")) == ("aligned", 1)


def test_tabix_window_requires_one_exact_allele_pair():
    header = ["chr", "pos", "ref", "alt", "beta_CSA"]
    text = "1\t4\tT\tG\t0.1\n1\t5\tT\tC\t-0.1\n1\t5\tT\tG\t0.9\n1\t6\tA\tG\t0.0"
    assert read_site(header, text, 1, 5, {"T","C"})["beta_CSA"] == "-0.1"
    with pytest.raises(ValueError, match="cover"):
        read_site(header, text, 1, 10, {"T","C"})
    with pytest.raises(ValueError, match="one"):
        read_site(header, text, 1, 5, {"A","C"})
