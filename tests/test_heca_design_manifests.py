"""Check donor exchangeability when technically repeated pairs share a capture design."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_heca_design_manifests import summarize_graph


def test_repeated_pairs_restrict_joint_assignments():
    graph = [{'case-A', 'control-A'}, {'case-A', 'control-A'},
             {'case-B', 'control-B'}, {'case-A'}]
    result = summarize_graph(graph, {'case-A', 'case-B'})
    assert result['eligible_donor_library_links'] == 7
    assert result['mixed_libraries'] == 3
    assert result['unique_mixed_memberships'] == 2
    assert result['libraries_per_distinct_mixed_pair'] == [1, 2]
    assert result['allocations_preserving_mixed_library_case_totals'] == 4
    assert result['allocations_preserving_all_eligible_library_case_totals'] == 2
