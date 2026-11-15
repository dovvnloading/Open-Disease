from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from melbourne_signature_phase import infer


def test_batch_effect_is_removed_and_known_stage_effect_retained():
    stages = np.array([3,3,3,6,6,6] * 2)
    batches = np.array(['a']*6 + ['b']*6)
    y = (stages == 6)*2. + (batches == 'b')*80. + np.array([-1.,0.,1.] * 4)
    result = infer(y, stages, batches, draws=399)
    assert result['stage6_minus_stage3'] == pytest.approx(2.)
    assert result['phase_partial_r2'] == pytest.approx(.6)
    assert result['by_batch']['a']['stage6_minus_stage3'] == pytest.approx(2.)
    assert result['by_batch']['b']['stage6_minus_stage3'] == pytest.approx(2.)
    assert result['contrast_hc3_95_interval'][0] < 2 < result['contrast_hc3_95_interval'][1]


def test_phase_exactly_aliased_with_batch_raises():
    stages = np.array([3,3,6,6]); batches = np.array(['a','a','b','b'])
    with pytest.raises(ValueError, match='aliased'):
        infer(np.array([1.,2.,3.,4.]), stages, batches, draws=9)
