"""Sharp array-omission bounds for fixed scores with no case–comparator ties."""

from fractions import Fraction


def _stratum_options(cases, comparators, higher):
    """Generate possible omission counts, remaining pairs, and extreme win counts.

    If R pairs remain out of P, at most P-R wins can be lost and at most R
    can remain. Both bounds are attainable by rankings consistent with U wins.
    To retain as many wins as possible, start with all omitted cases below
    both retained groups and all omitted comparators above them. Any U up to
    R is realizable within the retained group alone. For U > R, keep that
    group perfectly ordered and swap adjacent omitted/other items toward a
    full ordering, adding one win at a time up to P. Reversing score order
    gives the lower bound. This argument concerns possible rankings, not the
    identity of the arrays that would attain a bound in the original study.
    """
    for removed_cases in range(cases):
        for removed_comparators in range(comparators):
            remaining = (cases - removed_cases) * (comparators - removed_comparators)
            lost = cases * comparators - remaining
            yield (removed_cases + removed_comparators, remaining,
                   max(0, higher - lost), min(higher, remaining))


def _fraction(value):
    return {"numerator": value.numerator, "denominator": value.denominator}


def omission_summary(strata):
    """Bound the pooled fraction of within-stratum case wins after array omissions.

    Supply (cases, comparators, higher) integer triples. The caller must establish
    there are no case–comparator score ties. At least one array in each arm of
    every stratum is retained. Bounds include all score orders consistent with
    the counts, use fixed scores, and are for exactly the stated omission count.
    """
    strata = list(strata)
    if not strata:
        raise ValueError("At least one stratum is required")
    # For given totals of omitted arrays and remaining pairs, only the smallest
    # and largest surviving win count can affect any later stratum's bounds.
    states = {(0, 0): (0, 0)}
    total_pairs = total_higher = 0
    for cases, comparators, higher in strata:
        if any(type(value) is not int for value in (cases, comparators, higher)):
            raise ValueError("Stratum counts must be integers")
        if min(cases, comparators) < 1 or not 0 <= higher <= cases * comparators:
            raise ValueError("Strata need both arms and a valid count of case wins")
        total_pairs += cases * comparators
        total_higher += higher
        combined = {}
        for removed, pairs, low, high in _stratum_options(cases, comparators, higher):
            for (prior_removed, prior_pairs), (prior_low, prior_high) in states.items():
                key = (prior_removed + removed, prior_pairs + pairs)
                low_sum, high_sum = prior_low + low, prior_high + high
                if key in combined:
                    low_sum = min(low_sum, combined[key][0])
                    high_sum = max(high_sum, combined[key][1])
                combined[key] = (low_sum, high_sum)
        states = combined

    extrema = {}
    for (removed, pairs), (low, high) in states.items():
        lower, upper = Fraction(low, pairs), Fraction(high, pairs)
        if removed in extrema:
            lower = min(lower, extrema[removed][0])
            upper = max(upper, extrema[removed][1])
        extrema[removed] = (lower, upper)

    def bounds(removed):
        if removed not in extrema:
            return None
        lower, upper = extrema[removed]
        return {"omissions": removed, "lower": _fraction(lower), "upper": _fraction(upper)}

    baseline = Fraction(total_higher, total_pairs)
    half = Fraction(1, 2)
    if baseline > half:
        direction = "above_half"
        threshold = min((k for k, (lower, _) in extrema.items() if lower <= half), default=None)
    elif baseline < half:
        direction = "below_half"
        threshold = min((k for k, (_, upper) in extrema.items() if upper >= half), default=None)
    else:
        direction, threshold = "at_half", 0
    return {
        "baseline": {"higher": total_higher, "pairs": total_pairs, "fraction": _fraction(baseline)},
        "direction": direction,
        "one_omission": bounds(1),
        "minimum_omissions_to_lose_direction": threshold,
        "just_before_threshold": bounds(threshold - 1) if threshold is not None else None,
        "at_threshold": bounds(threshold),
    }
