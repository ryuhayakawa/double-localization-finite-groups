#!/usr/bin/env python3
"""Reproduce the exact fusion-aware S3 star-block certificate.

The computation uses only integer polynomial arithmetic and fractions from
the Python standard library.  It exhausts all boundary comparisons needed in
the paper on the interval 0 <= t <= 2.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Callable

HERE = Path(__file__).resolve().parent
OUT = HERE / "s3_star_certificate.json"
TARGETS = {"S3": Fraction(2)}
# Polynomial tuples list coefficients in ascending powers: (a0, a1, ...).
Poly = tuple[int, ...]


@dataclass(frozen=True)
class GroupData:
    """Representation data used to construct the vertex weights."""

    name: str
    labels: tuple[str, ...]
    dimensions: dict[str, int]
    dual: dict[str, str]
    multiplicity: Callable[[str, str, str], int]


LABELS = ("trivial", "sign", "standard")
# Columns are the identity, transposition, and three-cycle conjugacy classes.
CHARACTERS = {
    "trivial": (1, 1, 1),
    "sign": (1, -1, 1),
    "standard": (2, 0, -1),
}
CLASS_SIZES = (1, 3, 2)


def invariant_multiplicity(left: str, middle: str, right: str) -> int:
    """Return dim Inv(V_left tensor V_middle tensor V_right) for S3."""
    numerator = sum(
        size
        * CHARACTERS[left][class_index]
        * CHARACTERS[middle][class_index]
        * CHARACTERS[right][class_index]
        for class_index, size in enumerate(CLASS_SIZES)
    )
    if numerator % 6:
        raise ArithmeticError("the S3 character average is not integral")
    return numerator // 6


GROUPS = (
    GroupData(
        name="S3",
        labels=LABELS,
        dimensions={"trivial": 1, "sign": 1, "standard": 2},
        dual={label: label for label in LABELS},
        multiplicity=invariant_multiplicity,
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-only", action="store_true")
    return parser.parse_args()


def trim(poly: tuple[int | Fraction, ...]) -> tuple[int | Fraction, ...]:
    """Remove trailing zero coefficients, retaining (0,) for zero."""
    values = list(poly)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def add(first, second):
    size = max(len(first), len(second))
    return trim(
        tuple(
            (first[index] if index < len(first) else 0)
            + (second[index] if index < len(second) else 0)
            for index in range(size)
        )
    )


def scale(poly, scalar):
    return trim(tuple(scalar * value for value in poly))


def subtract(first, second):
    return add(first, scale(second, -1))


def multiply(first, second):
    result = [0] * (len(first) + len(second) - 1)
    for left_index, left in enumerate(first):
        for right_index, right in enumerate(second):
            result[left_index + right_index] += left * right
    return trim(tuple(result))


def derivative(poly):
    if len(poly) == 1:
        return (0,)
    return trim(tuple(index * poly[index] for index in range(1, len(poly))))


def evaluate(poly, value: Fraction) -> Fraction:
    """Evaluate by Horner's rule with exact rational arithmetic."""
    result = Fraction(0)
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def power_coefficients_on_interval(poly, left: Fraction, right: Fraction):
    """Return power-basis coefficients of p(left + (right-left)*x)."""
    width = right - left
    degree = len(poly) - 1
    result = [Fraction(0)] * (degree + 1)
    for source_degree, coefficient in enumerate(poly):
        for target_degree in range(source_degree + 1):
            result[target_degree] += (
                coefficient
                * math.comb(source_degree, target_degree)
                * left ** (source_degree - target_degree)
                * width**target_degree
            )
    return trim(tuple(result))


def bernstein_coefficients(poly, left: Fraction, right: Fraction):
    """Expand the rescaled polynomial in binom(n,k)*x^k*(1-x)^(n-k)."""
    power = power_coefficients_on_interval(poly, left, right)
    degree = len(power) - 1
    return tuple(
        sum(
            power[index]
            * Fraction(math.comb(order, index), math.comb(degree, index))
            for index in range(order + 1)
        )
        for order in range(degree + 1)
    )


def fixed_sign(poly, right: Fraction) -> int | None:
    """Certify a weak sign on [0, right]; None means inconclusive."""
    # Bernstein basis functions are nonnegative on [0, 1].
    coefficients = bernstein_coefficients(poly, Fraction(0), right)
    if all(value >= 0 for value in coefficients):
        return 1
    if all(value <= 0 for value in coefficients):
        return -1
    return None


def nonnegative_by_subdivision(
    poly,
    left: Fraction,
    right: Fraction,
    depth: int = 0,
    maximum_depth: int = 12,
) -> tuple[bool, int, Fraction]:
    """Certify nonnegativity by rational bisection; False is not a disproof."""
    coefficients = bernstein_coefficients(poly, left, right)
    minimum = min(coefficients)
    if minimum >= 0:
        return True, depth, minimum
    if max(coefficients) < 0 or depth >= maximum_depth:
        return False, depth, minimum
    midpoint = (left + right) / 2
    left_result = nonnegative_by_subdivision(
        poly, left, midpoint, depth + 1, maximum_depth
    )
    if not left_result[0]:
        return left_result
    right_result = nonnegative_by_subdivision(
        poly, midpoint, right, depth + 1, maximum_depth
    )
    return (
        right_result[0],
        max(left_result[1], right_result[1]),
        min(left_result[2], right_result[2]),
    )


def local_weight_poly(data: GroupData, triple: tuple[str, str, str]) -> Poly:
    """Return z_t = product of dimensions + t * invariant multiplicity."""
    dimension = math.prod(data.dimensions[label] for label in triple)
    return (dimension, data.multiplicity(*triple))


def star_packet(data: GroupData, boundary: tuple[str, ...]):
    """Return all 27 unnormalized star weights and their partition sum."""
    # Boundary order is (b_1,1, b_1,2, b_2,1, b_2,2, b_3,1, b_3,2).
    weights = []
    for state in itertools.product(data.labels, repeat=3):
        weight = local_weight_poly(data, state)
        for coordinate, label in enumerate(state):
            weight = multiply(
                weight,
                local_weight_poly(
                    data,
                    (
                        data.dual[label],
                        boundary[2 * coordinate],
                        boundary[2 * coordinate + 1],
                    ),
                ),
            )
        weights.append(weight)
    partition = (0,)
    for weight in weights:
        partition = add(partition, weight)
    return tuple(weights), partition


def adjacent_boundary_pairs(data: GroupData):
    """Enumerate all 6 * 3^5 * binom(3,2) unordered one-label changes."""
    labels = data.labels
    for coordinate in range(6):
        other = [index for index in range(6) if index != coordinate]
        for fixed in itertools.product(labels, repeat=5):
            base: list[str | None] = [None] * 6
            for index, value in zip(other, fixed):
                base[index] = value
            for left_index, left in enumerate(labels):
                for right in labels[left_index + 1 :]:
                    first = list(base)
                    second = list(base)
                    first[coordinate] = left
                    second[coordinate] = right
                    yield coordinate, tuple(first), tuple(second)


def fraction_record(value: Fraction) -> dict[str, object]:
    # Decimal values are for display only; verification uses the exact field.
    return {
        "exact": f"{value.numerator}/{value.denominator}",
        "decimal": float(value),
    }


def comparison_polynomial(
    first_packet,
    second_packet,
    changed_star_coordinate: int,
    data: GroupData,
    t_star: Fraction,
):
    """Build N, D, and N'D-ND' for TV(full) + 2*TV(other two labels)."""
    first_weights, first_partition = first_packet
    second_weights, second_partition = second_packet
    state_space = tuple(itertools.product(data.labels, repeat=3))
    # (F1): clear positive partition denominators and certify numerator signs.
    full_differences = tuple(
        subtract(
            multiply(first_weight, second_partition),
            multiply(second_weight, first_partition),
        )
        for first_weight, second_weight in zip(first_weights, second_weights)
    )
    full_signs = tuple(fixed_sign(poly, t_star) for poly in full_differences)
    if None in full_signs:
        return None, "full-sign"
    full_numerator = (0,)
    for sign, poly in zip(full_signs, full_differences):
        full_numerator = add(full_numerator, scale(poly, sign))

    keep = [index for index in range(3) if index != changed_star_coordinate]
    # Sum out the star label incident to the changed boundary coordinate.
    first_marginals = {
        key: (0,) for key in itertools.product(data.labels, repeat=2)
    }
    second_marginals = {
        key: (0,) for key in itertools.product(data.labels, repeat=2)
    }
    for state, first_weight, second_weight in zip(
        state_space, first_weights, second_weights
    ):
        key = (state[keep[0]], state[keep[1]])
        first_marginals[key] = add(first_marginals[key], first_weight)
        second_marginals[key] = add(second_marginals[key], second_weight)
    marginal_differences = tuple(
        subtract(
            multiply(first_marginals[key], second_partition),
            multiply(second_marginals[key], first_partition),
        )
        for key in first_marginals
    )
    marginal_signs = tuple(
        fixed_sign(poly, t_star) for poly in marginal_differences
    )
    if None in marginal_signs:
        return None, "marginal-sign"
    marginal_numerator = (0,)
    for sign, poly in zip(marginal_signs, marginal_differences):
        marginal_numerator = add(marginal_numerator, scale(poly, sign))

    # Fixed signs remove absolute values; the common denominator is 2*Z_b*Z_b'.
    numerator = add(full_numerator, scale(marginal_numerator, 2))
    denominator = scale(multiply(first_partition, second_partition), 2)
    derivative_numerator = subtract(
        multiply(derivative(numerator), denominator),
        multiply(numerator, derivative(denominator)),
    )
    return (numerator, denominator, derivative_numerator), None


def internal_influence_interval(data: GroupData, t_star: Fraction):
    """(F4): certify monotonicity and endpoint maxima for 81 internal TVs."""
    maximum = Fraction(0)
    witnesses = []
    comparisons = 0
    sign_failures = 0
    derivative_failures = 0
    maximum_subdivision_depth = 0
    # outer=(b_1,1,b_1,2), fixed=lambda_3; left/right are two lambda_2 values.
    for outer in itertools.product(data.labels, repeat=2):
        for fixed in data.labels:
            for left_index, left in enumerate(data.labels):
                for right in data.labels[left_index + 1 :]:
                    first_weights = []
                    second_weights = []
                    for output in data.labels:
                        # Only the two vertices incident to the updated edge remain.
                        outer_weight = local_weight_poly(
                            data, (data.dual[output], outer[0], outer[1])
                        )
                        first_weights.append(
                            multiply(
                                local_weight_poly(data, (output, left, fixed)),
                                outer_weight,
                            )
                        )
                        second_weights.append(
                            multiply(
                                local_weight_poly(data, (output, right, fixed)),
                                outer_weight,
                            )
                        )
                    first_total = (0,)
                    second_total = (0,)
                    for weight in first_weights:
                        first_total = add(first_total, weight)
                    for weight in second_weights:
                        second_total = add(second_total, weight)
                    differences = tuple(
                        subtract(
                            multiply(first_weight, second_total),
                            multiply(second_weight, first_total),
                        )
                        for first_weight, second_weight in zip(
                            first_weights, second_weights
                        )
                    )
                    signs = tuple(fixed_sign(poly, t_star) for poly in differences)
                    comparisons += 1
                    if None in signs:
                        sign_failures += 1
                        continue
                    numerator = (0,)
                    for sign, poly in zip(signs, differences):
                        numerator = add(numerator, scale(poly, sign))
                    denominator = scale(multiply(first_total, second_total), 2)
                    derivative_numerator = subtract(
                        multiply(derivative(numerator), denominator),
                        multiply(numerator, derivative(denominator)),
                    )
                    nonnegative, depth, _ = nonnegative_by_subdivision(
                        derivative_numerator, Fraction(0), t_star
                    )
                    maximum_subdivision_depth = max(
                        maximum_subdivision_depth, depth
                    )
                    if not nonnegative:
                        derivative_failures += 1
                        continue
                    # Certified monotonicity makes this endpoint value the maximum.
                    tv = evaluate(numerator, t_star) / evaluate(
                        denominator, t_star
                    )
                    witness = {
                        "outer_pair": list(outer),
                        "fixed_internal_label": fixed,
                        "first_changed_label": left,
                        "second_changed_label": right,
                    }
                    if tv > maximum:
                        maximum = tv
                        witnesses = [witness]
                    elif tv == maximum:
                        witnesses.append(witness)
    return {
        "maximum": maximum,
        "comparison_count": comparisons,
        "fixed_sign_failures": sign_failures,
        "derivative_nonnegativity_failures": derivative_failures,
        "maximum_subdivision_depth": maximum_subdivision_depth,
        "witnesses": witnesses,
    }


def audit_group(data: GroupData, t_star: Fraction) -> dict[str, object]:
    """Exhaust the boundary inventory and assemble the finite certificate."""
    packets = {
        boundary: star_packet(data, boundary)
        for boundary in itertools.product(data.labels, repeat=6)
    }
    maximum = Fraction(0)
    maximum_witnesses = []
    comparison_count = 0
    sign_failures = {"full-sign": 0, "marginal-sign": 0}
    derivative_failures = 0
    maximum_subdivision_depth = 0
    minimum_accepted_bernstein = None
    comparison_cache = {}

    for coordinate, first_boundary, second_boundary in adjacent_boundary_pairs(data):
        comparison_count += 1
        first_packet = packets[first_boundary]
        second_packet = packets[second_boundary]
        # Reuse identical polynomial data, but still count every boundary pair.
        # Integer division maps each pair of boundary labels to its star edge.
        cache_key = (coordinate // 2, first_packet, second_packet)
        if cache_key not in comparison_cache:
            result, failure = comparison_polynomial(
                first_packet,
                second_packet,
                coordinate // 2,
                data,
                t_star,
            )
            if failure is None:
                numerator, denominator, derivative_numerator = result
                # (F2): D>0, so N'D-ND' >= 0 certifies monotonicity of N/D.
                nonnegative, depth, local_minimum = nonnegative_by_subdivision(
                    derivative_numerator, Fraction(0), t_star
                )
                # (F3): exact endpoint value; accepted below only if F2 passes.
                value = evaluate(numerator, t_star) / evaluate(
                    denominator, t_star
                )
                comparison_cache[cache_key] = (
                    failure,
                    nonnegative,
                    depth,
                    local_minimum,
                    value,
                )
            else:
                comparison_cache[cache_key] = (failure, False, 0, Fraction(0), None)
        failure, nonnegative, depth, local_minimum, value = comparison_cache[
            cache_key
        ]
        if failure is not None:
            sign_failures[failure] += 1
            continue
        maximum_subdivision_depth = max(maximum_subdivision_depth, depth)
        if not nonnegative:
            derivative_failures += 1
            continue
        if minimum_accepted_bernstein is None:
            minimum_accepted_bernstein = local_minimum
        else:
            minimum_accepted_bernstein = min(
                minimum_accepted_bernstein, local_minimum
            )
        witness = {
            "changed_boundary_coordinate": coordinate,
            "first_boundary": list(first_boundary),
            "second_boundary": list(second_boundary),
        }
        if value > maximum:
            maximum = value
            maximum_witnesses = [witness]
        elif value == maximum:
            maximum_witnesses.append(witness)

    if minimum_accepted_bernstein is None:
        raise ValueError("no boundary comparison passed verification")
    internal_audit = internal_influence_interval(data, t_star)
    internal_influence = internal_audit["maximum"]
    internal_witnesses = internal_audit["witnesses"]
    # Apply the paper's analytic bounds, not numerical eigenvalue estimates.
    internal_floor = 1 - 2 * internal_influence
    block_curvature = 2 - 4 * maximum
    comparison_factor = internal_floor / 2
    edge_gap_floor = block_curvature * comparison_factor
    return {
        "group": data.name,
        "t_interval": f"0<=t<={t_star}",
        "state_count": len(data.labels) ** 3,
        "boundary_count": len(packets),
        "boundary_comparison_count": comparison_count,
        "distinct_polynomial_comparison_count": len(comparison_cache),
        "fixed_sign_failures": sign_failures,
        "derivative_nonnegativity_failures": derivative_failures,
        "maximum_subdivision_depth": maximum_subdivision_depth,
        "minimum_accepted_bernstein_coefficient": fraction_record(
            minimum_accepted_bernstein
        ),
        "maximum_residual_coupling_bound": fraction_record(maximum),
        "maximum_witness_count": len(maximum_witnesses),
        "maximum_witnesses": maximum_witnesses[:4],
        "star_block_curvature_floor": fraction_record(block_curvature),
        "internal_influence_comparison_count": internal_audit[
            "comparison_count"
        ],
        "internal_influence_fixed_sign_failures": internal_audit[
            "fixed_sign_failures"
        ],
        "internal_influence_derivative_nonnegativity_failures": internal_audit[
            "derivative_nonnegativity_failures"
        ],
        "internal_influence_maximum_subdivision_depth": internal_audit[
            "maximum_subdivision_depth"
        ],
        "maximum_internal_one_neighbor_influence": fraction_record(
            internal_influence
        ),
        "one_site_endpoint_negative_control": {
            "maximum_single_neighbor_influence": fraction_record(
                internal_influence
            ),
            "degree_four_row_sum": fraction_record(4 * internal_influence),
        },
        "internal_witness_count": len(internal_witnesses),
        "internal_witnesses": internal_witnesses[:4],
        "conditional_star_heat_bath_floor": fraction_record(internal_floor),
        "edge_to_star_comparison_factor": fraction_record(comparison_factor),
        "edge_label_gap_floor": fraction_record(edge_gap_floor),
    }


def build_receipt() -> dict[str, object]:
    """Recompute the certificate and package its results for JSON output."""
    groups = []
    for data in GROUPS:
        if data.name in TARGETS:
            groups.append(audit_group(data, TARGETS[data.name]))
    return {
        "schema_version": 1,
        "certificate": "fusion-aware S3 star-block interval certificate",
        "generated_by": "verify_s3_star_certificate.py",
        "software_requirements": {
            "python": ">=3.9",
            "third_party_packages": [],
        },
        "star_weight": "z_t(l1,l2,l3) prod_i z_t(li^*,b_i,c_i)",
        "residual_transport_bound": "W1 <= TV(full)+2 TV(opposite marginal)",
        "path_coupling_geometry": "kappa_star >= 2-4 W1",
        "groups": groups,
    }


def verify(receipt: dict[str, object]) -> None:
    """Reject failed checks and changes to the published values or inventory."""
    expected_targets = {"S3": "0<=t<=2"}
    groups = receipt["groups"]
    if len(groups) != 1 or groups[0]["group"] != "S3":
        raise ValueError("the public certificate must contain exactly S3")
    for group in groups:
        if group["t_interval"] != expected_targets[group["group"]]:
            raise ValueError("target interval changed")
        if any(group["fixed_sign_failures"].values()):
            raise ValueError("a probability difference changes unresolved sign")
        if group["derivative_nonnegativity_failures"]:
            raise ValueError("a residual coupling bound is not certified monotone")
        if group["internal_influence_fixed_sign_failures"]:
            raise ValueError("an internal influence difference changes unresolved sign")
        if group["internal_influence_derivative_nonnegativity_failures"]:
            raise ValueError("an internal influence is not certified monotone")
        if Fraction(group["maximum_residual_coupling_bound"]["exact"]) >= Fraction(
            1, 2
        ):
            raise ValueError("star path-coupling curvature is not positive")
        if Fraction(group["conditional_star_heat_bath_floor"]["exact"]) <= 0:
            raise ValueError("conditional star comparison floor is not positive")
        if Fraction(group["edge_label_gap_floor"]["exact"]) <= 0:
            raise ValueError("edge label floor is not positive")
        if Fraction(
            group["one_site_endpoint_negative_control"]["degree_four_row_sum"][
                "exact"
            ]
        ) <= 1:
            raise ValueError("target endpoint does not exceed the one-site region")

        expected_exact = {
            "maximum_residual_coupling_bound": Fraction(2818, 7395),
            "star_block_curvature_floor": Fraction(3518, 7395),
            "maximum_internal_one_neighbor_influence": Fraction(12, 35),
            "conditional_star_heat_bath_floor": Fraction(11, 35),
            "edge_to_star_comparison_factor": Fraction(11, 70),
            "edge_label_gap_floor": Fraction(19349, 258825),
        }
        for key, expected in expected_exact.items():
            if Fraction(group[key]["exact"]) != expected:
                raise ValueError(f"unexpected exact value for {key}")

        expected_counts = {
            "state_count": 27,
            "boundary_count": 729,
            "boundary_comparison_count": 4374,
            "distinct_polynomial_comparison_count": 288,
            "internal_influence_comparison_count": 81,
        }
        for key, expected in expected_counts.items():
            if group[key] != expected:
                raise ValueError(f"unexpected inventory count for {key}")


def main() -> None:
    args = parse_args()
    receipt = build_receipt()
    verify(receipt)
    rendered = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.check_only:
        # Recompute and compare without modifying the stored receipt.
        if not OUT.exists():
            raise SystemExit(f"missing receipt: {OUT}")
        if OUT.read_text() != rendered:
            raise SystemExit("stored receipt differs from recomputation")
    else:
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
