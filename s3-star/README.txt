Exact S3 Star-Block Certificate
===============================

This directory is the self-contained ancillary package for the finite
calculation in the appendix "Exact Finite Certificates for the S3 Star
Block."  It may be copied verbatim into a submission archive or a separate
public repository.

Scope
-----

The package verifies the exhaustive finite polynomial statement used in the
paper.  For the three irreducible representations of S3 and 0 <= t <= 2, it
checks all conditioned-star boundary comparisons and the internal one-site
influence inventory.  These are checks (F1)-(F4) in Appendix C.2, as
described below.  It does not reprove the analytic path-coupling,
edge-to-star comparison, or the full Davies-gap theorem; those arguments
remain in the paper.

Files
-----

verify_s3_star_certificate.py
    Standalone exact-arithmetic verifier.  The S3 character data and fusion
    multiplicities are included in the file.

s3_star_certificate.json
    Machine-readable output produced by the verifier.

SHA256SUMS
    SHA-256 hashes of this README, the verifier, the stored output, and LICENSE.

LICENSE
    MIT License for this ancillary package.

License
-------

This ancillary package is distributed under the MIT License; see LICENSE.
Copyright (c) 2026 Ryu Hayakawa, Angus Southwell, Caesnan M. G. Leditto,
Kuo-Chin Chen, and Min-Hsiu Hsieh.
This license applies to the package, not to the accompanying paper.

Requirements
------------

Python 3.9 or newer.  No third-party Python packages, network access, or
external data files are required.

Reproduction
------------

From this directory, run

    python3 verify_s3_star_certificate.py --check-only

The command recomputes the complete inventory using integer polynomials and
fractions, checks the claimed exact values, and compares the result byte for
byte with s3_star_certificate.json.  It exits with status zero exactly when
all checks pass.  The JSON record is also printed to standard output.
The --check-only mode does not overwrite the stored record.

A successful run confirms that the implemented checks pass and reproduce
the stored record.  Checking that the implementation represents the
manuscript's distributions and covers its required cases is a separate
part of the mathematical review.

To regenerate the stored record intentionally, run

    python3 verify_s3_star_certificate.py

To verify the distributed files before running the computation, use

    shasum -a 256 -c SHA256SUMS

This checks file integrity, not mathematical correctness.  After editing
any listed file, finalize the package and update SHA256SUMS before release.

Certified inventory
-------------------

The verifier checks:

* 27 internal states;
* 729 boundary conditions;
* all 4374 one-boundary-label comparisons;
* 81 internal one-neighbor influence comparisons;
* no unresolved sign or monotonicity certificates on 0 <= t <= 2.

All boundary comparisons are enumerated.  Comparisons with identical
weight and partition polynomials on both sides and the same coordinate
summed out for the two-label marginal reuse a single verification.  This
leaves 288 distinct polynomial checks covering all 4374 comparisons; no
boundary comparison is omitted from the coverage.

Checks corresponding to Appendix C.2
------------------------------------

(F1) Fixed signs.
     For each boundary comparison, the numerator of each state-probability
     difference, and of each two-label marginal difference, has a fixed
     sign throughout [0,2].  Exact Bernstein coefficients certify this,
     allowing the absolute values in the total-variation sums to be removed.

(F2) Boundary-response derivatives.
     Write the residual transport bound as N(t)/D(t), with D(t) > 0.
     The Bernstein coefficients of N'(t)D(t) - N(t)D'(t) are nonnegative
     on [0,2] for every distinct boundary comparison.  Thus each bound is
     nondecreasing, and its interval maximum is attained at t = 2.

(F3) Boundary endpoint maximum.
     Exact evaluation at t = 2 gives the maximum residual transport bound
     2818/7395 over all boundary comparisons.  Together with (F1)-(F2),
     this determines B_2 on the full interval.  The coupling argument in
     the paper turns this transport bound into a Wasserstein upper bound.

(F4) Internal-response checks.
     For all 81 internal comparisons, the probability-difference
     numerators have fixed signs and the resulting total-variation
     derivative numerators have nonnegative Bernstein coefficients on
     [0,2].  Exact endpoint evaluation gives the maximum internal
     one-neighbor influence c_int,2 = 12/35 on the full interval.

Inputs to the label-chain proposition
------------------------------------

The paper uses the certificate with T = 2 and just two inputs:

    B_T = B_2                 = 2818/7395 < 1/2
    delta_int,T = delta_int,2 = 11/35 > 0

The first bounds the Hamming Wasserstein distance between conditioned-star
distributions whose boundary labels differ in one coordinate.  The second
is a lower bound on the internal single-edge heat-bath gap for every fixed
boundary condition, uniformly for 0 <= t <= 2.

The internal influence c_int,2 = 12/35 is an auxiliary quantity used to
obtain delta_int,2 = 1 - 2*c_int,2 via the paper's Dobrushin argument.
The record also reports the following arithmetic consequences, not
additional inputs or independently computed spectral gaps:

    star-coupling contraction margin  2 - 4*B_2                 = 3518/7395
    edge-to-star comparison factor   delta_int,2 / 2           = 11/70
    edge-label spectral-gap floor    (delta_int,2/2)*(2-4*B_2) = 19349/258825

The validity of these gap bounds uses the analytic arguments in the paper;
the verifier computes their numerical values from the certified constants.

Method
------

Every conditional weight is represented as an integer polynomial in t.  The
signs of the probability differences and the monotonicity of the resulting
rational bounds are certified with exact Bernstein coefficients, with exact
subdivision available when needed; no subdivision is needed for this S3
inventory.  All endpoint evaluations and reported
constants use fractions rather than floating-point arithmetic.  Decimal
values in the JSON file are included only for readability.

Provenance
----------

This public package is the S3-only, dependency-free extraction of the
calculation used to produce the constants in the accompanying manuscript.
The machine-readable receipt records the required Python version and the
absence of third-party dependencies.
