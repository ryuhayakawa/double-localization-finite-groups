# Double Localization for Finite Groups

Exact computational certificates for double-localization gap bounds in
finite-group models.

## Current Certificate

The [S3 star certificate](s3-star/README.txt) accompanies:

**Double Localization for Quantum Gibbs Sampler Gaps:
From an Abstract Framework to Finite-Group Models**

Ryu Hayakawa, Angus Southwell, Caesnan M. G. Leditto, Kuo-Chin Chen,
and Min-Hsiu Hsieh.

The self-contained verifier reproduces the finite checks (F1)-(F4) in
Appendix C.2 on the full interval `0 <= t <= 2`. It uses exact integer and
rational arithmetic, checks all 4374 boundary comparisons and 81 internal
comparisons, and compares its output with the stored JSON record.
The analytic path-coupling and spectral-gap arguments remain in the paper.

## Reproduce

Python 3.9 or newer is required. No third-party packages are needed.
From the repository root, run:

```sh
cd s3-star
python3 verify_s3_star_certificate.py --check-only
```

A successful run exits with status zero and prints the verification record.
It does not modify the stored record. To check the package's file integrity:

```sh
shasum -a 256 -c SHA256SUMS
```

Run this command from `s3-star` as well. Hashes check file integrity, not
mathematical correctness. See the [package README](s3-star/README.txt) for
the inventory, exact constants, and scope of the verification.

## Attribution and License

When using this certificate in research, please cite the accompanying paper
and identify the repository version used.

The code and associated documentation are distributed under the
[MIT License](LICENSE). The paper is not covered by this repository license.
