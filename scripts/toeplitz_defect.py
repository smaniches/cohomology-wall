#!/usr/bin/env python3
"""
toeplitz_defect.py -- canonical evaluator for the v0.4 source-target defect law.

This is the single source of truth for the closed-form quantities of Theorem 5.1
(general Toeplitz determinantal saturation) and Corollary 5.2 (the corrected
equal-charge onset rho(b)). It is pure arithmetic: no numpy, no scipy, no
python-flint. It is what the verifier checks finite-field ranks against, and what
the b=7 runner compares its measured rank to.

For D = (-a_1, -a_2, -a_3, b) with a_i >= 4 and b >= 2, the Koszul connecting map
identifies (Kunneth + Serre duality, factor by factor) with

    mu_{b,a} : H^0(Y, O(a-2))^{(b+1)} -> H^0(Y, O(a))^{(b-1)},   Y = (P^1)^3,

whose source and target dimensions are

    E = (b+1) * prod_i (a_i - 1),     F = (b-1) * prod_i (a_i + 1).

In the acyclic range a_i >= 2b+1 the global-section complex of the Buchsbaum-Rim
resolution is exact, the cokernel length is C_b = 48 * C(b+1,3), and therefore

    def(mu_{b,a}) = C_b                 if E >= F
                  = C_b + (E - F)       if E <  F.

rho(b) is the first equal-charge depth d >= 2b+1 at which E >= F, i.e. the first
depth at which the defect locks onto the ceiling C_b. tau(b) is the OLD
cubic-crossing reference; it agrees with rho(b) through b=6 and fails from b=7
(tau(7)=20, rho(7)=21).
"""
from math import comb, ceil

__all__ = ["E", "F", "C_b", "min_EF", "defect_acyclic", "in_acyclic_range",
           "rho", "tau", "ceiling_holds"]


def E(b, a):
    """Source dimension (b+1) * prod_i (a_i - 1)."""
    p = 1
    for x in a:
        p *= (x - 1)
    return (b + 1) * p


def F(b, a):
    """Target dimension (b-1) * prod_i (a_i + 1)."""
    p = 1
    for x in a:
        p *= (x + 1)
    return (b - 1) * p


def C_b(b):
    """Local Buchsbaum-Rim cokernel length 48 * C(b+1, 3)."""
    return 48 * comb(b + 1, 3)


def min_EF(b, a):
    return min(E(b, a), F(b, a))


def in_acyclic_range(b, a):
    """Theorem 5.1 hypothesis: every charge at least 2b+1."""
    return all(x >= 2 * b + 1 for x in a)


def defect_acyclic(b, a):
    """Theorem 5.1 piecewise defect. Defined only for a_i >= 2b+1."""
    if not in_acyclic_range(b, a):
        raise ValueError(
            f"defect_acyclic is the Theorem 5.1 formula, valid only in the acyclic "
            f"range a_i >= 2b+1; got b={b}, a={tuple(a)} (need a_i >= {2*b+1})."
        )
    e, f, c = E(b, a), F(b, a), C_b(b)
    return c if e >= f else c + e - f


def ceiling_holds(b, a):
    """True iff the defect equals the Thom-Porteous ceiling C_b in the acyclic range:
    a_i >= 2b+1 for all i and (b+1) prod(a_i-1) >= (b-1) prod(a_i+1)."""
    return in_acyclic_range(b, a) and E(b, a) >= F(b, a)


def rho(b):
    """Corrected equal-charge onset: least d >= 2b+1 with (b+1)(d-1)^3 >= (b-1)(d+1)^3."""
    d = 2 * b + 1
    while not ((b + 1) * (d - 1) ** 3 >= (b - 1) * (d + 1) ** 3):
        d += 1
    return d


def tau(b):
    """Old cubic-crossing reference (superseded). tau(b) = (2b-1) + ceil(2*(b(b-1))^(1/3))."""
    return (2 * b - 1) + ceil(2 * (b * (b - 1)) ** (1.0 / 3.0))


if __name__ == "__main__":
    print("b | rho(b) | tau(b) | agree | C_b")
    print("-" * 40)
    for b in range(2, 9):
        print(f"{b} | {rho(b):6d} | {tau(b):6d} | {str(rho(b) == tau(b)):5s} | {C_b(b)}")
    print()
    for b, d in [(7, 20), (7, 21)]:
        a = (d, d, d)
        e, f = E(b, a), F(b, a)
        print(f"b={b} d={d}: E={e} F={f} C_b={C_b(b)} "
              f"branch={'C_b' if e >= f else 'C_b+E-F'} defect={defect_acyclic(b, a)}")
