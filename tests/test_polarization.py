"""Non-HFS polarization regressions; references never call reconstruction helpers."""
from pathlib import Path
import os
import sys

os.environ.setdefault("MPMATH_NOGMPY", "1")
os.environ.setdefault("SYMPY_GROUND_TYPES", "python")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pytest
import AtomSpect

# Zero-field He I fine-structure energies, in cm^-1, copied from the
# E_ground / E_excited inputs in Examples/He_ex1/He_Example.py (lines 88-89).
# EG: lower 1s2p triplet P term (L=1, S=1), ordered by J=0,1,2.
# EE: upper 1s4d triplet D term (L=2, S=1), ordered by J=1,2,3.
# These are the repository's example data, not independently validated here.
EG = np.array([169087.8308131, 169086.8428979, 169086.7664725])
EE = np.array([191444.49804915, 191444.47952984, 191444.47832914])

# Numerical crossing field (tesla) for the lower term with the EG inputs above
# and AtomSpect's existing g_L, g_S and muBcm constants. In the coupled basis,
# build H(B) = E0 + B*muBcm*Psi @ diag(g_L*mL + g_S*mS) @ Psi.T,
# where Psi = stateprocess(1, 1, EG, 0)[0]. Diagonalize each fixed-M block.
# The second-lowest eigenvalue of M=-1 crosses the lowest of M=+1 here.
# Reproduce by subtracting EG[-1]*identity from E0 for numerical stability,
# then use scipy.optimize.brentq on their energy difference over [0.05, 0.065]
# tesla (xtol=1e-15). This gives 0.05784143464329161 with the current inputs.
# The decimal places specify a regression location, not experimental precision;
# tests probe this location and its neighborhood to exercise crossing behavior.
CROSSING = .05784143464329161
GEOMETRIES = [dict(Bangle=90), dict(Bangle=45), dict(Bangle=0),
              dict(Bangle=37, gamma=23, Filter=True, T_Alpha=.65, T_Beta=.35)]


def he_levels(field, sort=False):
    g = AtomSpect.stateprocess(1, 1, EG, field, sortE=sort)
    e = AtomSpect.stateprocess(2, 1, EE, field, sortE=sort)
    return g, AtomSpect.Zeeman_func(g, field), e, AtomSpect.Zeeman_func(e, field)


def assert_e1_labels(g, e, result, geometry=None):
    """Check returned transformation against M operators and separate q channels."""
    geometry = geometry or {}
    transform = np.asarray(result[4])
    lower_m = np.tile([s[1] for s in g[-1]], len(e[-1]))
    upper_m = np.repeat([s[1] for s in e[-1]], len(g[-1]))
    probabilities = np.abs(transform)**2
    np.testing.assert_allclose(probabilities.sum(axis=0), 1, atol=1e-13)
    means = []
    for magnetic_m in [lower_m, upper_m]:
        mean = magnetic_m @ probabilities
        variance = (magnetic_m[:, None] - mean)**2
        assert np.max(np.sum(probabilities * variance, axis=0)) < 1e-20
        means.append(mean)
    expected = np.rint(means[1] - means[0]).astype(int)
    expected[np.abs(expected) > 1] = 3
    np.testing.assert_array_equal(np.asarray(result[2])[:, 1], expected)
    assert not np.any(expected[result[0] > 1e-8] == 3)

    # Build three dipole operators independently of the signal's basis factors.
    raw = np.zeros(len(expected))
    for q in [-1, 0, 1]:
        operator = np.array([[AtomSpect.dipolestr(u, l) if u[1] - l[1] == q else 0
                              for l in g[-1]] for u in e[-1]])
        amplitude = transform.T @ operator.ravel()
        assert np.max(np.abs(amplitude[expected != q]), initial=0) < 1e-12
        factor = AtomSpect.polarization(q, 0, geometry.get("Bangle", 90),
                                geometry.get("gamma", 0),
                                Pol_Filter=geometry.get("Filter", False),
                                TAlpha=geometry.get("T_Alpha", 1),
                                TBeta=geometry.get("T_Beta", 1))[0]
        raw += np.abs(amplitude * factor)**2
    np.testing.assert_allclose(result[8], raw, rtol=1e-12, atol=1e-15)


@pytest.mark.parametrize("field", [.5, 1.5, 3., 5.])
@pytest.mark.parametrize("sort", [False, True])
@pytest.mark.parametrize("geometry", GEOMETRIES)
def test_he_labels_follow_input_eigenpairs(field, sort, geometry):
    g, zg, e, ze = he_levels(field, sort)
    result = AtomSpect.Zeeman_signal(g, zg, e, ze, **geometry)
    mg = np.array([s[1] for s in g[-1]]) @ np.abs(np.asarray(zg[1]))**2
    me = np.array([s[1] for s in e[-1]]) @ np.abs(np.asarray(ze[1]))**2
    expected = np.rint((me[:, None] - mg).ravel()).astype(int)
    expected[np.abs(expected) > 1] = 3
    np.testing.assert_array_equal(np.asarray(result[2])[:, 1], expected)
    assert_e1_labels(g, e, result, geometry)


@pytest.mark.parametrize("field", [0., 1e-9, 1e-6, CROSSING-1e-13,
                                   CROSSING, CROSSING+1e-13])
@pytest.mark.parametrize("sort", [False, True])
def test_crossing_and_small_field_spherical_channels(field, sort):
    g, zg, e, ze = he_levels(field, sort)
    assert_e1_labels(g, e, AtomSpect.Zeeman_signal(g, zg, e, ze))


def test_nondegenerate_eigenpair_and_basis_permutations():
    g, zg, e, ze = he_levels(.5)
    original = AtomSpect.Zeeman_signal(g, zg, e, ze)
    rng = np.random.default_rng(42)
    pg, pe = rng.permutation(len(zg[0])), rng.permutation(len(ze[0]))
    reordered = AtomSpect.Zeeman_signal(g, [zg[0][pg], np.asarray(zg[1])[:, pg]],
                               e, [ze[0][pe], np.asarray(ze[1])[:, pe]])
    indices = (pe[:, None] * len(pg) + pg).ravel()
    np.testing.assert_array_equal(reordered[2], np.asarray(original[2])[indices])
    np.testing.assert_allclose(reordered[0], original[0][indices], rtol=0, atol=1e-11)
    gp, ep = list(g), list(e)
    gp[-1], ep[-1] = [g[-1][i] for i in pg], [e[-1][i] for i in pe]
    rows = AtomSpect.Zeeman_signal(gp, [zg[0], np.asarray(zg[1])[pg]],
                          ep, [ze[0], np.asarray(ze[1])[pe]])
    np.testing.assert_array_equal(rows[2], original[2])
    np.testing.assert_allclose(rows[0], original[0], rtol=0, atol=1e-11)


@pytest.mark.parametrize("field", [0., 1e-9, .5])
def test_j_zero_to_one_ratio(field):
    g = AtomSpect.stateprocess(0, 0, np.array([0.]), field)
    e = AtomSpect.stateprocess(1, 0, np.array([20000.]), field)
    result = AtomSpect.Zeeman_signal(g, AtomSpect.Zeeman_func(g, field), e, AtomSpect.Zeeman_func(e, field))
    by_q = {int(r[1]): s for r, s in zip(result[2], result[0])}
    np.testing.assert_allclose([by_q[-1], by_q[0], by_q[1]], [.25, .5, .25], atol=1e-14)


def test_degenerate_subspace_is_resolved_without_mutating_input():
    g = AtomSpect.stateprocess(0, 0, np.array([0.]), 0.)
    e = AtomSpect.stateprocess(1, 0, np.array([20000.]), 0.)
    zg, ze = AtomSpect.Zeeman_func(g, 0.), AtomSpect.Zeeman_func(e, 0.)
    rotation = np.linalg.qr(np.random.default_rng(42).normal(size=(3, 3)))[0]
    for vectors in [np.asarray(ze[1])[:, [2, 0, 1]], np.asarray(ze[1]) @ rotation]:
        supplied = [ze[0].copy(), vectors.copy()]
        result = AtomSpect.Zeeman_signal(g, zg, e, supplied)
        assert_e1_labels(g, e, result)
        np.testing.assert_array_equal(supplied[0], ze[0])
        np.testing.assert_array_equal(supplied[1], vectors)
        # Equal-energy columns have no preserved individual identity. Test the
        # physically observable partition and total at their common wavelength.
        by_q = {int(r[1]): s for r, s in zip(result[2], result[0])}
        np.testing.assert_allclose([by_q[-1], by_q[0], by_q[1]], [.25, .5, .25], atol=1e-14)
        np.testing.assert_allclose(result[7], 2/3, atol=1e-14)
        np.testing.assert_allclose(result[1], result[1][0], atol=1e-14)


def test_film_partition_and_exported_sigma_signs():
    deck = dict(s_ground=1, s_excited=1, l_ground=1, l_excited=2,
                E_ground=EG, E_excited=EE, Bmag=3.2, b_angle=90,
                TAlpha=.65, TBeta=.35, amu=4, Temp=300,
                Convfxn="GaussianInstrum", specstep=.005,
                spec_window=[446., 448.], fxnwindow=.1)
    result = AtomSpect.AtomSpect_Main(deck)
    labels = np.asarray(result["rad"])[:, 1]
    for key, q, total in [("pi_sticks", 0, .65), ("sigmaplus_sticks", -1, .175),
                          ("sigmaminus_sticks", 1, .175), ("other_sticks", 3, 0.)]:
        np.testing.assert_array_equal(result[key][1], result["signal"][labels == q])
        np.testing.assert_allclose(np.sum(result[key][1]), total, atol=1e-12)


@pytest.mark.parametrize("geometry", [dict(Bangle=90),
    dict(Bangle=90, Filter=True, gamma=0),
    dict(Bangle=90, Filter=True, gamma=90),
    dict(Bangle=37, Filter=True, gamma=23, T_Alpha=.65, T_Beta=.35)])
def test_zero_field_he_preserves_strength_at_each_degenerate_line(geometry):
    g, zg, e, ze = he_levels(0.)
    result = AtomSpect.Zeeman_signal(g, zg, e, ze, **geometry)
    reference = np.zeros(len(result[0]))
    for q in [-1, 0, 1]:
        operator = np.array([[AtomSpect.dipolestr(u, l) if u[1]-l[1] == q else 0
                              for l in g[-1]] for u in e[-1]])
        amplitude = np.asarray(ze[1]).T @ operator @ np.asarray(zg[1])
        factor = AtomSpect.polarization(q, 0, geometry["Bangle"], geometry.get("gamma", 0),
                                Pol_Filter=geometry.get("Filter", False),
                                TAlpha=geometry.get("T_Alpha", 1),
                                TBeta=geometry.get("T_Beta", 1))[0]
        reference += np.abs(factor * amplitude.ravel())**2
    # Eight decimal places in nm group only coincident zero-field components;
    # distinct fine-structure lines are much farther apart in this fixture.
    line = np.round(result[6], 8)
    assert len(np.unique(line)) == 9
    for wavelength in np.unique(line):
        mask = line == wavelength
        np.testing.assert_allclose(result[8][mask].sum(), reference[mask].sum(),
                                   rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(result[0][mask].sum(), reference[mask].sum()/reference.sum(),
                                   rtol=1e-12, atol=1e-14)
