"""Fast unit tests: hand-checkable relations in the model."""
import numpy as np
import pytest
from urea_nitrate_cell import model as m


def test_standard_potentials():
    th = m.thermo_check()
    assert th["NO3RR"] == pytest.approx(0.708, abs=1e-3)
    assert th["UOR_calc"] == pytest.approx(0.072, abs=1e-3)
    assert th["UOR_CO3"] == pytest.approx(-0.025, abs=2e-3)


def test_cell_potential_and_gibbs_energy():
    E, G = m.cell_E0(0.37)
    assert E == pytest.approx(0.338, abs=1e-3)
    assert G == pytest.approx(-8 * 96.485 * E, rel=1e-6)          # ΔG = −nFE, n = 8


def test_nernst_shift_for_1000x_dilution():
    assert m.nernst_shift(1e-3) * 1e3 == pytest.approx(-22.18, abs=0.01)


def test_tafel_anchor_is_10_mA():
    for p in m.SCEN.values():
        assert m.j_uor(p["E10_UOR"], p) == pytest.approx(10.0, rel=1e-12)


def test_limiting_current_matches_hand_calculation():
    hand = 8 * 96485 * 1.902e-5 * (100e-3 / 1000) / 43e-4 * 1e3
    assert m.j_lim_no3(100, 43) == pytest.approx(hand, rel=1e-12)


def test_s1_cell_voltage_at_10_mA():
    r = m.solve_cell(10.0, m.S1)
    assert r["V"] == pytest.approx(-1.309, abs=2e-3)
    assert m.solve_cell(10.0, m.OER_ONLY)["V"] - r["V"] == pytest.approx(-0.180, abs=2e-3)


def test_faradaic_efficiency_bounds():
    for j in np.linspace(0.5, 40, 20):
        r = m.solve_cell(j, m.S1, 4.3)
        assert 0 <= r["FE_NO3"] <= 1 + 1e-9 and 0 <= r["FE_UOR"] <= 1 + 1e-9


def test_batch_mass_balance():
    r = m.simulate_batch(m.S1, "cc")
    assert r["n_N_no3"] * 8 * m.F / r["Q_no3"] == pytest.approx(1.0, rel=1e-9)
