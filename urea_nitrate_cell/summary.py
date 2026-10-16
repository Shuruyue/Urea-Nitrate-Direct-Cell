"""主要結果：全部由模型重算，回傳一個 dict（scripts/summary.py 會印出並寫成 JSON）。"""
import numpy as np
from . import model as m
from . import calibration as cal
from .yao2025 import REPORTED, FIG14D_VCELL


def compute():
    N = {}
    th = m.thermo_check()
    N["thermo"] = dict(th)
    for k in ("UOR_lit", "UOR_calc", "UOR_CO3"):
        E, G = m.cell_E0(m.E0[k])
        N["thermo"][f"Ecell_{k}"] = E
        N["thermo"][f"dG_{k}"] = G
    N["thermo"]["dH_urea_electrolysis_kJ"] = -393.51 - (-317.65) - (-285.83)
    N["thermo"]["nernst_1000x_mV"] = float(m.nernst_shift(1e-3) * 1e3)
    N["thermo"]["NiOOH_RHE"] = 0.49 + 0.828          # 0.49 V vs SHE（Anuratha 2022）→ vs RHE, pH 14
    N["jlim"] = {f"{C:g}": m.j_lim_no3(C, m.S1["delta_um"]) for C in (1, 60 / 14.007, 10, 100)}
    op = {k: m.solve_cell(10.0, p) for k, p in list(m.SCEN.items()) + [("OER_only", m.OER_ONLY)]}
    N["op_10mA"] = op
    js = np.linspace(0.5, 60, 1200)
    P3 = [(j, j * m.solve_cell(j, m.SCEN["S3"])["V"]) for j in js]
    jmax, pmax = max(P3, key=lambda t: t[1])
    N["S3_peak_power"] = dict(j=float(jmax), P_mW_cm2=float(pmax))
    N["overpotential_S1"] = dict(anode=op["S1"]["Ea"] - m.E0["UOR_lit"],
                                 anode_vs_0072=op["S1"]["Ea"] - m.E0["UOR_calc"],
                                 cathode=m.E0["NO3RR"] - op["S1"]["Ec"])
    N["urea_saving"] = dict(dV=op["S1"]["V"] - op["OER_only"]["V"],
                            pct=(1 - op["S1"]["V"] / op["OER_only"]["V"]) * 100)
    an, ca, cf = cal.fit_anode(), cal.fit_cathode(gamma=1.0), cal.fit_cathode()
    N["calibration"] = dict(anode_Ru=an["used"]["Ru"], anode_rms=an["used"]["rms"],
                            cathode=ca["params"], cathode_std=ca["std"], cathode_rms_by_C=ca["rms_by_C"],
                            gamma_free=cf["params"]["gamma"], gamma_std=cf["std"]["gamma"])
    V0, rows = m.sensitivity()
    N["sensitivity_20mA"] = dict(V0=V0, rows={r[0]: dict(lo=r[1], hi=r[2], dV_lo=r[3], dV_hi=r[4]) for r in rows})
    B = {}
    for n, p in (("OER_only", m.OER_ONLY), ("S1", m.S1), ("S3", m.SCEN["S3"])):
        for s in ("cc", "follow"):
            r = m.simulate_batch(p, s)
            B[f"{n}_{s}"] = dict(t_h=r["t_end"], FE_avg=r["FE_avg"], FE_start=float(r["FE"][0]),
                                 FE_end=float(r["FE"][-1]), kWh_per_kgN_nitrate=r["kWh_per_kgN_no3"],
                                 kWh_per_kgN_total=r["kWh_per_kgN_all"])
    N["batch_60mgN_per_L"] = B
    lam_OH, lam_NO3 = 198.0, 71.4                     # 極限莫耳導電度 S cm2 mol-1（CRC Handbook）
    N["AEM_nitrate_transference"] = {f"{c:g} mM": c * 1e-3 * lam_NO3 / (lam_OH + c * 1e-3 * lam_NO3)
                                     for c in (60 / 14.007, 20, 100)}
    N["yao2025"] = dict(reported=REPORTED, cell_voltage_fig14d=FIG14D_VCELL)
    return N
