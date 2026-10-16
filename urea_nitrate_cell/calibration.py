"""
用 Yao 2025 的 LSV 數據校正 S1 參數（最小平方法）
------------------------------------------------
量測到的電位包含未補償的歐姆壓降 iR_u，所以模型寫成
    陽極：j = j_UOR(E_meas − j R_u) + j_OER(E_meas − j R_u)
    陰極：j = j_bg + j_HER(E_true) + j_NO3(E_true, C)，E_true = E_meas + j R_u
其中 j_bg 為 LSV 的背景電流（電容充電等，只在擬合時使用，電池模型不含）。
每個量測點用 brentq 解出模型電流，再以相對誤差 (j_model − j_data)/(j_data + 3) 做最小平方。

輸出：擬合參數 ± 標準誤差、各曲線的 RMS 相對誤差（scripts/calibrate.py 會寫入 outputs/calibration_result.json）。
"""
import numpy as np
from scipy.optimize import least_squares, brentq
from . import model as m
from .yao2025 import FIG3B, FIG11A_UOR, FIG11A_OER, REPORTED


# ---------------------------------------------------------------------------
# 陽極（Fig. 11a）
# ---------------------------------------------------------------------------
def anode_meas(E, p, Ru, with_uor=True):
    q = dict(p, UOR_on=with_uor)
    f = lambda j: j - (m.j_uor(E - j * 1e-3 * Ru, q) + m.j_oer(E - j * 1e-3 * Ru, q))
    return 0.0 if f(0.0) >= 0 else brentq(f, 0.0, 1e4)


def anode_resid(p, Ru):
    r = []
    for (E, J), uor in ((FIG11A_UOR, True), (FIG11A_OER, False)):
        for e, jj in zip(E, J):
            if jj >= (2.0 if uor else 5.0):      # OER 曲線 5 mA 以下為 Ni(OH)2 氧化前峰，不擬合
                r.append((anode_meas(e, p, Ru, uor) - jj) / (jj + 3))
    return np.array(r)


def fit_anode():
    base = dict(m.S1)
    out = {}
    # (A) 模型實際使用的參數（論文報告值 + OER b = 60 mV/dec），只擬合 R_u
    s = least_squares(lambda x: anode_resid(base, x[0]), [1.0], bounds=([0], [5]))
    out["used"] = dict(Ru=s.x[0], rms=float(np.sqrt(np.mean(s.fun ** 2))))
    # (B) 全部放開：E10_UOR, b_UOR, E10_OER, b_OER, R_u
    def res(x):
        return anode_resid(dict(base, E10_UOR=x[0], b_UOR=x[1], E10_OER=x[2], b_OER=x[3]), x[4])
    s = least_squares(res, [1.33, 0.044, 1.51, 0.06, 1.0],
                      bounds=([1.2, 0.01, 1.4, 0.03, 0], [1.45, 0.2, 1.6, 0.25, 5]))
    cov = np.linalg.inv(s.jac.T @ s.jac) * np.mean(s.fun ** 2)
    out["free"] = dict(zip(["E10_UOR", "b_UOR", "E10_OER", "b_OER", "Ru"], s.x.tolist()))
    out["free_std"] = dict(zip(out["free"], np.sqrt(np.diag(cov)).tolist()))
    out["free"]["rms"] = float(np.sqrt(np.mean(s.fun ** 2)))
    # (C) b_UOR 固定為論文值 44.2 mV/dec
    def res_c(x):
        return anode_resid(dict(base, E10_UOR=x[0], E10_OER=x[1], b_OER=x[2]), x[3])
    s = least_squares(res_c, [1.33, 1.51, 0.06, 1.0], bounds=([1.2, 1.4, 0.03, 0], [1.45, 1.6, 0.25, 5]))
    out["fixed_b"] = dict(zip(["E10_UOR", "E10_OER", "b_OER", "Ru"], s.x.tolist()))
    out["fixed_b"]["rms"] = float(np.sqrt(np.mean(s.fun ** 2)))
    return out


# ---------------------------------------------------------------------------
# 陰極（Fig. 3b，0 / 1 / 10 / 100 mM 用於擬合，50 mM 留作獨立檢驗）
# ---------------------------------------------------------------------------
CATH_KEYS = ["E10_HER", "b_HER", "E10_NO3", "b_NO3", "gamma", "delta_um", "Ru", "j_bg"]


def cathode_meas(E, C, q):
    p = dict(m.S1, **{k: q[k] for k in CATH_KEYS if k not in ("Ru", "j_bg")})
    f = lambda j: j - (q["j_bg"] + m.j_her(E + j * 1e-3 * q["Ru"], p)
                       + (m.j_no3(E + j * 1e-3 * q["Ru"], p, C) if C > 0 else 0.0))
    return brentq(f, 0.0, 1e4)


def cathode_resid(q, concs=(0, 1, 10, 100), jmin=3.0):
    r = []
    for C in concs:
        E, J = FIG3B[C]
        for e, jj in zip(E, J):
            if np.isfinite(jj) and jj >= jmin:
                r.append((cathode_meas(e, C, q) - jj) / (jj + 3))
    return np.array(r)


def fit_cathode(gamma=None):
    keys = [k for k in CATH_KEYS if not (k == "gamma" and gamma is not None)]
    x0 = dict(E10_HER=-0.21, b_HER=0.12, E10_NO3=0.05, b_NO3=0.12, gamma=0.8, delta_um=40,
              Ru=2.0, j_bg=2.0)
    lo = dict(E10_HER=-0.5, b_HER=0.04, E10_NO3=-0.3, b_NO3=0.03, gamma=0.0, delta_um=5, Ru=0, j_bg=0)
    hi = dict(E10_HER=0.0, b_HER=0.3, E10_NO3=0.4, b_NO3=0.3, gamma=1.5, delta_um=300, Ru=10, j_bg=5)

    def unpack(x):
        q = dict(zip(keys, x))
        if gamma is not None:
            q["gamma"] = gamma
        return q
    s = least_squares(lambda x: cathode_resid(unpack(x)), [x0[k] for k in keys],
                      bounds=([lo[k] for k in keys], [hi[k] for k in keys]))
    q = unpack(s.x)
    cov = np.linalg.inv(s.jac.T @ s.jac) * np.mean(s.fun ** 2)
    std = dict(zip(keys, np.sqrt(np.diag(cov)).tolist()))
    rms = {C: float(np.sqrt(np.mean(cathode_resid(q, (C,)) ** 2))) for C in (0, 1, 10, 50, 100)}
    return dict(params=q, std=std, rms_fit=float(np.sqrt(np.mean(s.fun ** 2))), rms_by_C=rms,
                n_points=len(s.fun))
