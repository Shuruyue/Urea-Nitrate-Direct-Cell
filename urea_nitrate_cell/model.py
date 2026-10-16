"""
尿素–硝酸鹽直接電池（UOR || NO3RR）簡化模型 —— 核心計算（不含繪圖）
S1 參數由 Yao 2025 原始圖數據擬合（含反應級數；HER 由 0 mM 曲線擬合）

這個檔案只做四件事，每件事都可以用手算驗證（見 validation.py）：
  1. 熱力學：由標準生成自由能算出各半反應的標準電位
  2. 單一電極：給電位 → 算電流（Tafel 動力學 + 擴散極限）
  3. 全電池：給電流 → 反解兩極電位 → 電池電壓
  4. 批次處理：硝酸鹽濃度隨時間下降（質量守恆 + 法拉第定律）

單位一律：電位 V vs RHE（pH 14）、電流密度 mA cm-2、濃度 mM、擴散層 µm、電阻 Ω cm2
符號：電池電壓 V_cell = E_c − E_a − jR_Ω；> 0 代表自發產電，< 0 代表需外加電壓
"""
import numpy as np
from scipy.optimize import brentq

F = 96485.0           # 法拉第常數 C mol-1
RT_F = 8.314 * 298.15 / F   # 0.02569 V
D_NO3 = 1.902e-5      # NO3- 在水中的擴散係數 cm2 s-1（25 °C，CRC Handbook）
N_E_NO3 = 8           # NO3- -> NH3 的電子數
N_E_UOR = 6           # CO(NH2)2 -> N2 + CO2 的電子數

# ===========================================================================
# 1. 熱力學
# ===========================================================================
# 標準生成自由能 ΔfG°（kJ mol-1, 298.15 K；NBS 熱力學表 Wagman et al. 1982）
DG_F = {"H2O": -237.13, "OH-": -157.24, "CO2(g)": -394.36, "CO3^2-": -527.81,
        "NO3-": -111.25, "NH3(aq)": -26.50,
        # 尿素(aq) = 尿素(s) -197.33 + 溶解自由能 (ΔH_sol 16.91 kJ, ΔS_sol 62.59 J/K)
        "urea(aq)": -197.33 + (16.91 - 298.15 * 62.59e-3)}

# 標準電位（V vs RHE, pH 14）；thermo_check() 會重算並與此表比對
E0 = {
    "HER": 0.000,
    "UOR_lit": 0.37,      # 文獻常用值（Boggs 2009：-0.46 V vs SHE）
    "UOR_calc": 0.072,    # 由 ΔfG° 重算，產物 CO2(g)（Protsenko 2023 同值）
    "UOR_CO3": -0.025,    # 產物為 CO3^2-（pH 14 實際情形；Protsenko 2024 報告 EMF 0.024 V）
    "NO3RR": 0.708,       # NO3- -> NH3(aq)
    "OER": 1.229,
}


def thermo_check(verbose=False):
    """由 ΔfG° 重算標準電位，回傳 dict（V vs RHE）。
    還原電位 E° = -ΔG_red / (nF)；vs RHE = vs SHE − E°(HER, pH 14)。"""
    G, k = DG_F, 96.485   # k = F in kJ V-1 mol-1
    e_her = -(2 * G["OH-"] - 2 * G["H2O"]) / (2 * k)                    # 2H2O+2e- -> H2+2OH-
    e_uor = -((G["urea(aq)"] + 6 * G["OH-"]) - (5 * G["H2O"] + G["CO2(g)"])) / (6 * k)
    e_uor3 = -((G["urea(aq)"] + 8 * G["OH-"]) - (6 * G["H2O"] + G["CO3^2-"])) / (6 * k)
    e_no3 = -((G["NH3(aq)"] + 9 * G["OH-"]) - (G["NO3-"] + 6 * G["H2O"])) / (8 * k)
    out = dict(HER_SHE=e_her, UOR_calc=e_uor - e_her, UOR_CO3=e_uor3 - e_her,
               NO3RR=e_no3 - e_her)
    # 焓基準：CO(NH2)2(aq) + H2O(l) -> N2 + CO2(g) + 3H2，ΔH°/6F（熱中性電壓）
    dH = -393.51 - (-317.65) - (-285.83)
    out["UOR_dH"] = dH / (6 * k)
    if verbose:
        for kk, v in out.items():
            print(f"  {kk:9s} {v:+.4f} V")
    return out


def cell_E0(E0_uor):
    """標準電池電位與 ΔG°（每 mol NO3-，n = 8）。"""
    E = E0["NO3RR"] - E0_uor
    return E, -N_E_NO3 * F * E / 1e3   # V, kJ mol-1


def nernst_shift(factor, n=N_E_NO3):
    """反應物（NO3-）活性變為 factor 倍時，平衡電位的變化（V），其他活性不變。"""
    return RT_F / n * np.log(factor)


# ===========================================================================
# 2. 單一電極的電流–電位關係
# ===========================================================================
def tafel_current(eta, eta10, b):
    """Tafel 型動力學電流（mA cm-2，取正值）。
    j = j0 (10^(η/b) − 1)，其中 j0 由「η = η10 時 j = 10 mA cm-2」反推：
        j0 = 10 / (10^(η10/b) − 1)
    η ≤ 0 時電流為 0（只算正向反應；在本研究的大過電位範圍，逆反應可忽略）。
    注意：j0 只是換算用的常數，不是可量測的交換電流密度。"""
    if eta <= 0:
        return 0.0
    return 10.0 * (10 ** (eta / b) - 1) / (10 ** (eta10 / b) - 1)


def j_lim_no3(C_mM, delta_um):
    """NO3- 擴散極限電流密度（mA cm-2）：j_lim = nFDC/δ。"""
    C = C_mM * 1e-6            # mol cm-3
    delta = delta_um * 1e-4    # cm
    return N_E_NO3 * F * D_NO3 * C / delta * 1e3


def j_uor(E, p):
    if not p.get("UOR_on", True):        # 對照組：陽極沒有尿素
        return 0.0
    return tafel_current(E - p["E0_UOR"], p["E10_UOR"] - p["E0_UOR"], p["b_UOR"])


def j_oer(E, p):
    return tafel_current(E - E0["OER"], p["E10_OER"] - E0["OER"], p["b_OER"])


def j_no3(E, p, C):
    """NO3RR 電流：動力學 × 濃度^γ，再與擴散極限串聯（Koutecký–Levich）。
    E10_NO3 是 C_ref（100 mM）時的值；γ = 反應級數。"""
    jk = tafel_current(E0["NO3RR"] - E, E0["NO3RR"] - p["E10_NO3"], p["b_NO3"])
    jk *= (C / p["C_ref"]) ** p["gamma"]
    if jk <= 0:
        return 0.0
    jl = j_lim_no3(C, p["delta_um"])
    return jk * jl / (jk + jl)


def j_her(E, p):
    return tafel_current(E0["HER"] - E, E0["HER"] - p["E10_HER"], p["b_HER"])


# ===========================================================================
# 3. 全電池
# ===========================================================================
def anode_potential(j, p):
    """解 j_UOR(E) + j_OER(E) = j。電流隨 E 單調遞增，brentq 必有唯一解。"""
    lo = min(p["E0_UOR"], E0["OER"])
    return brentq(lambda E: j_uor(E, p) + j_oer(E, p) - j, lo, 3.0, xtol=1e-12)


def cathode_potential(j, p, C):
    """解 j_NO3(E) + j_HER(E) = j。電流隨 E 降低而單調遞增。"""
    return brentq(lambda E: j_no3(E, p, C) + j_her(E, p) - j, -3.0, E0["NO3RR"], xtol=1e-12)


def solve_cell(j, p, C=None):
    """給定電流密度 j（mA cm-2），回傳 dict(V, Ea, Ec, FE_UOR, FE_NO3)。"""
    C = p["C_NO3"] if C is None else C
    Ea = anode_potential(j, p)
    Ec = cathode_potential(j, p, C)
    V = Ec - Ea - j * 1e-3 * p["R_ohm"]
    return dict(V=V, Ea=Ea, Ec=Ec, FE_UOR=j_uor(Ea, p) / j, FE_NO3=j_no3(Ec, p, C) / j)


# ===========================================================================
# 4. 參數（每一個都標明來源）
# ===========================================================================
# 來源代號：
#   [Y-報告]  Yao 2025 文中直接寫出的數值
#   [Y-擬合]  由 Yao 2025 Fig. 3b / 11a 的數位化數據擬合（calibration.py）
#   [文獻]    其他文獻 / 手冊
#   [假設]    本研究假設，敏感度分析中變動
S1 = dict(
    E0_UOR=E0["UOR_lit"],
    E10_UOR=1.335, b_UOR=0.0442,   # [Y-報告] ANP，1 M NaOH + 0.5 M urea（Fig. 11b, 11c）
    E10_OER=1.515, b_OER=0.060,    # E10 [Y-報告]；b [Y-擬合] 55–65 mV/dec（Fig. 11a）
    E10_NO3=0.037, b_NO3=0.099,    # [Y-擬合] ANP，100 mM NaNO3 + 1 M NaOH（Fig. 3a/3b），已扣 iR
    gamma=1.0, C_ref=100.0,        # 反應級數 [Y-擬合] 0.87 ± 0.13，取一級
    E10_HER=-0.202, b_HER=0.081,   # [Y-擬合] ANP，0 mM NaNO3（Fig. 3b 黑線）
    delta_um=43.0,                 # 擴散層 [Y-擬合]（Fig. 3b 1 mM、10 mM 曲線；攪拌條件）
    C_NO3=100.0,                   # 極化曲線比較用濃度 = 校正濃度
    R_ohm=1.0,                     # [假設] 零間隙 AEM 電池，1 M KOH
)
SCEN = {
    "S1": S1,                                                           # 現有觸媒（Yao 2025 ANP）
    "S2": dict(S1, E10_UOR=0.60, b_UOR=0.040, E10_NO3=0.30, b_NO3=0.060),   # [假設] 改良
    "S3": dict(S1, E10_UOR=0.47, b_UOR=0.040, E10_NO3=0.60, b_NO3=0.060),   # [假設] 理想
}
SCEN_LABEL = {"S1": "S1 現有觸媒（Yao 2025）", "S2": "S2 改良情境（假設）",
              "S3": "S3 理想情境（假設）"}
# 對照組：陰極相同，陽極沒有尿素（只有 OER）
OER_ONLY = dict(S1, UOR_on=False)
# 其他實測觸媒點（操作模式地圖用）
NF_POINT = dict(E10_UOR=1.423,    # [Y-報告] Ni 泡沫 UOR
                E10_NO3=-0.24)    # [Y-擬合] Ni 泡沫 NO3RR，100 mM（Fig. 3a 黑線，含 iR 約 0.02 V）


# ===========================================================================
# 5. 敏感度分析（一次改一個參數）
# ===========================================================================
J_SENS = 20.0
SENS = [
    # (key, 低值, 高值)
    ("E10_UOR", 1.30, 1.42),
    ("E10_NO3", -0.06, 0.14),
    ("b_NO3", 0.070, 0.130),
    ("b_UOR", 0.023, 0.080),
    ("R_ohm", 0.5, 5.0),
    ("C_NO3", 10.0, 200.0),
    ("gamma", 0.5, 1.0),
    ("delta_um", 10.0, 100.0),
    ("E10_HER", -0.30, -0.10),
    ("E0_UOR", E0["UOR_CO3"], E0["UOR_lit"]),
]


def sensitivity(base=None, j=J_SENS):
    base = base or S1
    V0 = solve_cell(j, base)["V"]
    rows = []
    for key, lo, hi in SENS:
        dv = [solve_cell(j, dict(base, **{key: v}))["V"] - V0 for v in (lo, hi)]
        rows.append((key, lo, hi, dv[0], dv[1]))
    rows.sort(key=lambda r: abs(r[4] - r[3]))
    return V0, rows


# ===========================================================================
# 6. 批次處理（陰極室硝酸鹽濃度隨時間下降）
#    dC/dt = − FE_NO3 · I / (8 F V)    （法拉第定律；陽極尿素視為過量）
# ===========================================================================
BATCH = dict(V_L=0.10, A_cm2=10.0, C0=60 / 14.007, dt=10.0, removal=0.90)  # C0 = 60 mg N/L


def simulate_batch(p, strategy, j_set=10.0, frac=0.8, **kw):
    """strategy = 'cc'：恆電流 j_set；'follow'：j = min(j_set, frac × j_lim(C))。
    以顯式 Euler 積分到去除率達 removal。能量 E_J 以產電為正。"""
    b = dict(BATCH, **kw)
    C, t = b["C0"], 0.0
    T, Cs, FEs, Vs, Js = [0.0], [C], [], [], []
    E_J = Q = Q_no3 = Q_uor = 0.0
    while C > (1 - b["removal"]) * b["C0"] and t < 72 * 3600:
        j = j_set if strategy == "cc" else min(j_set, frac * j_lim_no3(C, p["delta_um"]))
        r = solve_cell(j, p, C)
        I = j * 1e-3 * b["A_cm2"]                                   # A
        dC = r["FE_NO3"] * I / (N_E_NO3 * F) / b["V_L"] * 1e3 * b["dt"]   # mM
        C = max(C - dC, 0.0)
        t += b["dt"]
        E_J += r["V"] * I * b["dt"]
        Q += I * b["dt"]
        Q_no3 += r["FE_NO3"] * I * b["dt"]
        Q_uor += r["FE_UOR"] * I * b["dt"]
        T.append(t / 3600); Cs.append(C); FEs.append(r["FE_NO3"]); Vs.append(r["V"]); Js.append(j)
    n_N_no3 = (b["C0"] - C) * 1e-3 * b["V_L"]          # mol N（硝酸鹽）
    n_N_urea = Q_uor / (N_E_UOR * F) * 2               # mol N（尿素：每分子 2 個 N）
    kWh = E_J / 3.6e6
    return dict(t=np.array(T), C=np.array(Cs), FE=np.array(FEs), V=np.array(Vs), j=np.array(Js),
                t_end=t / 3600, E_kWh=kWh, Q=Q, Q_no3=Q_no3, n_N_no3=n_N_no3, n_N_urea=n_N_urea,
                FE_avg=Q_no3 / Q,
                kWh_per_kgN_no3=kWh / (n_N_no3 * 14.007e-3),
                kWh_per_kgN_all=kWh / ((n_N_no3 + n_N_urea) * 14.007e-3))


if __name__ == "__main__":
    print("熱力學重算（V vs RHE）："); thermo_check(verbose=True)
    for k in ("UOR_lit", "UOR_calc", "UOR_CO3"):
        E, G = cell_E0(E0[k])
        print(f"  E°cell（UOR = {k}）= {E:+.3f} V, ΔG° = {G:.0f} kJ/mol NO3-")
    print(f"  NO3- 稀釋 1000 倍的 Nernst 位移：{nernst_shift(1e-3)*1e3:+.1f} mV")
    print("\n10 mA cm-2 下各情境：")
    for k, p in SCEN.items():
        r = solve_cell(10.0, p)
        print(f"  {k}: V={r['V']:+.3f} Ea={r['Ea']:.3f} Ec={r['Ec']:+.3f} "
              f"FE_NO3={r['FE_NO3']*100:.1f}%")
    r = solve_cell(10.0, OER_ONLY)
    print(f"  OER 對照: V={r['V']:+.3f} Ea={r['Ea']:.3f}")
