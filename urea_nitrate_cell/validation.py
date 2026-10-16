"""
模型驗證：逐項檢查 model.py，每項都有「預期值（手算或文獻）」與「容許誤差」。
用法：run_all() 回傳所有檢查結果（list of dict）；命令列請用 scripts/validate.py。
"""
import re
import numpy as np
from . import model as m
from .yao2025 import REPORTED, FIG3B
from . import calibration as cal

ROWS = []
VERBOSE = False


def check(cat, name, value, expect, tol, how, note=""):
    ok = abs(value - expect) <= tol
    ROWS.append(dict(cat=cat, name=name, value=value, expect=expect, tol=tol, ok=bool(ok),
                     how=how, note=note))
    if VERBOSE:
        print(f"[{'PASS' if ok else 'FAIL'}] {cat:6s} {name:42s} 模型 {value:>11.5g}  預期 {expect:>11.5g}  (±{tol:g})")


def flag(cat, name, ok, how, note=""):
    ROWS.append(dict(cat=cat, name=name, value=None, expect=None, tol=None, ok=bool(ok), how=how,
                     note=note))
    if VERBOSE:
        print(f"[{'PASS' if ok else 'FAIL'}] {cat:6s} {name}")


def _run():
    # ---------------------------------------------------------------------------
    # A. 化學計量：原子數與電荷守恆
    # ---------------------------------------------------------------------------
    def atoms(formula):
        """簡易化學式解析：支援括號，例如 CO(NH2)2。"""
        def parse(s):
            cnt, i = {}, 0
            while i < len(s):
                if s[i] == "(":
                    depth, j = 1, i + 1
                    while depth:
                        depth += {"(": 1, ")": -1}.get(s[j], 0); j += 1
                    inner = parse(s[i + 1:j - 1])
                    k = re.match(r"\d*", s[j:]).group(); mult = int(k) if k else 1
                    for e, n in inner.items():
                        cnt[e] = cnt.get(e, 0) + n * mult
                    i = j + len(k)
                else:
                    mm = re.match(r"([A-Z][a-z]?)(\d*)", s[i:])
                    e, k = mm.group(1), mm.group(2)
                    cnt[e] = cnt.get(e, 0) + (int(k) if k else 1)
                    i += len(mm.group(0))
            return cnt
        return parse(formula)


    def balance(lhs, rhs):
        """lhs / rhs：[(係數, 化學式, 電荷)]；回傳 (原子差 dict, 電荷差)。"""
        tot = {}
        for sgn, side in ((1, lhs), (-1, rhs)):
            for n, f, z in side:
                for e, k in atoms(f).items():
                    tot[e] = tot.get(e, 0) + sgn * n * k
                tot["charge"] = tot.get("charge", 0) + sgn * n * z
        return {k: v for k, v in tot.items() if v != 0}


    E = ("e", "", -1)
    REACTIONS = {
        "UOR：CO(NH2)2 + 6OH- → N2 + 5H2O + CO2 + 6e-":
            ([(1, "CO(NH2)2", 0), (6, "OH", -1)], [(1, "N2", 0), (5, "H2O", 0), (1, "CO2", 0), (6, "", -1)]),
        "NO3RR：NO3- + 6H2O + 8e- → NH3 + 9OH-":
            ([(1, "NO3", -1), (6, "H2O", 0), (8, "", -1)], [(1, "NH3", 0), (9, "OH", -1)]),
        "總反應：4CO(NH2)2 + 3NO3- → 4N2 + 4CO2 + 3NH3 + 2H2O + 3OH-":
            ([(4, "CO(NH2)2", 0), (3, "NO3", -1)],
             [(4, "N2", 0), (4, "CO2", 0), (3, "NH3", 0), (2, "H2O", 0), (3, "OH", -1)]),
        "UOR（碳酸根）：CO(NH2)2 + 8OH- → N2 + 6H2O + CO3^2- + 6e-":
            ([(1, "CO(NH2)2", 0), (8, "OH", -1)], [(1, "N2", 0), (6, "H2O", 0), (1, "CO3", -2), (6, "", -1)]),
    }
    for nm, (l, r) in REACTIONS.items():
        d = balance(l, r)
        flag("計量", nm, not d, "逐一計算 C、H、N、O 原子數與電荷，左右相減應為 0", str(d) if d else "")

    # ---------------------------------------------------------------------------
    # B. 熱力學
    # ---------------------------------------------------------------------------
    th = m.thermo_check()
    check("熱力學", "E°(NO3-/NH3) vs RHE", th["NO3RR"], m.E0["NO3RR"], 0.001, "ΔfG° 重算 vs 表列值")
    check("熱力學", "E°(UOR, CO2) vs RHE", th["UOR_calc"], 0.072, 0.001, "與 Protsenko 2023 的 0.072 V 比較")
    check("熱力學", "E°(UOR, CO3^2-) vs RHE", th["UOR_CO3"], -0.024, 0.002,
          "與 Protsenko 2024 的 EMF 0.024 V（UOR+HER）比較")
    check("熱力學", "ΔH°/6F（尿素電解）", th["UOR_dH"], 0.37, 0.01,
          "檢驗「0.37 V 可能是焓基準」的推測", "0.363 V；只是數值接近，未證實文獻來源")
    check("熱力學", "E°cell（UOR = 0.37 V）", m.cell_E0(0.37)[0], 0.338, 0.001, "0.708 − 0.37")
    check("熱力學", "ΔG°（UOR = 0.37 V），kJ/mol NO3-", m.cell_E0(0.37)[1], -8 * 96.485 * 0.338, 0.5,
          "ΔG = −nFE，n = 8")
    check("熱力學", "Nernst：NO3- 稀釋 1000 倍（mV）", m.nernst_shift(1e-3) * 1e3,
          -8.314 * 298.15 / 96485 / 8 * np.log(1000) * 1e3, 0.01, "(RT/8F) ln 1000 = 22.2 mV")

    # ---------------------------------------------------------------------------
    # C. 動力學與質傳（與手算比較）
    # ---------------------------------------------------------------------------
    for k, p in m.SCEN.items():
        check("動力學", f"{k}: j_UOR(E10) = 10 mA cm-2", m.j_uor(p["E10_UOR"], p), 10.0, 1e-9,
              "Tafel 式的定義：E10 時電流必為 10")
        check("動力學", f"{k}: j_NO3,k(E10) = 10 mA cm-2（C = C_ref）",
              m.tafel_current(m.E0["NO3RR"] - p["E10_NO3"], m.E0["NO3RR"] - p["E10_NO3"], p["b_NO3"]),
              10.0, 1e-9, "同上")
    # 數值斜率應等於 b
    p = m.S1
    E1, E2 = 1.40, 1.45
    slope = (E2 - E1) / np.log10(m.j_uor(E2, p) / m.j_uor(E1, p))
    check("動力學", "UOR 數值 Tafel 斜率（mV/dec）", slope * 1e3, p["b_UOR"] * 1e3, 0.01,
          "取兩點 E 與 log j 計算斜率，應等於輸入的 b")
    jl = m.j_lim_no3(100, 43.0)
    hand = 8 * 96485 * 1.902e-5 * (100e-3 / 1000) / (43e-4) * 1e3
    check("質傳", "j_lim(100 mM, 43 µm)（mA cm-2）", jl, hand, 1e-6,
          "手算 nFDC/δ：8 × 96485 × 1.902e-5 × 1e-4 / 4.3e-3")
    check("質傳", "j_lim 與濃度成正比（10 mM / 1 mM）", m.j_lim_no3(10, 43) / m.j_lim_no3(1, 43), 10, 1e-9,
          "線性關係")
    q = dict(m.S1, delta_um=43.0)
    check("質傳", "K–L：極負電位時 j_NO3 → j_lim（20 mM）", m.j_no3(-1.5, q, 20) / m.j_lim_no3(20, 43),
          1.0, 1e-3, "動力學極快時，電流由擴散決定")
    check("質傳", "K–L：j ≪ j_lim 時 j_NO3 ≈ j_k（100 mM, 0.65 V）",
          m.j_no3(0.65, q, 100) / m.tafel_current(m.E0["NO3RR"] - 0.65, m.E0["NO3RR"] - q["E10_NO3"], q["b_NO3"]),
          1.0, 1e-3, "擴散極快時，電流由動力學決定")

    # ---------------------------------------------------------------------------
    # D. 數值解法
    # ---------------------------------------------------------------------------
    js = np.linspace(0.5, 40, 80)
    worst_a = worst_c = 0.0
    mono = bounds = True
    for k, p in list(m.SCEN.items()) + [("OER", m.OER_ONLY)]:
        prev = None
        for j in js:
            r = m.solve_cell(j, p)
            worst_a = max(worst_a, abs(m.j_uor(r["Ea"], p) + m.j_oer(r["Ea"], p) - j))
            worst_c = max(worst_c, abs(m.j_no3(r["Ec"], p, p["C_NO3"]) + m.j_her(r["Ec"], p) - j))
            if prev and not (r["Ea"] > prev["Ea"] and r["Ec"] < prev["Ec"]):
                mono = False
            if not (0 <= r["FE_NO3"] <= 1 + 1e-9 and 0 <= r["FE_UOR"] <= 1 + 1e-9
                    and r["Ec"] <= m.E0["NO3RR"] and r["V"] < m.cell_E0(p["E0_UOR"])[0]):
                bounds = False
            prev = r
    check("數值", "陽極電流平衡最大殘差（mA cm-2）", worst_a, 0, 1e-8, "代回方程式 j_UOR + j_OER − j")
    check("數值", "陰極電流平衡最大殘差（mA cm-2）", worst_c, 0, 1e-8, "代回方程式 j_NO3 + j_HER − j")
    flag("數值", "電流增加時 E_a 單調上升、E_c 單調下降（4 組 × 80 點）", mono, "極化曲線不可回折")
    flag("數值", "0 ≤ FE ≤ 1、E_c ≤ E°(NO3RR)、V_cell < E°cell（全部點）", bounds, "物理上限")
    # 與解析反函數比較（S3、10 mA，副反應可忽略時）
    p = m.SCEN["S3"]
    eta10 = p["E10_UOR"] - p["E0_UOR"]
    Ea_exact = p["E0_UOR"] + p["b_UOR"] * np.log10(10 / 10 * (10 ** (eta10 / p["b_UOR"]) - 1) + 1)
    check("數值", "S3 陽極電位：brentq vs 解析反函數（V）", m.anode_potential(10, p), Ea_exact, 1e-6,
          "Tafel 式可直接反解：E = E° + b log10[(j/j0) + 1]")

    # ---------------------------------------------------------------------------
    # E. 與 Yao 2025 數據比較（校正與獨立檢驗）
    # ---------------------------------------------------------------------------
    an = cal.fit_anode()
    check("校正", "陽極 LSV：模型參數重現 Fig. 11a（RMS 相對誤差）", an["used"]["rms"], 0, 0.10,
          "UOR + OER 兩條曲線、2–175 mA cm-2，只擬合 R_u", f"R_u = {an['used']['Ru']:.2f} Ω cm2")
    check("校正", "陽極 iR 校正後 E10(UOR)（V）", an["fixed_b"]["E10_UOR"], REPORTED["UOR_E10_ANP"], 0.01,
          "b 固定 44.2 mV/dec 擬合 vs 論文報告 1.335 V")
    check("校正", "UOR 與 OER 在 10 mA 的電位差（mV）",
          (m.anode_potential(10, m.OER_ONLY) - m.anode_potential(10, m.S1)) * 1e3, 180, 5,
          "Fig. 11a 標註 180 mV")
    ca = cal.fit_cathode(gamma=1.0)
    check("校正", "陰極 LSV：0/1/10/100 mM 擬合 RMS", ca["rms_fit"], 0, 0.15,
          "4 條曲線共 127 點、7 個參數", "")
    check("校正", "陰極 LSV：100 mM 曲線 RMS", ca["rms_by_C"][100], 0, 0.10, "模型主要使用的濃度")
    check("檢驗", "陰極 LSV：50 mM 曲線（未參與擬合）RMS", ca["rms_by_C"][50], 0, 0.35,
          "獨立預測", "Yao 的 50 mM 曲線在 −0.85 ~ −1.0 V 高於 100 mM，與濃度單調性矛盾，屬原始數據問題")
    check("校正", "擬合 E10(NO3RR) 與模型採用值（V）", ca["params"]["E10_NO3"], m.S1["E10_NO3"], 0.002, "")
    check("校正", "擬合 b(NO3RR) 與模型採用值（V/dec）", ca["params"]["b_NO3"], m.S1["b_NO3"], 0.001, "")
    check("校正", "擬合 E10(HER) 與模型採用值（V）", ca["params"]["E10_HER"], m.S1["E10_HER"], 0.002, "")
    check("校正", "擬合 δ 與模型採用值（µm）", ca["params"]["delta_um"], m.S1["delta_um"], 1.0, "")
    cf = cal.fit_cathode()
    check("校正", "反應級數 γ（放開擬合）", cf["params"]["gamma"], 1.0, 2 * cf["std"]["gamma"],
          "γ = 1 是否在 2 個標準誤差內", f"γ = {cf['params']['gamma']:.2f} ± {cf['std']['gamma']:.2f}")

    # ---------------------------------------------------------------------------
    # F. 批次模擬：質量守恆、時間步長、能量符號
    # ---------------------------------------------------------------------------
    r = m.simulate_batch(m.S1, "cc")
    check("批次", "質量守恆：去除的 NO3- × 8F = 用於 NO3RR 的電量", r["n_N_no3"] * 8 * m.F / r["Q_no3"], 1.0,
          1e-9, "兩邊獨立累加")
    flag("批次", "濃度隨時間單調下降", bool(np.all(np.diff(r["C"]) <= 0)), "")
    # FE = 1 時的解析時間：t = 8F ΔC V / I
    p_ideal = dict(m.S1, delta_um=1e-3, E10_HER=-5.0)       # 讓 FE ≡ 1
    ri = m.simulate_batch(p_ideal, "cc")
    t_hand = 8 * m.F * (0.9 * m.BATCH["C0"] * 1e-3 * m.BATCH["V_L"]) / (10e-3 * m.BATCH["A_cm2"]) / 3600
    check("批次", "FE = 1 時的處理時間（h）", ri["t_end"], t_hand, m.BATCH["dt"] / 3600,
          "法拉第定律手算 t = nFΔCV / I")
    r2 = m.simulate_batch(m.S1, "cc", dt=2.0)
    check("批次", "時間步長 10 s vs 2 s：處理時間相對差", abs(r["t_end"] / r2["t_end"] - 1), 0, 0.005, "收斂性")
    check("批次", "時間步長 10 s vs 2 s：能耗相對差", abs(r["E_kWh"] / r2["E_kWh"] - 1), 0, 0.005, "收斂性")
    rf = m.simulate_batch(m.S1, "follow")
    r2f = m.simulate_batch(m.S1, "follow", dt=2.0)
    check("批次", "時間步長（跟隨策略）：能耗相對差", abs(rf["E_kWh"] / r2f["E_kWh"] - 1), 0, 0.005, "收斂性")
    flag("批次", "S1 全程 V_cell < 0（耗電），能量為負", bool(np.all(r["V"] < 0) and r["E_kWh"] < 0),
         "符號慣例：產電為正")
    check("批次", "S3 於 10 mA cm-2 為正電壓（產電）", float(m.solve_cell(10, m.SCEN["S3"])["V"] > 0), 1, 0, "")

    # ---------------------------------------------------------------------------
    # G. 操作模式地圖的直線近似
    # ---------------------------------------------------------------------------
    worst, n_in, n_out = 0.0, 0, 0
    for xa in np.linspace(0.45, 1.6, 24):
        for yc in np.linspace(-0.4, 0.68, 24):
            rr = m.solve_cell(10, dict(m.S1, E10_UOR=xa, E10_NO3=yc))
            if rr["FE_NO3"] >= 0.9 and rr["FE_UOR"] >= 0.9:
                n_in += 1
                worst = max(worst, abs(rr["V"] - (yc - xa - 0.01)))
            else:
                n_out += 1
    check("地圖", "主反應區（兩極 FE ≥ 90 %）：直線 V = y − x − jR 與完整模型最大差（mV）", worst * 1e3, 0, 5,
          f"576 個格點中 {n_in} 點位於主反應區", f"其餘 {n_out} 點副反應（HER 或 OER）佔 >10 % 電流，"
          "直線失效；因此地圖改用完整模型計算等電壓線，並標出副反應區")


def run_all(verbose=False):
    """執行全部檢查，回傳 list of dict（cat, name, value, expect, tol, ok, how, note）。"""
    global VERBOSE
    VERBOSE = verbose
    ROWS.clear()
    _run()
    return list(ROWS)
