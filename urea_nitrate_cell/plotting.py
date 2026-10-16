"""
繪圖：極化曲線、操作模式地圖、電位階梯、法拉第效率、批次處理、敏感度、示意圖、參數校正。
預設版面：寬 14 cm、300 dpi、白底、襯線字型（Liberation Serif／Times 系）、數學符號 STIX；
刻度與圖例 9.5 pt、軸標題 10.5 pt（插入文件時以 14 cm 寬放置，圖中字級即實際字級）。
每個函式的 (w, h) 只決定長寬比。輸出目錄：OUT（預設 docs/figures）。
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Circle, Polygon
from . import model as m
from . import calibration as cal
from .yao2025 import FIG3B, FIG11A_UOR, FIG11A_OER

OUT = "docs/figures"
W_CM = 14.0
FS, FS_LAB = 9.5, 10.5
CARD = "white"
INK, GRAY = "#000000", "#5a5a5a"
C_S1, C_S2, C_S3 = "#e0602a", "#1f9d6a", "#2a6fd6"
plt.rcParams.update({
    "font.family": ["Liberation Serif", "DejaVu Serif"], "font.size": FS, "mathtext.fontset": "stix",
    "axes.labelsize": FS_LAB, "xtick.labelsize": FS, "ytick.labelsize": FS, "legend.fontsize": FS - 0.5,
    "axes.linewidth": 0.8, "lines.linewidth": 1.5, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 4, "ytick.major.size": 4,
    "legend.frameon": True, "legend.fancybox": False, "legend.edgecolor": INK, "legend.framealpha": 1,
    "legend.facecolor": CARD, "figure.facecolor": CARD, "axes.facecolor": CARD,
    "savefig.facecolor": CARD, "axes.unicode_minus": True,
})
J_UNIT = r"$\mathrm{mA\ cm^{-2}}$"
XJ = r"Current density $j$ (" + J_UNIT + ")"


def newfig(w, h, **kw):
    """寬度固定 W_CM，(w, h) 決定長寬比。"""
    wi = W_CM / 2.54
    return plt.figure(figsize=(wi, wi * h / w), **kw)


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(f"{OUT}/{name}.png", dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------------------
def fig_polarization(w=790, h=560):
    js = np.linspace(0.2, 30, 200)
    fig = newfig(w, h)
    ax = fig.add_axes([0.13, 0.15, 0.83, 0.80])
    ax.axhline(0, color=INK, lw=1.0, ls=(0, (4, 3)))
    sty = [("S1", "S1 current catalyst (fitted)", C_S1, "-"), ("S2", "S2 improved (assumed)", C_S2, "--"),
           ("S3", "S3 ideal (assumed)", C_S3, "-.")]
    for k, lab, c, ls in sty:
        ax.plot(js, [m.solve_cell(j, m.SCEN[k])["V"] for j in js], color=c, ls=ls, label=lab)
    ax.plot(js, [m.solve_cell(j, m.OER_ONLY)["V"] for j in js], color=GRAY, ls=":", lw=2.2,
            label="Control: S1 with OER anode")
    v1 = m.solve_cell(10, m.S1)["V"]
    ax.plot(10, v1, "o", color=C_S1, mec=INK, ms=7, zorder=5)
    ax.annotate(f"{v1:.2f} V".replace("-", "−"), (10, v1), (10, v1 + 0.28), ha="center", fontsize=FS,
                arrowprops=dict(arrowstyle="-", color=INK, lw=0.8))
    ax.text(29.3, 0.07, r"Power out ($V_\mathrm{cell} > 0$)", ha="right", va="bottom")
    ax.text(29.3, -0.07, r"Bias needed ($V_\mathrm{cell} < 0$)", ha="right", va="top")
    ax.set_xlabel(XJ)
    ax.set_ylabel(r"Cell voltage $V_\mathrm{cell}$ (V)")
    ax.set_xlim(0, 30); ax.set_ylim(-1.8, 0.45)
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.43))
    save(fig, "f_polarization")


# ---------------------------------------------------------------------------
def fig_mode_map(w=860, h=560, j=10.0):
    XL, XR, YB, YT = -0.1, 1.7, -0.45, 0.85
    xs = np.linspace(XL, XR, 91)
    ys = np.linspace(YB, YT, 66)
    V = np.zeros((len(ys), len(xs))); FEN = np.zeros_like(V); FEU = np.zeros_like(V)
    for i, yc in enumerate(ys):
        for k, xa in enumerate(xs):
            if xa <= m.E0["UOR_CO3"] + 0.01 or yc >= m.E0["NO3RR"] - 0.01:
                V[i, k] = FEN[i, k] = FEU[i, k] = np.nan
                continue
            r = m.solve_cell(j, dict(m.S1, E10_UOR=xa, E10_NO3=yc,
                                     E0_UOR=min(m.E0["UOR_lit"], xa - 0.005)))
            V[i, k], FEN[i, k], FEU[i, k] = r["V"], r["FE_NO3"], r["FE_UOR"]
    fig = newfig(w, h)
    ax = fig.add_axes([0.105, 0.14, 0.86, 0.83])
    eU, eUc, eN, Rj = m.E0["UOR_lit"], m.E0["UOR_calc"], m.E0["NO3RR"], j * 1e-3 * m.S1["R_ohm"]
    # 熱力學不可能區
    ax.axhspan(eN, YT, color="#dddddd", lw=0, zorder=0)
    ax.axvspan(XL, eUc, color="#dddddd", lw=0, zorder=0)
    # 自發區（兩種 E°UOR）
    ax.add_patch(Polygon([[eUc, eUc + Rj], [eUc, eN], [eN - Rj, eN]], closed=True, fc="#d6e6fb", lw=0, zorder=1))
    ax.add_patch(Polygon([[eU, eU + Rj], [eU, eN], [eN - Rj, eN]], closed=True, fc="#8fbcf0", lw=0, zorder=1))
    # 副反應區（FE < 90 %）
    X, Y = np.meshgrid(xs, ys)
    side = np.where((FEN < 0.9) | (FEU < 0.9), 1.0, 0.0)
    ax.contourf(X, Y, side, levels=[0.5, 1.5], colors="none", hatches=["////"], zorder=2)
    ax.contour(X, Y, side, levels=[0.5], colors=GRAY, linewidths=1.0, zorder=2)
    # 等電壓線（完整模型）
    cs = ax.contour(X, Y, V, levels=[-1.5, -1.0, -0.5], colors=GRAY, linewidths=1.2,
                    linestyles=[(0, (3, 3))] * 3, zorder=3)
    ax.clabel(cs, fmt=lambda v: f"{v:+.1f} V".replace("-", "−"), fontsize=FS, inline_spacing=4,
              manual=[(1.32, -0.19), (1.15, 0.14), (0.95, 0.44)])
    for t in cs.labelTexts:
        t.set_bbox(dict(fc=CARD, ec="none", pad=0.5))
    cs0 = ax.contour(X, Y, V, levels=[0.0], colors=INK, linewidths=1.8, zorder=3)
    ax.clabel(cs0, fmt=lambda v: r"$V_\mathrm{cell}$ = 0", fontsize=FS, manual=[(0.22, 0.23)])
    for xx, ls in ((eU, "-"), (eUc, (0, (4, 3)))):
        ax.axvline(xx, color=INK, lw=1.0, ls=ls, zorder=3)
    ax.axhline(eN, color=INK, lw=1.0, zorder=3)
    ax.text(eUc - 0.015, 0.33, "UOR E° = 0.072 V", rotation=90, ha="right", va="center", zorder=4)
    ax.text(eU - 0.015, 0.13, "UOR E° = 0.37 V", rotation=90, ha="right", va="center", zorder=4)
    ax.text(XR - 0.03, (eN + YT) / 2, "NO$_3$RR E° = 0.71 V", ha="right", va="center", zorder=4)
    ax.text(0.09, 0.62, "Spontaneous", zorder=4)
    ax.text(1.05, -0.30, "Side reaction > 10 %\n(HER or OER)", ha="center", fontsize=FS - 1, zorder=4,
            bbox=dict(fc=CARD, ec="none", pad=1.5))
    pts = [("S1", 1.335, m.S1["E10_NO3"], "S1 current (ANP)", C_S1, "o", (0.10, -0.085), "right"),
           ("S2", 0.60, 0.30, "S2 improved", C_S2, "^", (0.04, -0.02), "left"),
           ("S3", 0.47, 0.60, "S3 ideal", C_S3, "D", (0.05, 0.0), "left")]
    for k, xp, yp, lab, c, mk, (dx, dy), ha in pts:
        ax.scatter(xp, yp, s=70, color=c, marker=mk, edgecolor=INK, linewidth=0.8, zorder=6)
        ax.text(xp + dx, yp + dy, lab, ha=ha, va="center", zorder=6,
                bbox=dict(fc=CARD, ec="none", pad=1.0, alpha=0.9))
    ax.set_xlim(XL, XR); ax.set_ylim(YB, YT)
    ax.set_xlabel(r"$E_{10}$ of UOR anode (V vs RHE)")
    ax.set_ylabel(r"$E_{10}$ of NO$_3$RR cathode (V vs RHE)")
    save(fig, "f_mode_map")


# ---------------------------------------------------------------------------
def fig_ladder(w=790, h=610):
    X1 = 1.55
    fig = newfig(w, h)
    ax = fig.add_axes([0.14, 0.12, 0.84, 0.86])
    Ec10 = m.cathode_potential(10, m.S1, 100.0)
    items = [(0, m.E0["OER"], "OER  1.23", C_S1, "-"),
             (0, m.E0["NO3RR"], "NO$_3$RR  0.71", C_S3, "-"),
             (0, m.E0["UOR_lit"], "UOR  0.37 (cited)", C_S1, "-"),
             (0, m.E0["UOR_calc"], "UOR  0.07 (ΔfG°)", C_S1, "-"),
             (0, m.E0["HER"], "HER  0", C_S3, "-"),
             (X1, m.S1["E10_OER"], "OER  1.52", C_S1, "-"),
             (X1, m.S1["E10_UOR"], "UOR  1.34", C_S1, "-"),
             (X1, Ec10, f"NO$_3$RR  +{Ec10:.2f}", C_S3, "-"),
             (X1, m.S1["E10_HER"], "HER  −0.20", C_S3, (0, (3, 2)))]
    for x, y, lab, c, ls in items:
        ax.hlines(y, x - 0.2, x + 0.2, color=c, lw=3.0, ls=ls)
        dy = {"UOR  0.07 (ΔfG°)": 0.035, "HER  0": -0.035}.get(lab, 0)
        ax.text(x + 0.25, y + dy, lab, va="center")
    ax.annotate("", xy=(-0.29, m.E0["NO3RR"]), xytext=(-0.29, m.E0["UOR_lit"]),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=1.2))
    ax.text(-0.33, (m.E0["NO3RR"] + m.E0["UOR_lit"]) / 2, "+0.34 V\nspontaneous", ha="right", va="center")
    ax.annotate("", xy=(X1 - 0.29, m.S1["E10_UOR"]), xytext=(X1 - 0.29, Ec10),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=1.2))
    ax.text(X1 - 0.33, 0.97, f"−{m.S1['E10_UOR'] - Ec10:.2f} V\nbias needed", ha="right", va="center")
    ax.set_xticks([0, X1])
    ax.set_xticklabels(["Equilibrium (E°)", "At 10 mA cm$^{-2}$ (ANP, LSV)"])
    ax.tick_params(axis="x", which="both", top=False, bottom=False)
    ax.set_xlim(-1.3, X1 + 1.05); ax.set_ylim(-0.35, 1.68)
    ax.set_ylabel("Potential (V vs RHE, pH 14)")
    h1, = ax.plot([], [], color=C_S1, lw=3, label="Anode reactions")
    h2, = ax.plot([], [], color=C_S3, lw=3, label="Cathode reactions")
    ax.legend(handles=[h1, h2], loc="upper left")
    save(fig, "f_ladder")


# ---------------------------------------------------------------------------
def fig_fe(w=760, h=434):
    js = np.linspace(0.2, 40, 220)
    LBL = {10: (32, 67), 60 / 14.007: (14, 58), 1: (27, 26)}
    fig = newfig(w, h)
    ax = fig.add_axes([0.12, 0.19, 0.85, 0.77])
    for C, lab, c, ls in ((10, "10 mM", C_S3, "-."), (60 / 14.007, "4.3 mM (60 mg N L$^{-1}$)", C_S2, "--"),
                          (1, "1 mM", C_S1, "-")):
        fe = [m.solve_cell(j, m.S1, C)["FE_NO3"] * 100 for j in js]
        jl = m.j_lim_no3(C, m.S1["delta_um"])
        ax.plot(js, fe, color=c, ls=ls)
        ax.text(*LBL[C], f"{lab}\n" + r"$j_\mathrm{lim}$" + f" = {jl:.1f} " + J_UNIT, color=c, ha="center", va="center",
                bbox=dict(fc=CARD, ec="none", pad=1))
    ax.set_xlabel(XJ)
    ax.set_ylabel("NO$_3$RR Faradaic eff. (%)")
    ax.set_xlim(0, 40); ax.set_ylim(0, 105)
    save(fig, "f_fe_conc")


def fig_batch(w=760, h=434):
    runs = {k: m.simulate_batch(m.S1, k) for k in ("cc", "follow")}
    fig = newfig(w, h)
    a1 = fig.add_axes([0.12, 0.585, 0.85, 0.385])
    a2 = fig.add_axes([0.12, 0.19, 0.85, 0.385], sharex=a1)
    sty = {"cc": ("Constant 10 " + J_UNIT, C_S1, "-"),
           "follow": (r"Following: $j = \min(10,\ 0.8\,j_\mathrm{lim})$", C_S3, "--")}
    for k, r in runs.items():
        lab, c, ls = sty[k]
        a1.plot(r["t"], r["C"] / m.BATCH["C0"] * 100, color=c, ls=ls, label=lab)
        a2.plot(r["t"][1:], r["FE"] * 100, color=c, ls=ls)
    a1.axhline(10, color=GRAY, lw=1.0, ls=":")
    a1.text(0.05, 13, "90 % removed", color=GRAY, fontsize=FS - 1)
    a1.set_ylabel("$C/C_0$ (%)"); a1.set_ylim(0, 105)
    a1.legend(loc="upper right", fontsize=FS - 1)
    plt.setp(a1.get_xticklabels(), visible=False)
    a2.set_ylabel("FE (%)"); a2.set_ylim(0, 105)
    a2.set_xlabel("Time $t$ (h)")
    a2.set_xlim(0, max(r["t"][-1] for r in runs.values()) * 1.02)
    save(fig, "f_batch")
    return runs


# ---------------------------------------------------------------------------
SENS_LAB = {"E10_UOR": ("UOR $E_{10}$", "{:.2f}–{:.2f} V", 1),
            "E10_NO3": ("NO$_3$RR $E_{10}$", "{:+.2f} to {:+.2f} V", 1),
            "b_NO3": ("NO$_3$RR Tafel $b$", "{:.0f}–{:.0f} mV/dec", 1e3),
            "b_UOR": ("UOR Tafel $b$", "{:.0f}–{:.0f} mV/dec", 1e3),
            "R_ohm": (r"Ohmic $R_\Omega$", "{:.1f}–{:.0f} Ω cm²", 1),
            "C_NO3": ("NO$_3^-$ conc. $C$", "{:.0f}–{:.0f} mM", 1),
            "gamma": ("Reaction order γ", "{:.1f}–{:.1f}", 1),
            "delta_um": (r"Diffusion layer $\delta$", "{:.0f}–{:.0f} µm", 1),
            "E10_HER": ("HER $E_{10}$", "{:+.2f} to {:+.2f} V", 1),
            "E0_UOR": ("UOR E°", "{:+.3f} / {:+.2f} V", 1)}


def fig_sensitivity(w=880, h=570):
    V0, rows = m.sensitivity()
    rows = [r for r in rows if r[0] != "gamma"]
    fig = newfig(w, h)
    ax = fig.add_axes([0.36, 0.13, 0.61, 0.84])
    labs = []
    for k, (key, lo, hi, dlo, dhi) in enumerate(rows):
        ax.barh(k, dlo * 1e3, color=C_S1, height=0.62, edgecolor=INK, lw=0.6)
        ax.barh(k, dhi * 1e3, color=C_S3, height=0.62, edgecolor=INK, lw=0.6)
        if max(abs(dlo), abs(dhi)) < 1e-3:
            ax.text(4, k, "≈ 0 (no effect)", va="center", color=GRAY)
        name, fmt, sc = SENS_LAB[key]
        labs.append(f"{name}\n" + fmt.format(lo * sc, hi * sc).replace("-", "−"))
    ax.axvline(0, color=INK, lw=1.0)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(labs, fontsize=FS - 1, linespacing=1.0)
    ax.tick_params(axis="y", which="both", left=False, right=False)
    ax.set_xlabel(r"Change in $V_\mathrm{cell}$, $\Delta V$ (mV)")
    ax.set_xlim(-150, 150)
    h1 = plt.Rectangle((0, 0), 1, 1, fc=C_S1, ec=INK, lw=0.6)
    h2 = plt.Rectangle((0, 0), 1, 1, fc=C_S3, ec=INK, lw=0.6)
    ax.legend([h1, h2], ["Parameter low", "Parameter high"], loc="lower left")
    save(fig, "f_sensitivity")
    return V0, rows


# ---------------------------------------------------------------------------
def fig_schematic(w=900, h=553):
    fig = newfig(w, h)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 900); ax.set_ylim(553, 0); ax.axis("off")
    TA, TC = "#fde6d6", "#dde9fa"
    ax.add_patch(Rectangle((60, 150), 360, 260, fc=TA, ec=INK, lw=1.2))
    ax.add_patch(Rectangle((480, 150), 360, 260, fc=TC, ec=INK, lw=1.2))
    ax.add_patch(Rectangle((420, 150), 60, 260, fc="#d9d9d9", ec=INK, lw=1.0, hatch="///"))
    ax.text(450, 430, "Anion-exchange membrane (AEM)", ha="center", va="top")
    ax.add_patch(Rectangle((170, 120), 22, 250, fc=C_S1, ec=INK, lw=1.0))
    ax.add_patch(Rectangle((708, 120), 22, 250, fc=C_S3, ec=INK, lw=1.0))
    ax.text(181, 112, "Anode (Ni-based)", ha="center", va="bottom")
    ax.text(719, 112, "Cathode", ha="center", va="bottom")
    ax.plot([181, 181, 360], [82, 52, 52], color=INK, lw=1.2)
    ax.plot([540, 719, 719], [52, 52, 82], color=INK, lw=1.2)
    ax.add_patch(Rectangle((360, 28), 180, 48, fc="white", ec=INK, lw=1.2))
    ax.text(450, 52, "Load / supply", ha="center", va="center")
    for x0, x1 in ((240, 330), (610, 700)):
        ax.add_patch(FancyArrowPatch((x0, 28), (x1, 28), arrowstyle="-|>", mutation_scale=16, color=INK, lw=1.2))
        ax.text(x0 - 8, 28, "e$^-$", ha="right", va="center")
    ax.add_patch(FancyArrowPatch((560, 300), (350, 300), arrowstyle="-|>", mutation_scale=18, color=INK, lw=1.6))
    ax.text(570, 300, "OH$^-$", ha="left", va="center")
    for (x, y, r) in [(215, 250, 7), (232, 215, 5.5), (218, 190, 4.5), (240, 265, 5)]:
        ax.add_patch(Circle((x, y), r, fc="white", ec=INK, lw=0.8))
    ax.text(255, 205, "N$_2$↑", va="center")
    ax.text(255, 245, "CO$_2$ → CO$_3^{2-}$", va="center")
    ax.text(160, 190, "CO(NH$_2$)$_2$", ha="right", va="center")
    ax.text(160, 222, "urea", ha="right", va="center", color=GRAY)
    ax.text(825, 190, "NO$_3^-$", ha="right", va="center")
    ax.text(825, 222, "nitrate", ha="right", va="center", color=GRAY)
    ax.text(695, 345, "NH$_3$ (recovered)", ha="right", va="center")
    ax.add_patch(FancyArrowPatch((8, 385), (58, 385), arrowstyle="-|>", mutation_scale=14, color=INK, lw=1.2))
    ax.add_patch(FancyArrowPatch((892, 385), (842, 385), arrowstyle="-|>", mutation_scale=14, color=INK, lw=1.2))
    ax.text(70, 392, "urea wastewater + KOH", va="center", color=GRAY)
    ax.text(830, 392, "nitrate wastewater", va="center", ha="right", color=GRAY)
    ax.text(240, 490, "CO(NH$_2$)$_2$ + 6OH$^-$ → N$_2$ + CO$_2$ + 5H$_2$O + 6e$^-$", ha="center", va="center")
    ax.text(665, 490, "NO$_3^-$ + 6H$_2$O + 8e$^-$ → NH$_3$ + 9OH$^-$", ha="center", va="center")
    ax.text(240, 525, "anode (oxidation)", ha="center", va="center", color=GRAY)
    ax.text(665, 525, "cathode (reduction)", ha="center", va="center", color=GRAY)
    save(fig, "f_schematic")


# ---------------------------------------------------------------------------
def fig_calibration(w=800, h=470):
    an = cal.fit_anode()
    ca = cal.fit_cathode(gamma=1.0)
    fig = newfig(w, h)
    a1 = fig.add_axes([0.105, 0.17, 0.36, 0.71])
    a2 = fig.add_axes([0.61, 0.17, 0.36, 0.71])
    Ru = an["used"]["Ru"]
    for (E, J), uor, c, mk, lab in ((FIG11A_UOR, True, INK, "o", "UOR"), (FIG11A_OER, False, C_S1, "s", "OER")):
        a1.plot(E[::2], J[::2], mk, color=c, mfc="none", mew=1.4, ms=6, label=f"{lab} data")
        Eg = np.linspace(E.min(), E.max(), 120)
        a1.plot(Eg, [cal.anode_meas(e, m.S1, Ru, uor) for e in Eg], color=c, lw=2)
    a1.set_xlabel("E (V vs RHE)"); a1.set_ylabel("j (" + J_UNIT + ")")
    a1.set_xlim(1.30, 1.62); a1.set_ylim(-5, 185)
    a1.set_title("(a) Anode — Yao Fig. 11a", fontsize=FS)
    a1.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), fontsize=FS - 1, handlelength=1.0)
    a1.text(0.04, 0.50, "symbols: data\nlines: model", transform=a1.transAxes, fontsize=FS - 1, color=GRAY)
    q = ca["params"]
    cols = {0: GRAY, 1: C_S1, 10: C_S2, 50: "#8a5ac7", 100: C_S3}
    for C in (0, 1, 10, 50, 100):
        E, J = FIG3B[C]
        ok = np.isfinite(J)
        a2.plot(E[ok][::3], J[ok][::3], "o", color=cols[C], mfc="none" if C == 50 else cols[C],
                ms=5, mew=1.2)
        Eg = np.linspace(E[ok].min(), E[ok].max(), 100)
        a2.plot(Eg, [cal.cathode_meas(e, C, q) for e in Eg], color=cols[C], lw=2,
                ls="--" if C == 50 else "-", label=f"{C} mM")
    a2.set_xlabel("E (V vs RHE)"); a2.set_ylabel("−j (" + J_UNIT + ")")
    a2.set_xlim(-0.28, 0.30); a2.set_ylim(0, 95)
    a2.set_title("(b) Cathode — Yao Fig. 3b", fontsize=FS)
    a2.legend(loc="upper right", bbox_to_anchor=(1.0, 1.0), fontsize=FS - 1, handlelength=1.2,
              title="NO$_3^-$", title_fontsize=FS - 1)
    save(fig, "f_calibration")
    return an, ca


def make_all():
    """產生全部的圖，回傳檔名列表。"""
    fig_ladder(790, 610)
    fig_schematic(900, 553)
    fig_calibration(800, 470)
    fig_polarization(790, 540)
    fig_mode_map(860, 580)
    fig_sensitivity(880, 600)
    fig_fe(760, 450)
    fig_batch(760, 520)
    return sorted(os.listdir(OUT))
