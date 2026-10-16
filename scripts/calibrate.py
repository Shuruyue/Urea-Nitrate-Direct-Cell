"""以 Yao et al. 2025 的 LSV 擬合模型參數，印出結果並寫入 outputs/calibration_result.json。"""
import json
import os
import _path  # noqa: F401
from urea_nitrate_cell.calibration import fit_anode, fit_cathode
from urea_nitrate_cell.yao2025 import REPORTED

if __name__ == "__main__":
    an = fit_anode()
    print("=== 陽極（Fig. 11a）===")
    print(f"(A) 模型採用值（UOR 1.335 V / 44.2 mV/dec，OER 1.515 V / 60 mV/dec）："
          f"R_u = {an['used']['Ru']:.2f} Ω cm2，RMS 相對誤差 {an['used']['rms']*100:.1f} %")
    print("(B) 全部放開：", {k: round(v, 4) for k, v in an["free"].items()})
    print("    標準誤差：", {k: round(v, 4) for k, v in an["free_std"].items()})
    print("(C) b_UOR 固定 44.2：", {k: round(v, 4) for k, v in an["fixed_b"].items()})
    print(f"    論文報告：UOR {REPORTED['UOR_E10_ANP']} V，OER {REPORTED['OER_E10_ANP']} V")

    print("\n=== 陰極（Fig. 3b）===")
    ca_free = fit_cathode()
    ca_g1 = fit_cathode(gamma=1.0)
    for lab, ca in (("γ 放開", ca_free), ("γ = 1（模型採用）", ca_g1)):
        print(f"[{lab}] N = {ca['n_points']}，RMS = {ca['rms_fit']*100:.1f} %")
        for k, v in ca["params"].items():
            sd = ca["std"].get(k)
            print(f"   {k:9s} {v:+.4f}" + (f" ± {sd:.4f}" if sd is not None else "（固定）"))
        print("   各濃度 RMS：", {C: f"{v*100:.1f}%" for C, v in ca["rms_by_C"].items()},
              "（50 mM 未參與擬合）")
    os.makedirs(_path.OUTPUTS, exist_ok=True)
    json.dump(dict(anode=an, cathode_free=ca_free, cathode_g1=ca_g1),
              open(os.path.join(_path.OUTPUTS, "calibration_result.json"), "w"), indent=1, ensure_ascii=False,
              default=lambda o: float(o))
