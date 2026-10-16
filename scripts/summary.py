"""印出主要結果並寫入 outputs/results.json。"""
import json
import os
import _path  # noqa: F401
from urea_nitrate_cell.summary import compute

if __name__ == "__main__":
    N = compute()
    t, op, B = N["thermo"], N["op_10mA"], N["batch_60mgN_per_L"]
    print(f"E°cell = {t['Ecell_UOR_lit']:+.2f} / {t['Ecell_UOR_calc']:+.2f} / {t['Ecell_UOR_CO3']:+.2f} V"
          f"（UOR E° = 0.37 / 0.072 / −0.025 V）")
    print("10 mA cm-2, 100 mM NO3-:")
    for k, r in op.items():
        print(f"  {k:9s} V_cell = {r['V']:+.2f} V  (E_a = {r['Ea']:.3f}, E_c = {r['Ec']:+.3f})")
    print(f"  urea saves {N['urea_saving']['dV']:.2f} V ({N['urea_saving']['pct']:.0f} %); "
          f"S3 peak {N['S3_peak_power']['P_mW_cm2']:.1f} mW cm-2 at {N['S3_peak_power']['j']:.0f} mA cm-2")
    print("Batch, 60 mg N/L → 90 % removed:")
    for k, b in B.items():
        print(f"  {k:16s} t = {b['t_h']:.2f} h, FE = {b['FE_avg']*100:.0f} %, "
              f"{b['kWh_per_kgN_nitrate']:+.1f} kWh/kg NO3-N, {b['kWh_per_kgN_total']:+.1f} kWh/kg N")
    os.makedirs(_path.OUTPUTS, exist_ok=True)
    json.dump(N, open(os.path.join(_path.OUTPUTS, "results.json"), "w"), ensure_ascii=False, indent=1, default=float)
