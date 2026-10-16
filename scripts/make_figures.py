"""產生全部的圖到 docs/figures/（或第一個參數指定的目錄）。"""
import os
import sys
import _path  # noqa: F401
from urea_nitrate_cell import plotting

if __name__ == "__main__":
    plotting.OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_path.ROOT, "docs", "figures")
    print(plotting.make_all())
