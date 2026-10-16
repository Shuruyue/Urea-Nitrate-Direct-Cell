"""讓 scripts/ 在未安裝套件時也能直接執行（把專案根目錄加入 sys.path）。"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUTPUTS = os.path.join(ROOT, "outputs")
