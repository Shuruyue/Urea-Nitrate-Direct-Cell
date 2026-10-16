"""Urea–nitrate direct cell (UOR ‖ NO3RR): thermodynamics, kinetics and mass-transport model.

Modules
-------
model        熱力學、單電極動力學、全電池求解、敏感度分析、批次模擬
yao2025      Yao et al., Sustainability 17, 3835 (2025) 圖中數位化的數據（CC BY 4.0）
calibration  以 yao2025 的 LSV 擬合模型參數
validation   47 項自動檢查（run_all）
summary      主要結果（compute）
plotting     繪圖
"""
__version__ = "0.1.0"
