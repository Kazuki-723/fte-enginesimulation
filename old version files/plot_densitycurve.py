import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# CSV読み込み（行: 分率, 列: 圧力）
df = pd.read_csv("mixture_density_matrix.csv", index_col=0)
f_vals = df.index.values.astype(float)           # 液相質量分率（縦軸）
P_vals = df.columns.values.astype(float)         # 圧力 [MPa]（横軸）
rho_vals = df.values                              # 混合密度 [kg/m³]

# 等密度線を描画する密度値を指定
target_rho = 800  # 例：800 kg/m³

# メッシュグリッド作成
P_grid, f_grid = np.meshgrid(P_vals, f_vals)

# プロット
plt.figure(figsize=(8, 6))
contour = plt.contour(P_grid, f_grid, rho_vals, levels=[target_rho], colors='red')
plt.clabel(contour, fmt={target_rho: f'{target_rho} kg/m³'}, inline=True, fontsize=10)

# 背景に等高線マップを追加（任意）
plt.contourf(P_grid, f_grid, rho_vals, levels=50, cmap='viridis')
plt.colorbar(label="Mixture Density [kg/m³]")

# 軸ラベルとタイトル
plt.xlabel("Pressure [MPa]")
plt.ylabel("Liquid Mass Fraction")
plt.title(f"Isodensity Curve for ρ = {target_rho} kg/m³")
plt.grid(True)
plt.tight_layout()
plt.show()