import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from inputprograms.interp_density import OxidizerDatabase

# 圧力と分率の範囲設定
P_vals = np.linspace(0.1, 5.0, 490)  # MPa（横軸）
f_vals = np.linspace(0.0, 1.0, 400)  # 液相質量分率（縦軸）

rho_mix = np.zeros((len(f_vals), len(P_vals)))  # 行: 分率, 列: 圧力
P_grid, f_grid = np.meshgrid(P_vals, f_vals)

def compute_mixture_density(P_MPa: float, liquid_mass_fraction: float) -> float:
    db = OxidizerDatabase()
    rho_l = db.li_interp(P_MPa)
    rho_g = db.va_interp(P_MPa)
    f_l = liquid_mass_fraction
    f_g = 1 - f_l
    rho_mix = f_l * rho_l + f_g * rho_g
    return rho_mix

# 密度計算
for i in range(len(f_vals)):
    for j in range(len(P_vals)):
        rho_mix[i][j] = compute_mixture_density(P_vals[j], f_vals[i])

# 3Dプロット
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(P_grid, f_grid, rho_mix, cmap='viridis', edgecolor='none')
ax.set_xlabel("Pressure [MPa]")
ax.set_ylabel("Liquid Mass Fraction")
ax.set_zlabel("Mixture Density [kg/m³]")
ax.set_title("Mixture Density vs Pressure and Liquid Fraction")
plt.tight_layout()
plt.show()

# CSV出力（行: 分率, 列: 圧力）
df = pd.DataFrame(rho_mix, index=f_vals, columns=P_vals)
df.index.name = "LiquidFraction"
df.columns.name = "Pressure[MPa]"
df.to_csv("mixture_density_matrix.csv")