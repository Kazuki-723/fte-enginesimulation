import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from inputprograms.interp_density import OxidizerDatabase
db = OxidizerDatabase()

# 初期条件
V_total = 0.002    # タンク容積 [m^3]
Pt_init = 4.2      # 初期内圧 [MPa] 
m_gas = 0          # 初期気相質量 [kg]
dt = 0.001           # 微小時間ステップ [s]
mdot_out = 0.2     # 排出質量流量 [kg/s]
v_liquid = V_total # 初期液相体積 [m^3]
v_gas = 0          # 初期気相体積 [m^3]
density_ave = 0    # 平均圧力(計算用)
Pt = Pt_init       # 圧力(計算用)
liquid_frac = 1    # 体積分率

# 初期液相質量 [kg]
m_liquid = V_total * db.li_interp(Pt_init)
print(m_liquid * 1000)

# CSV読み込み（行: 分率, 列: 圧力）
df = pd.read_csv("mixture_density_matrix.csv", index_col=0)
f_vals = df.index.values.astype(float)           # 液相質量分率（縦軸）
P_vals = df.columns.values.astype(float)         # 圧力 [MPa]（横軸）
rho_vals = df.values                              # 混合密度 [kg/m³]

# 転置して補間用に整形（scipyは [x軸, y軸] = [圧力, 分率] の順）
rho_vals_T = rho_vals.T  # shape: (len(P_vals), len(f_vals))

# 補間関数作成（密度 → 分率・圧力 → 補間値）
interp_func = RegularGridInterpolator(
    (P_vals, f_vals), rho_vals_T, bounds_error=False, fill_value=None
)

def estimate_pressure_from_density_2d(f_l_input: float, rho_input: float) -> float:
    """
    液相質量分率と混合密度から圧力を2次元補間で推定する
    :param f_l_input: 液相質量分率（0〜1）
    :param rho_input: 混合密度 [kg/m³]
    :return: 推定圧力 [MPa]
    """
    # 圧力探索範囲
    P_search = np.linspace(P_vals.min(), P_vals.max(), 500)

    # 誤差最小の圧力を探索
    errors = []
    for P in P_search:
        rho_est = interp_func((P, f_l_input))
        errors.append(abs(rho_est - rho_input))

    P_best = P_search[np.argmin(errors)]
    return P_best

# 記録用リスト
time_list = []
pressure_list = []
liquid_fraction_list = []
m_liquid_list = []
m_gas_list = []
m_out_list = []
m_total_list = []
m_out_total = 0

# 時間ステップで更新
t = 0
i = 0
while liquid_frac > 0.01:
    m_out = mdot_out * dt
    m_out_total += m_out
    m_liquid = max(m_liquid - m_out, 0)
    v_out = m_out / db.li_interp(Pt)

    # 各相体積計算
    v_gas += v_out
    v_gas = v_gas * (1 + ((db.va_interp(Pt))/(db.li_interp(Pt) - db.va_interp(Pt))))
    v_liquid = V_total - v_gas
    liquid_frac = v_liquid / V_total
    print("liquid_frac =",liquid_frac)

    # 平均密度計算
    m_total = m_gas + m_liquid
    density_ave = m_total / V_total
    print("average density = ", density_ave)

    # 平衡圧力を更新
    P = estimate_pressure_from_density_2d(liquid_frac, density_ave)
    print("Pt_guess = ", P)

    # 密度計算
    rho_l = db.li_interp(P)
    rho_g = db.va_interp(P)
    m_liquid = rho_l * v_liquid
    m_gas = rho_g * v_gas

    # 記録
    time_list.append(t)
    pressure_list.append(P)
    liquid_fraction_list.append(liquid_frac)
    m_liquid_list.append(m_liquid)
    m_gas_list.append(m_gas)
    m_out_list.append(m_out_total)
    m_total_list.append(m_liquid + m_gas + m_out_total)

    print("------")
    print("Pt = ", Pt)
    print("out_mass = ", m_out_total * 1000)
    print("gas_mass = ", m_gas * 1000)
    print("liquid_mass = ", m_liquid * 1000)
    print("gas_v = ", v_gas * 1000)
    print("liquid_v = ", v_liquid * 1000)
    print("------")

    Pt = P
    t += dt
    i += 1

# 結果のプロット
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.plot(time_list, pressure_list)
plt.xlabel("Time [s]")
plt.ylabel("Pressure [MPa]")
plt.title("Tank Pressure Over Time")

plt.subplot(1, 2, 2)
plt.plot(time_list, liquid_fraction_list)
plt.xlabel("Time [s]")
plt.ylabel("Liquid Mass Fraction")
plt.title("Liquid Fraction Over Time")

plt.tight_layout()
plt.show()

# 質量履歴のプロット
plt.figure(figsize=(10, 6))
plt.plot(time_list, m_liquid_list, label="Liquid Mass [kg]")
plt.plot(time_list, m_gas_list, label="Gas Mass [kg]")
plt.plot(time_list, m_out_list, label="Mass Out [kg]")
plt.plot(time_list, m_total_list, label="Total Mass [kg]")
plt.xlabel("Time [s]")
plt.ylabel("Mass [kg]")
plt.title("Mass Evolution in Tank")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

