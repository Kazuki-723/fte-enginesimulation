import numpy as np
import math
import matplotlib.pyplot as plt

# シミュレーションパラメータ
P0 = 1.697e6             # 初期内圧 [Pa]
Pa = 101325              # 大気圧 [Pa]
Df = 45.86               # 最終ポート径 [mm]
Vt = 2000 * 1e3          # タンク容積 [mm^3], 累乗はccからの補正
# 燃焼室容積：ポート径 Df により、直径の2乗×長さ（360 mm）で計算（mm^3単位）
V_mm3 = (Df ** 2 * 360 * math.pi) / 4  
V = (V_mm3 + Vt) * 1e-9  # m^3に変換
T0 = 3096.55             # 温度 [K]（定常温度と仮定）
gamma = 1.2606           # 比熱比
Mole = 28.442
R = 8314 / Mole          # 気体定数 [J/(kg*K)]

# ユーザー指定の初期推力 [N]
F_init = 308.2

# ノズル幾何学
# 燃焼室断面積 Ac：ポート径Dfを利用（m単位）
Ac = (math.pi * (Df / 1000)**2) / 4    # [m^2]
At = Ac / 3.0                          # 収束部（ノズル喉部）面積 [m^2]（面積比で3分の1）
Ae = 3.16 * At                         # 拡散部（出口）面積 [m^2]

# ノズルの質量流率に含む定数 F_const（チョークド流の条件）
F_const = np.sqrt(gamma / R) * (2 / (gamma + 1))**((gamma + 1) / (2 * (gamma - 1)))

# 初期ノズル出口速度（等エントロピー膨張の理論式）
v_e0 = math.sqrt((2 * gamma / (gamma - 1)) * R * T0 *
                 (1 - (Pa / P0)**((gamma - 1) / gamma)))

# 初期状態における質量流率は
# mdot0 = C_d * At * (P0/sqrt(T0)) * F_const
# 初期推力 F_init = mdot0 * v_e0 から、C_d を求める
Cd = F_init / (At * (P0 / math.sqrt(T0)) * F_const * v_e0)
print(f"Calculated Cd (for continuous initial thrust {F_init} N): {Cd:.4f}")

# 初期質量（室内の理想気体： m = P * V / (R T)）
m0 = P0 * V / (R * T0)

# シミュレーションの設定
dt = 1e-5            # タイムステップ [s]
t_max = 0.1          # 最大シミュレーション時間 [s]
N_steps = int(t_max / dt)

# 状態保存用リスト
time_arr = []
P_arr = []
thrust_arr = []
impulse_arr = []

# 初期状態の設定
P = P0               # 内圧 [Pa]
m = m0               # ガス質量 [kg]
t = 0.0
cumulative_impulse = 0.0

# シミュレーションループ：内圧が大気圧以下になるまでシミュレーションする
for i in range(N_steps):
    if P <= Pa:
        break

    # チョークド流の仮定下での質量流率計算
    mdot = Cd * At * P / np.sqrt(T0) * F_const   # [kg/s]

    # ノズル出口速度 v_e の計算（等エントロピー膨張条件）
    v_e = math.sqrt((2 * gamma / (gamma - 1)) * R * T0 *
                    (1 - (Pa / P)**((gamma - 1) / gamma)))

    # 推力の計算： F = mdot * v_e (出口圧力補正項は、理想膨張ならゼロに近い)
    thrust = mdot * v_e  # [N]

    # 積分してインパルス（推力積分値）を更新
    cumulative_impulse += thrust * dt  # [N·s]

    # 内部ガス質量の更新（オイラー積分）
    m = m - mdot * dt
    if m <= 0:
        break

    # 内部圧力の更新： 定常温度を仮定して P = m R T / V
    P = m * R * T0 / V

    # 結果の記録
    time_arr.append(t)
    P_arr.append(P)
    thrust_arr.append(thrust)
    impulse_arr.append(cumulative_impulse)

    t += dt

# 結果表示
print(f"Simulation end time: {t:.6f} s")
print(f"Final chamber pressure: {P:.2f} Pa")
print(f"Total impulse: {cumulative_impulse:.4f} N·s")
print(f"Initial mass: {m0:.6f} kg, Final mass: {m:.6f} kg")

# プロット
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(time_arr, np.array(P_arr) / 1e6)  # 圧力をMPa表示
plt.xlabel("Time [s]")
plt.ylabel("Chamber Pressure [MPa]")
plt.title("Time vs Chamber Pressure")

plt.subplot(1, 2, 2)
plt.plot(time_arr, thrust_arr)
plt.xlabel("Time [s]")
plt.ylabel("Thrust [N]")
plt.title("Time vs Thrust")

plt.tight_layout()
plt.show()
