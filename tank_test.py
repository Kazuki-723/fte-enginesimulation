import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from scipy.optimize import brentq

# liquid csv load
df = pd.read_csv("inputdatas\\N2O_liquid.csv")
li_pressure = df['Pres'].values
li_density = df['Dens'].values
li_temperture = df['Temp'].values

# 補完関数の作成（線形補完またはスプライン補完）
li_interp = interp1d(li_pressure, li_density, kind='cubic', fill_value='extrapolate')

# 補完用の圧力範囲を生成
li_pressure_new = np.linspace(min(li_pressure), max(li_pressure), 500)
li_density_new = li_interp(li_pressure_new)

# vapor csv load
df = pd.read_csv("inputdatas\\N2O_vapor.csv")
va_pressure = df['Pres'].values
va_density = df['Dens'].values
va_temperture = df['Temp'].values

# 補完関数の作成（線形補完またはスプライン補完）
va_interp = interp1d(va_pressure, va_density, kind='cubic', fill_value='extrapolate')

# 補完用の圧力範囲を生成
va_pressure_new = np.linspace(min(va_pressure), max(va_pressure), 500)
va_density_new = va_interp(va_pressure_new)

# 定数系ロード
tank_v = 0.002 # m^3
Pt_init = 4.0 #MPa
Pc_init = 2.0 #MPa

Kstar = 4.480801535925207e-06
dt = 0.001      # 時間ステップ [s]
t_max = 10.0   # 最大シミュレーション時間 [s]

# --- 初期質量計算 ---
rho_l = li_interp(Pt_init)
m_all = rho_l * tank_v
m_liquid = m_all
m_gas = 0.0
print(m_liquid)

# --- ログ用リスト ---
time_log = []
pressure_log = []
m_liquid_log = []
v_liquid_log = []
m_gas_log = []

Pt = Pt_init
t = 0.0
tol = 1e-6
max_iter = 1000
dp = 1e-4

# 初期値
time_log.append(t)
pressure_log.append(Pt)
m_liquid_log.append(m_liquid)
m_gas_log.append(m_gas)
print("init rho = ", rho_l)

while t < t_max and m_liquid > 0 and Pt - Pc_init > tol:
    # 流出質量計算
    rho_l =li_interp(Pt)
    mdot = Kstar * np.sqrt(2 * rho_l * (Pt - Pc_init) * 1e6)
    dm = mdot * dt
    print("dm = ", dm)

    # 液相の質量減少
    m_liquid -= dm
    m_remain = m_liquid + m_gas
    rho_average = m_remain / tank_v
    print("remain mass = ", m_remain)
    print("average density =", rho_average)

    # 排出分を空白相として三相分の体積算出
    v_empty = dm / li_interp(Pt)
    v_liquid = m_liquid / li_interp(Pt)
    v_vapor = m_gas / va_interp(Pt)
    print("liquid volume = ", v_liquid)
    print("vapor volume = ", v_vapor)
    print("empty volume = ", v_empty)

    # empty, vaporを気相、liquidを液相として系の平均密度が前のものと一致する密度を探査
    # Ptをdpずつ落として最初に符号が逆転するところで検知
    i = 1
    Pt_calc = Pt
    rho_guess = li_interp(Pt_calc) * (v_liquid/tank_v) + va_interp(Pt_calc) * ((v_vapor + v_empty)/tank_v)
    if rho_guess > rho_average:
        while i < max_iter:
            Pt_calc -= i * dp
            rho_guess = li_interp(Pt_calc) * (v_liquid/tank_v) + va_interp(Pt_calc) * ((v_vapor + v_empty)/tank_v)
            print("rho_guess = ",rho_guess)
            print("current Pt = ", Pt_calc)
            if rho_guess < rho_average:
                print("break")
                break
            Pt_calc = Pt
            Pt_calc += i * dp
            rho_guess = li_interp(Pt_calc) * (v_liquid/tank_v) + va_interp(Pt_calc) * ((v_vapor + v_empty)/tank_v)
            print("rho_guess = ",rho_guess)
            print("current Pt = ", Pt_calc)
            if rho_guess < rho_average:
                print("break")
                break
            Pt_calc = Pt
            i += 1
    else:
         while i < max_iter:
            Pt_calc -= i * dp
            rho_guess = li_interp(Pt_calc) * (v_liquid/tank_v) + va_interp(Pt_calc) * ((v_vapor + v_empty)/tank_v)
            print("rho_guess = ",rho_guess)
            print("current Pt = ", Pt_calc)
            if rho_guess > rho_average:
                print("break")
                break
            Pt_calc = Pt
            Pt_calc += i * dp
            rho_guess = li_interp(Pt_calc) * (v_liquid/tank_v) + va_interp(Pt_calc) * ((v_vapor + v_empty)/tank_v)
            print("rho_guess = ",rho_guess)
            print("current Pt = ", Pt_calc)
            if rho_guess > rho_average:
                print("break")
                break
            Pt_calc = Pt
            i += 1

    Pt = Pt_calc
    m_gas = (v_vapor + v_empty) * va_interp(Pt)
    m_liquid = v_liquid * li_interp(Pt)
    t += dt
    print("Pt = ",Pt)

    # ログ記録
    time_log.append(t)
    pressure_log.append(Pt)
    m_liquid_log.append(m_liquid)
    m_gas_log.append(m_gas)

# --- 結果プロット ---
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.plot(time_log, pressure_log)
plt.xlabel('Time [s]')
plt.ylabel('Tank Pressure [Pa]')
plt.title('Pressure Evolution')

plt.subplot(1, 2, 2)
plt.plot(time_log, m_liquid_log, label='Liquid Mass')
plt.plot(time_log, m_gas_log, label='Gas Mass')
plt.xlabel('Time [s]')
plt.ylabel('Mass [kg]')
plt.title('Phase Mass Evolution')
plt.legend()
plt.tight_layout()
plt.show()