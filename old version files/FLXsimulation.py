from pythonfiles import NASACEA
import numpy as np
import matplotlib.pyplot as plt
import math
import csv

#------------------------#
# 不変定数の定義
#------------------------#

R_univ = 8314  # 一般気体定数 [J/mol-K]

Pe_init = 0.1113  # 出口圧の定義（初期状態での適正膨張・過膨張を評価するため） [MPa]

F_req = 650  # 初期推力 [N]

Pa = 0.1013  # 地上雰囲気圧力の定義
Pa_max = 0.1013
Pa_min = 0.0950  # 作動終了時高度の定義


# ---設計変数定義---#
# タンク内圧
# 梟の燃焼試験時はstem部にいれた圧力計が充填時4.6MPa,カウント0で4MPaくらいまで下がってstemが抜けてた
Ptank_init = 4  # 初期タンク内圧[MPa]
Ptank_fin = 2.3  # 作動終了時タンク内圧[MPa]
# 燃焼室圧力
Pc_init = 2.0  # Pc_def
Pc_fin = 1.3

# 燃料の性質
rho_ox_init = 883  # 酸化剤密度[kg/m^3]（圧力依存性を実装したい）
rho_f_start = 1190  # 燃料密度[kg/m^3]（const.）
a_ox = 0.000131  # 酸化剤流束係数（PMMA, SUBARU）
n_ox = 0.34  # 酸化剤流束指数（PMMA, SUBARU）

Vol_ox = 2.000  # [L]

# グレイン形状
Df_init = 0.034  # 燃料内径[m](（入力, const）
Lf = 0.4  # 燃料長さ[m]（とりあえず出力）

# 性能値
eta_cstar = 0.80
eta_nozzle = (1 + math.cos(math.pi / 12)) / 2
eta = eta_cstar * eta_nozzle
F_req = F_req * eta

#------------------------#
# CEAの仮パラメータの入力
# Pc, OF, epsilon
#------------------------#

Pc_def = 2.0  # 初期燃焼室圧 [MPa](初期状態の収束計算においてconst. 設計変数)

OF_def = 6.5  # 仮OF [-]
# 酸化剤流量に対してK*（discharge coef * crosssectional area orifice）が定まる

epsilon_start = 3.3  # 仮開口比 [-]
epsilon_new = epsilon_start

# 質量流量
mdot_new = 0.32
mdot_old = mdot_new

# 仮入力でCEAを回し、出力を得る（気体の性質はノズル出口において）
gamma_tmp1, Cstar_tmp1, CF_tmp1, T_c_tmp1, T_t_tmp1, T_e_tmp1, Mole_tmp1, Pthroat_tmp1, Pe_tmp1, Mach_tmp1 = NASACEA.CEA(Pc_def, OF_def, epsilon_start)

# 収束用パラメータ
Pe_old = Pe_tmp1
diff_exit = 2
i = 0

# グラフ用
Pc_arr = np.zeros(0)
Pc_arr = np.append(Pc_arr, Pc_def)
Pe_arr = np.zeros(0)
Pe_arr = np.append(Pe_arr, Pe_tmp1)
eps_arr = np.zeros(0)
eps_arr = np.append(eps_arr, epsilon_start)

R_tmp1 = R_univ / Mole_tmp1  # ガス定数
a_tmp1 = np.sqrt(gamma_tmp1 * R_tmp1 * T_e_tmp1)  # 音速

Ve_tmp1 = a_tmp1 * Mach_tmp1
F = mdot_new * Ve_tmp1
diff_F = F_req - F

#--------------------------#
# 初期条件の収束
#--------------------------#

while abs(diff_F) > 0.1:

    if diff_F >= 0.1:
        mdot_new = mdot_old + 0.0001
    else:
        mdot_new = mdot_old - 0.0001

    mdot_old = mdot_new

    # CEA計算
    gamma_tmp1, Cstar_tmp1, CF_tmp1, T_c_tmp1, T_t_tmp1, T_e_tmp1, Mole_tmp1, Pthroat_tmp1, Pe_tmp1, Mach_tmp1 = \
        NASACEA.CEA(Pc_def, OF_def, epsilon_new)

    #------------------------#
    # epsilonの計算
    #------------------------#

    # 燃焼条件に依存する気体の性質定義（@ノズル出口）
    R_tmp1 = R_univ / Mole_tmp1  # ガス定数
    a_tmp1 = np.sqrt(gamma_tmp1 * R_tmp1 * T_e_tmp1)  # 音速

    # 出口速度の計算
    Ve_tmp1 = a_tmp1 * Mach_tmp1

    # スロース面積の計算
    At_new = (mdot_new * np.sqrt((R_tmp1 * T_t_tmp1) / gamma_tmp1)) / (Pthroat_tmp1 * 1000000)

    # 出口圧力が大気圧となる出口マッハ数の計算
    Me_new = np.sqrt((2 * ((Pc_def / Pe_init) ** ((gamma_tmp1 - 1) / gamma_tmp1) - 1)) / (gamma_tmp1 - 1)) + 0.01

    # 出口面積の計算
    Ae_new = At_new * (((1 + (((gamma_tmp1 - 1) * Me_new ** 2) / 2)) / ((gamma_tmp1 + 1) / 2)) ** ((gamma_tmp1 + 1) / (2 * (gamma_tmp1 - 1)))) / Me_new

    # 膨張比の計算
    epsilon_new = Ae_new / At_new

    # 推力の計算（要求推力に対する推力の差分を計算）
    # F = Ve_tmp1 * mdot_new

    # CFで推力計算
    CF_tmp1 = CF_tmp1 + (Pe_tmp1 - Pa) * epsilon_new / Pc_def
    F = CF_tmp1 * Cstar_tmp1 * eta * mdot_new
    diff_F = F_req - F
    diff_exit = Pe_tmp1 - Pe_init

    # スロース径、出口径の計算
    Dt = 2 * np.sqrt(At_new / math.pi)
    De = 2 * np.sqrt(Ae_new / math.pi)

    # 収束失敗時に脱出
    if diff_exit == 0.0:
        i = i + 1
        if i > 3:
            print("収束失敗")
            break

print("Thrust = ", F, "[N]")
print("diff_F = ", diff_F, "[N]")
print("mdot = ", mdot_new, "[kg/s]")
print("Pe = ", Pe_tmp1, "[MPa]")
print("diff_atm = ", abs(diff_exit))
print("epsilon_new = ", epsilon_new)
print("Dt, De = ", Dt, De, "[m]")
print("END")
print("-------------")

#------------------------#
# 積分計算の開始
#------------------------#

"""
# ---定数定義---#
# タンク内圧の定義
Pt_init = 3.5 # 初期タンク内圧[MPa]
Pt_fin = 2.5 # 作動終了時タンク内圧[MPa]
Pc_init = Pc_def
Pc_fin = 1.3

# 燃料の性質
rho_ox_start = 883 # 酸化剤密度[kg/m^3]
rho_f_start = 1190 # 燃料密度[kg/m^3]（const.）
a_ox = 0.000131 #酸化剤流束係数
n_ox = 0.34 #酸化剤流束指数


# グレイン形状
Df_init = 0.034 #燃料内径[m](（とりあえず出力）
Lf = 0.39 #燃料長さ[m]（とりあえず出力）
"""

mdot_start = mdot_new
mdot_ox_init = (OF_def / (OF_def + 1)) * mdot_start  # 初期酸化剤流量[kg/s]
mdot_f_init = (1 / (OF_def + 1)) * mdot_start  # 初期燃料流量[kg/s]

Kstar = 0  # 手打ち用
Kstar = mdot_ox_init / np.sqrt(2 * rho_ox_init * ((Ptank_init - Pc_init) * 1000000))
print("mdot_ox", mdot_ox_init)
print("Kstar = ", Kstar)

# 燃焼面積の計算
Ap_req = mdot_f_init / (rho_f_start * a_ox * ((4 * mdot_ox_init) / (math.pi * Df_init ** 2)) ** n_ox)  # 定義したOFを実現するのに必要な燃焼面積

Lf = Ap_req / (Df_init * math.pi)
print("Lf = ", Lf)
print("Ap = ", Ap_req)
print("mdot_f = ", mdot_f_init)
print("O/F = ", mdot_ox_init / mdot_f_init)

print("---------------START INTEGRATION---------------")
#====================#
# 積分計算
#====================#

# 積分計算のための定数定義（上の出力に非依存）
Mass_ox = Vol_ox * rho_ox_init  # 酸化剤質量[L] * [kg/m^3] ... [g] 2000cc
Mass_ox_remain = Mass_ox  # 酸化剤残量
print("Mass_ox = ", Mass_ox)
delta_t = 0.001  # 微小時間[s]

Ptank_tmp1 = Ptank_init
Pc_tmp1 = Pc_def
Df = Df_init

# 計算ごとに再定義が必要な値
# 作動開始からの経過時間
t = 9000  # 収束前の作動時間を勘で入力[ms]：[s]だと変な値が出た。たぶん小数点以下の桁数の問題
k = 0

# memo : NASA-CEAの入力はPc, O/F, epsilon
# 初回入力は
# Pt_tmp1 = Pt_init
# Pc_tmp1 = Pc_init
# O/F = OF_def
# epsilon = const.
print("Ap_req = ", Ap_req)
print(Ptank_tmp1)
print(Pc_tmp1)
print("------------")

OX = Mass_ox * 1000
Pt_arr = np.zeros(0)
Pt_arr = np.append(Pt_arr, Ptank_tmp1)
Pc_arr = np.zeros(0)
Pc_arr = np.append(Pc_arr, Pc_tmp1)
F_arr = np.zeros(0)
F_arr = np.append(F_arr, F)
OF_arr = np.zeros(0)
OF_arr = np.append(OF_arr, OF_def)
Ap_arr = np.zeros(0)
Ap_arr = np.append(Ap_arr, Ap_req)
mdot_arr = np.zeros(0)
mdot = mdot_ox_init + mdot_f_init
mdot_arr = np.append(mdot_arr, mdot)
Cstar_arr = np.zeros(0)
Cstar_arr = np.append(Cstar_arr, Cstar_tmp1)
CF_arr = np.zeros(0)
CF_arr = np.append(CF_arr, CF_tmp1)
F_fte_arr = np.zeros(0)
F_fte_arr = np.append(F_fte_arr, F)
M_ox_arr = np.zeros(0)
M_ox_arr = np.append(M_ox_arr, Mass_ox_remain)

Ap = Ap_req
It = 0

Pa_tmp1 = Pa_max  # 暫定[Pa]
print("epsilon_new = ", epsilon_new)

# while(diff_t > 0.1): # 作動時間の収束
while Mass_ox_remain >= 1:  # 酸化剤残量が0になるまで計算を続ける
    delta_p = (Ptank_tmp1 - Pc_tmp1) * 1000000
    mdot_ox = (Kstar * np.sqrt(2 * rho_ox_init * delta_p))  # 微小時間における流量[g/ms]
    mdot_f = (Ap * rho_f_start * a_ox * ((4 * mdot_ox) / (math.pi * Df ** 2)) ** n_ox)  # 微小時間における燃料流量[g/ms]
    print("Df = ", Df)

    # rdotによる評価
    rdot = a_ox * ((4 * mdot_ox) / (math.pi * Df ** 2)) ** n_ox
    print(rdot)
    Df = Df + (2 * rdot / 1000)
    Ap = Df * math.pi * Lf

    print("mdot_ox = ", mdot_ox, "[g/ms]")
    print("mdot_f = ", mdot_f, "[g/ms]")

    OF_tmp1 = mdot_ox / mdot_f
    print("OF_tmp1", OF_tmp1)

    gamma_tmp1, Cstar_tmp1, CF_tmp1, T_c_tmp1, T_t_tmp1, T_e_tmp1, Mole_tmp1, Pthroat_tmp1, Pe_tmp1, Mach_tmp1 = \
        NASACEA.CEA(Pc_tmp1, OF_tmp1, epsilon_new)

    R_tmp1 = R_univ / Mole_tmp1  # ガス定数
    a_tmp1 = np.sqrt(gamma_tmp1 * R_tmp1 * T_e_tmp1)  # 音速

    # 推力の計算
    F_fte = eta * ((mdot_ox + mdot_f) * a_tmp1 * Mach_tmp1) + (Pe_tmp1 - Pa_tmp1) * Ae_new
    # CFを実装する
    CF_tmp1 = CF_tmp1 + (Pe_tmp1 - Pa) * epsilon_new / Pc_tmp1
    F_new = eta * Cstar_tmp1 * (mdot_ox + mdot_f) * CF_tmp1

    # Pcの上書き
    print("F = ", F_new)
    print("Pe = ", Pe_tmp1)
    Pt = ((mdot_ox + mdot_f) * np.sqrt((R_tmp1 * T_t_tmp1) / gamma_tmp1)) / At_new
    Pc_tmp1 = Pt * (1 + (((gamma_tmp1 - 1)) / 2)) ** (gamma_tmp1 / (gamma_tmp1 - 1)) / 1000000

    # 外気圧力計算
    Pa_tmp1 = Pa_tmp1 - ((Pa_max - Pa_min) / t)

    # print("Pc_tmp1 = ", Pc_tmp1)
    Mass_ox_remain = Mass_ox_remain - mdot_ox

    # Ptの計算
    Ptank_tmp1 = Ptank_fin + (Ptank_init - Ptank_fin) * (Mass_ox_remain / Mass_ox)

    k = k + 1

    print("Mass_ox = ", Mass_ox)
    print("Pc_tmp1 = ", Pc_tmp1)
    print("Ptank_tmp1 = ", Ptank_tmp1)
    print("Mass_ox = ", Mass_ox)
    print("Lf = ", Lf)
    print("k = ", k)
    # print(Ap)
    print("---------------")

    Pt_arr = np.append(Pt_arr, Ptank_tmp1)
    Pc_arr = np.append(Pc_arr, Pc_tmp1)
    F_arr = np.append(F_arr, F_new)
    OF_arr = np.append(OF_arr, OF_tmp1)
    Ap_arr = np.append(Ap_arr, Ap)
    mdot = mdot_ox + mdot_f
    mdot_arr = np.append(mdot_arr, mdot)
    Cstar_arr = np.append(Cstar_arr, Cstar_tmp1)
    CF_arr = np.append(CF_arr, CF_tmp1)
    F_fte_arr = np.append(F_fte_arr, F_fte)
    M_ox_arr = np.append(M_ox_arr, Mass_ox_remain)

    It = It + F_new * 0.001

# ---------------- #
#      result
# ---------------- #

print("----------RESULT----------")
print("O/F_init = ", OF_def, "[-]")
print("It = ", It, "[Ns]")
print("Lf = ", Lf * 1000, "[mm]")
print("Dt = ", Dt * 1000, "[mm]")
print("Df_init = ", Df_init * 1000, "[mm]")
print("Df_final = ", Df * 1000, "[mm]")
print("F_ave =", It * 1000 / k, "[N]")
print("F_init = ", F_req, "[N]")
print(Ae_new)

x = np.arange(len(F_arr))
plt.plot(x, F_arr, label="CF")
plt.plot(x, F_fte_arr, label="mdot*Pe")
# plt.plot(x,OX_arr)
# plt.plot(x,eps_arr)
plt.show()

filename = f"FLXresult.csv"
F_values = np.stack([F_arr, F_fte_arr, Pt_arr, Pc_arr, OF_arr, mdot_arr, Cstar_arr, CF_arr, M_ox_arr]).T
with open(filename, 'w', newline='', encoding='utf-8') as file:
    writer = csv.writer(file, quoting=csv.QUOTE_NONE)
    writer.writerows(F_values)