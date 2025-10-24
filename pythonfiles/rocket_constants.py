import math

#------------------------#
# 不変定数の定義
#------------------------#

R_univ = 8314  # 一般気体定数 [J/mol-K]
Pe_init = 0.1113  # 出口圧の定義（初期状態での適正膨張・過膨張を評価するため） [MPa]

F_req = 650  # 初期推力 [N]

Pa = 0.1013  # 地上雰囲気圧力の定義
Pa_max = 0.1013
Pa_min = 0.0950  # 作動終了時高度の定義

a_ox = 0.000131  # 酸化剤流束係数（PMMA, SUBARU）
n_ox = 0.34  # 酸化剤流束指数（PMMA, SUBARU）


# ---設計変数定義---#
# タンク内圧
# 梟の燃焼試験時はstem部にいれた圧力計が充填時4.6MPa,カウント0で4MPaくらいまで下がってstemが抜けてた
Ptank_init = 4.2  # 初期タンク内圧[MPa]
Ptank_fin = 2.3  # 作動終了時タンク内圧[MPa]
# 燃焼室圧力
Pc_init = 2.0  # Pc_def
Pc_fin = 1.3

# 燃料の性質
rho_ox_init = 852.24  # 酸化剤密度[kg/m^3]（圧力依存性を実装したい）
rho_f_start = 1190  # 燃料密度[kg/m^3]（const.）

Vol_ox = 2.000  # [L]

# グレイン形状
Df_init = 0.034  # 燃料内径[m](（入力, const）
Lf = 0.4  # 燃料長さ[m]（とりあえず出力）

# 性能値
# eta_cstar = 0.80
# eta_nozzle = (1 + math.cos(math.pi / 12)) / 2
# eta = eta_cstar * eta_nozzle
# F_req = F_req * eta

# #------------------------#
# # CEAの仮パラメータの入力
# # Pc, OF, epsilon
# #------------------------#

# Pc_def = 2.0  # 初期燃焼室圧 [MPa](初期状態の収束計算においてconst. 設計変数)

# OF_def = 6.5  # 仮OF [-]
# # 酸化剤流量に対してK*（discharge coef * crosssectional area orifice）が定まる

# epsilon_start = 10  # 仮開口比 [-]
# epsilon_new = epsilon_start

# # 質量流量
# mdot_new = 0.32
