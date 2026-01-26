import numpy as np
import math
import csv
import matplotlib.pyplot as plt
from pythonfiles.rocket_constants import *
from pythonfiles.cea_interface import CEAInterface

class RocketSimulation:
    def __init__(self):
        # 定数・初期パラメータのセットアップ
        self.R_univ = R_univ
        self.Pe_init = Pe_init
        self.F_req = F_req
        self.Pa = Pa
        self.Pa_max = Pa_max
        self.Pa_min = Pa_min
        self.Ptank_init = Ptank_init
        self.Ptank_fin = Ptank_fin
        self.Pc_init = Pc_init
        self.Pc_fin = Pc_fin
        self.rho_ox_init = rho_ox_init
        self.rho_f_start = rho_f_start
        self.a_ox = a_ox
        self.n_ox = n_ox
        self.Vol_ox = Vol_ox
        self.Df_init = Df_init
        self.Lf = Lf
        self.eta_cstar = eta_cstar
        self.eta_nozzle = eta_nozzle
        self.eta = eta
        self.F_req = F_req  # 更新後のF_req
        # CEA初期パラメータ
        self.Pc_def = Pc_def
        self.OF_def = OF_def
        self.epsilon_start = epsilon_start
        self.epsilon_new = epsilon_new
        
        self.mdot_new = mdot_new
        self.mdot_old = self.mdot_new
         
        # 仮入力でCEAを回し、出力を得る（気体の性質はノズル出口において）
        (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1, 
         self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1, 
         self.Pe_tmp1, self.Mach_tmp1) = CEAInterface.compute(self.Pc_def, self.OF_def, self.epsilon_start)

        # 収束用パラメータ
        self.Pe_old = self.Pe_tmp1
        self.diff_exit = 2
        self.i = 0
        
        # グラフ用
        self.Pc_arr = np.array([self.Pc_def])
        self.Pe_arr = np.array([self.Pe_tmp1])
        self.eps_arr = np.array([self.epsilon_start])
        
        self.R_tmp1 = self.R_univ / self.Mole_tmp1  # ガス定数
        self.a_tmp1 = np.sqrt(self.gamma_tmp1 * self.R_tmp1 * self.T_e_tmp1)  # 音速
        self.Ve_tmp1 = self.a_tmp1 * self.Mach_tmp1
        self.F = self.mdot_new * self.Ve_tmp1
        self.diff_F = self.F_req - self.F
        
        # 積分計算用の配列初期化
        self.Pt_arr = np.array([])
        self.Pc_int_arr = np.array([])  # 積分計算中のPcの配列
        self.F_arr = np.array([])
        self.OF_arr = np.array([])
        self.Ap_arr = np.array([])
        self.mdot_arr = np.array([])
        self.Cstar_arr = np.array([])
        self.CF_arr = np.array([])
        self.F_fte_arr = np.array([])
        self.M_ox_arr = np.array([])
        
    def initial_convergence(self):
        log = []

        # 入力パラメータの表示（古いバージョンのスタイル）
        print("-------------")
        print("input data")
        print("Requirement Thrust = ", self.F_req, "[N]")
        print("initial chamber pressure = ", self.Pc_def, "[MPa]")
        print("initial O/F = ", self.OF_def, "[-]")
        print("initial epsilon = ", self.epsilon_start, "[-]")
        print("initial Mdot propellant = ", self.mdot_new, "[kg/s]")
        print("Cstar efficient = ", self.eta_cstar, "[-]")
        print("Nozzle efficient = ", self.eta_nozzle, "[-]")

        # epsilon 調整
        (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1,
        self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1,
        self.Pe_tmp1, self.Mach_tmp1) = CEAInterface.compute(self.Pc_def, self.OF_def, epsilon=3)

        self.epsilon_new = \
            ((self.gamma_tmp1 + 1) / 2) ** (1 / (self.gamma_tmp1 - 1)) * \
            (self.Pa / self.Pc_def) ** (1 / self.gamma_tmp1) * \
            np.sqrt((self.gamma_tmp1 + 1) / (self.gamma_tmp1 - 1) *
                    (1 - (self.Pa / self.Pc_def) ** ((self.gamma_tmp1 - 1) / self.gamma_tmp1)))
        self.epsilon_new = 1 / self.epsilon_new
        print("calculated epsilon = ", self.epsilon_new, "[-]")

        # 初期CEA計算
        (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1,
        self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1,
        self.Pe_tmp1, self.Mach_tmp1) = CEAInterface.compute(self.Pc_def, self.OF_def, self.epsilon_new)

        self.Pe_old = self.Pe_tmp1
        self.diff_exit = 2
        self.i = 0
        self.j = 1

        self.R_tmp1 = self.R_univ / self.Mole_tmp1
        self.a_tmp1 = np.sqrt(self.gamma_tmp1 * self.R_tmp1 * self.T_e_tmp1)
        self.Ve_tmp1 = self.a_tmp1 * self.Mach_tmp1
        self.F = self.mdot_new * self.Ve_tmp1
        self.diff_F = self.F_req - self.F

        while abs(self.diff_F) > 0.1:
            self.mdot_new = self.mdot_old + 0.0001 if self.diff_F >= 0.1 else self.mdot_old - 0.0001
            self.mdot_old = self.mdot_new

            (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1,
            self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1,
            self.Pe_tmp1, self.Mach_tmp1) = CEAInterface.compute(self.Pc_def, self.OF_def, self.epsilon_new)

            # 出口速度計算
            self.R_tmp1 = self.R_univ / self.Mole_tmp1
            self.a_tmp1 = np.sqrt(self.gamma_tmp1 * self.R_tmp1 * self.T_e_tmp1)
            self.Ve_tmp1 = self.a_tmp1 * self.Mach_tmp1

            # スロート断面積計算
            self.At_new = 4 * self.eta_cstar * self.Cstar_tmp1 * self.mdot_new / (math.pi * self.Pc_def * 10 ** 6)

            # 出口マッハ数
            self.Me_new = np.sqrt((2 * ((self.Pc_def / self.Pe_init) ** ((self.gamma_tmp1 - 1) / self.gamma_tmp1) - 1)) / (self.gamma_tmp1 - 1)) + 0.01

            # 出口面積と膨張比
            self.Ae_new = self.At_new * (((1 + (((self.gamma_tmp1 - 1) * self.Me_new ** 2) / 2)) / ((self.gamma_tmp1 + 1) / 2)) ** ((self.gamma_tmp1 + 1) / (2 * (self.gamma_tmp1 - 1)))) / self.Me_new
            self.epsilon_new = self.Ae_new / self.At_new

            # 推力計算
            self.CF_tmp1 = self.CF_tmp1 + (self.Pe_tmp1 - self.Pa) * self.epsilon_new / self.Pc_def
            self.F = self.CF_tmp1 * self.Cstar_tmp1 * self.eta * self.mdot_new
            self.diff_F = self.F_req - self.F
            self.diff_exit = self.Pe_tmp1 - self.Pe_init

            # スロート径、出口径
            self.Dt = 2 * np.sqrt(self.At_new / math.pi)
            self.De = 2 * np.sqrt(self.Ae_new / math.pi)

            print("-------------")
            print("Thrust = ", self.F, "[N]")
            print("diff_F = ", self.diff_F, "[N]")
            print("mdot = ", self.mdot_new, "[kg/s]")
            print("Pe = ", self.Pe_tmp1, "[MPa]")
            print("diff_atm = ", abs(self.diff_exit))
            print("epsilon_new = ", self.epsilon_new)
            print("Dt, De = ", self.Dt, self.De, "[m]")

            if self.diff_exit == 0.0:
                self.i += 1
                if self.i > 3:
                    print("収束失敗")
                    break

        # Kstar計算（古いバージョンの定義）
        self.Kstar = (self.OF_def / (self.OF_def + 1)) * self.mdot_new / np.sqrt(2 * self.rho_ox_init * ((self.Ptank_init - self.Pc_init) * 1e6))

        print("-------------")
        print("Thrust = ", self.F, "[N]")
        print("diff_F = ", self.diff_F, "[N]")
        print("mdot = ", self.mdot_new, "[kg/s]")
        print("Pe = ", self.Pe_tmp1, "[MPa]")
        print("diff_atm = ", abs(self.diff_exit))
        print("epsilon_new = ", self.epsilon_new)
        print("Dt, De = ", self.Dt, self.De, "[m]")
        print("Kstar = ", self.Kstar)
        print("END")
        print("-------------")
    
    def integration_simulation(self):
        #------------------------#
        # 積分計算の開始（古い版に合わせた変数・配列運用）
        #------------------------#
        self.mdot_start = self.mdot_new
        self.mdot_ox_init = (self.OF_def / (self.OF_def + 1)) * self.mdot_start  # 初期酸化剤流量[kg/s]
        self.mdot_f_init = (1 / (self.OF_def + 1)) * self.mdot_start            # 初期燃料流量[kg/s]

        # K*（古い版の定義を使用）
        self.Kstar = self.mdot_ox_init / np.sqrt(2 * self.rho_ox_init * ((self.Ptank_init - self.Pc_init) * 1e6))
        print("mdot_ox", self.mdot_ox_init)
        print("Kstar = ", self.Kstar)

        # 燃焼面積と燃料長（古い版のロジック）
        self.Ap_req = self.mdot_f_init / (self.rho_f_start * self.a_ox * ((4 * self.mdot_ox_init) / (math.pi * self.Df_init ** 2)) ** self.n_ox)
        self.Lf = self.Ap_req / (self.Df_init * math.pi)
        print("Lf = ", self.Lf)
        print("Ap = ", self.Ap_req)
        print("mdot_f = ", self.mdot_f_init)
        print("O/F = ", self.mdot_ox_init / self.mdot_f_init)

        print("---------------START INTEGRATION---------------")

        # 質量・時間・初期状態
        self.Mass_ox = self.Vol_ox * self.rho_ox_init     # [kg]（古い版は Vol_ox[m^3] * rho[kg/m^3]）
        self.Mass_ox_remain = self.Mass_ox
        print("Mass_ox = ", self.Mass_ox)
        self.delta_t = 0.001  # [s]

        self.Ptank_tmp1 = self.Ptank_init
        self.Pc_tmp1 = self.Pc_def
        self.Df = self.Df_init

        # 経過時間など
        self.t = 9000
        self.k = 0

        print("Ap_req = ", self.Ap_req)
        print(self.Ptank_tmp1)
        print(self.Pc_tmp1)
        print("------------")

        # 初期CEA（古い版の epsilon_new を使用）
        (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1,
        self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1,
        self.Pe_tmp1, self.Mach_tmp1) = CEAInterface.compute(self.Pc_tmp1, self.OF_def, self.epsilon_new)

        # ログ配列初期化（古い版の配列名に合わせる）
        self.OX = self.Mass_ox * 1000  # 古い版の単位運用に合わせて g 表記
        self.Pt_arr = np.array([self.Ptank_tmp1])
        self.Pc_int_arr = np.array([self.Pc_tmp1])
        self.F_arr = np.array([self.F])
        self.OF_arr = np.array([self.OF_def])
        self.Ap_arr = np.array([self.Ap_req])
        self.mdot_arr = np.array([self.mdot_ox_init + self.mdot_f_init])
        self.Cstar_arr = np.array([self.Cstar_tmp1])
        self.CF_arr = np.array([self.CF_tmp1])
        self.F_fte_arr = np.array([self.F])
        self.M_ox_arr = np.array([self.Mass_ox_remain])
        self.mdot_ox_arr = np.array([self.mdot_ox_init])
        self.gamma_arr = np.array([self.gamma_tmp1])

        self.Ap = self.Ap_req
        self.It = 0

        # 外気圧（古い版では Pa_max を使用）
        self.Pa_tmp1 = self.Pa_max
        print("epsilon_new = ", self.epsilon_new)

        # 出口面積（Dt は初期収束で確定、epsilon_new 固定）
        self.At_new = math.pi * (self.Dt ** 2) / 4
        self.Ae_new = self.At_new * self.epsilon_new

        #====================#
        # 積分ループ
        #====================#
        while self.Mass_ox_remain >= 1:
            # Δp と流量（古い版の K* 流量式）
            self.delta_p = (self.Ptank_tmp1 - self.Pc_tmp1) * 1e6
            self.mdot_ox = self.Kstar * np.sqrt(2 * self.rho_ox_init * self.delta_p)  # g/ms 相当のスケール記載があるが、以降の質量更新と整合
            self.mdot_f = self.Ap * self.rho_f_start * self.a_ox * ((4 * self.mdot_ox) / (math.pi * self.Df ** 2)) ** self.n_ox

            print("Df = ", self.Df)

            # rdot と燃料径の更新（古い版の式）
            self.rdot = self.a_ox * ((4 * self.mdot_ox) / (math.pi * self.Df ** 2)) ** self.n_ox
            print(self.rdot)
            self.Df = self.Df + (2 * self.rdot / 1000)  # 両側燃焼で直径増加、mm→m 換算
            self.Ap = self.Df * math.pi * self.Lf

            print("mdot_ox = ", self.mdot_ox, "[g/ms]")
            print("mdot_f = ", self.mdot_f, "[g/ms]")

            # 瞬時 OF
            self.OF_tmp1 = self.mdot_ox / self.mdot_f
            print("OF_tmp1", self.OF_tmp1)

            # CEA（Pc, OF, epsilon 固定）
            (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1,
            self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1,
            self.Pe_tmp1, self.Mach_tmp1) = CEAInterface.compute(self.Pc_tmp1, self.OF_tmp1, self.epsilon_new)

            # 気体定数・音速
            self.R_tmp1 = self.R_univ / self.Mole_tmp1
            self.a_tmp1 = np.sqrt(self.gamma_tmp1 * self.R_tmp1 * self.T_e_tmp1)

            # 推力計算（古い版の2系統）
            self.F_fte = self.eta * ((self.mdot_ox + self.mdot_f) * self.a_tmp1 * self.Mach_tmp1) + (self.Pe_tmp1 - self.Pa_tmp1) * self.Ae_new
            self.CF_tmp1 = self.CF_tmp1 + (self.Pe_tmp1 - self.Pa) * self.epsilon_new / self.Pc_tmp1
            self.F_new = self.eta * self.Cstar_tmp1 * (self.mdot_ox + self.mdot_f) * self.CF_tmp1

            print("F = ", self.F_new)
            print("Pe = ", self.Pe_tmp1)

            # 残量更新（古い版の運用に合わせ、mdot_ox をそのまま差し引く）
            self.Mass_ox_remain = self.Mass_ox_remain - self.mdot_ox

            # タンク圧更新（線形スケーリング）
            self.Ptank_tmp1 = self.Ptank_fin + (self.Ptank_init - self.Ptank_fin) * (self.Mass_ox_remain / self.Mass_ox)

            # 燃焼室圧更新（C* ベース）
            self.Pc_tmp1 = 4 * self.Cstar_tmp1 * (self.mdot_ox + self.mdot_f) / (math.pi * self.Dt ** 2) / 1e6

            self.k = self.k + 1

            print("Pc_tmp1 = ", self.Pc_tmp1)
            print("Ptank_tmp1 = ", self.Ptank_tmp1)
            print("Pt = ", self.Ptank_tmp1)
            print("Mass_ox = ", self.Mass_ox)
            print("Remain ox = ", self.Mass_ox_remain)
            print("Lf = ", self.Lf)
            print("k = ", self.k)
            print("---------------")

            # ログ配列（古い版に準拠）
            self.Pt_arr = np.append(self.Pt_arr, self.Ptank_tmp1)
            self.Pc_int_arr = np.append(self.Pc_int_arr, self.Pc_tmp1)
            self.F_arr = np.append(self.F_arr, self.F_new)
            self.OF_arr = np.append(self.OF_arr, self.OF_tmp1)
            self.Ap_arr = np.append(self.Ap_arr, self.Ap)
            self.mdot = self.mdot_ox + self.mdot_f
            self.mdot_arr = np.append(self.mdot_arr, self.mdot)
            self.Cstar_arr = np.append(self.Cstar_arr, self.Cstar_tmp1)
            self.CF_arr = np.append(self.CF_arr, self.CF_tmp1)
            self.F_fte_arr = np.append(self.F_fte_arr, self.F_fte)
            self.M_ox_arr = np.append(self.M_ox_arr, self.Mass_ox_remain)
            self.mdot_ox_arr = np.append(self.mdot_ox_arr, self.mdot_ox)
            self.gamma_arr = np.append(self.gamma_arr, self.gamma_tmp1)

            # インパルス積分（古い版の delta_t）
            self.It = self.It + self.F_new * self.delta_t

    def plot_and_save_results(self):
        # ---------------- #
        #      result
        # ---------------- #
        print("----------RESULT----------")
        print("Kstar = ", self.Kstar)
        print("O/F_init = ", self.OF_def, "[-]")
        print("It = ", self.It, "[Ns]")
        print("Lf = ", self.Lf * 1000, "[mm]")
        print("Dt = ", self.Dt * 1000, "[mm]")
        print("Df_init = ", self.Df_init * 1000, "[mm]")
        print("Df_final = ", self.Df * 1000, "[mm]")
        print("F_ave =", self.It * 1000 / self.k, "[N]")
        print("F_init = ", self.F_req, "[N]")
        print(self.Ae_new)

        x = np.arange(len(self.F_arr))
        x = x/1000
        plt.plot(x, self.F_arr, label="CF")
        plt.plot(x, self.F_fte_arr, label="mdot*Pe")
        #plt.plot(x,OX_arr)
        #plt.plot(x,eps_arr)
        plt.show()

        plt.plot(x, self.Pt_arr, label="Pt")
        plt.plot(x, self.Pc_int_arr, label="Pc")
        plt.show()

        filename = f"FLXresult_forsimuration.csv"
        F_values = np.stack([self.F_arr, self.F_fte_arr, self.Pt_arr, self.Pc_int_arr, self.OF_arr, self.mdot_arr, self.Cstar_arr, self.CF_arr, self.M_ox_arr, self.mdot_ox_arr, self.gamma_arr]).T
        #F_values = np.stack([x,self.F_arr, self.F_fte_arr,self.mdot_ox_arr]).T
        with open(filename, 'w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file, quoting=csv.QUOTE_NONE)
            writer.writerows(F_values)
    
    # 初期状態のみやる子
    def run_initialstate(self):
        self.initial_convergence()

    # 全部やる子
    def run_simulation(self):
        self.initial_convergence()
        self.integration_simulation()
        self.plot_and_save_results()