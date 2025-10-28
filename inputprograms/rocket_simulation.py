import numpy as np
import math
import csv
import matplotlib.pyplot as plt
from inputprograms.rocket_constants import *
from inputprograms.cea_interface import CEAInterface
from inputprograms.iteration_logger import IterationLogger

# 定数定義
R_univ = 8314 # 一般気体定数 [J/mol-K]
Pa = 0.1013 # 大気圧 [MPa]

class RocketSimulation:
    def __init__(self):
        # 定数・初期パラメータのセットアップ
        self.R_univ = R_univ
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
        self.Ve_th = 0

        # 積分計算用の配列初期化
        self.Pt_arr = np.array([])
        self.Pc_int_arr = np.array([])
        self.F_arr = np.array([])
        self.OF_arr = np.array([])
        self.Ap_arr = np.array([])
        self.mdot_arr = np.array([])
        self.Cstar_arr = np.array([])
        self.CF_arr = np.array([])
        self.F_fte_arr = np.array([])
        self.M_ox_arr = np.array([])
        self.kstar_cd_list = np.array([])


    def initial_convergence(self, F_req, Pc_def, OF_def, mdot_new, Df_init, eta_cstar, eta_nozzle):
        log = []

        # 入力パラメータの設定
        self.F_req = F_req * eta_cstar * eta_nozzle
        self.Pc_def = Pc_def
        self.OF_def = OF_def
        self.mdot_new = mdot_new
        self.mdot_old = mdot_new
        self.Df_init = Df_init
        self.eta_cstar = eta_cstar
        self.eta_nozzle = eta_nozzle
        self.eta = eta_cstar * eta_nozzle

        #epsilon 調整
        (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1,
         self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1,
         self.Pe_tmp1, self.Mach_tmp1) = CEAInterface.compute(self.Pc_def, self.OF_def, epsilon=3)
        
        # epsilonの計算式, 怪しいので調べる
        print(self.gamma_tmp1)
        print(self.Pc_def)
        print(self.Pa)
        self.epsilon_new = \
        ((self.gamma_tmp1 + 1) / 2) ** (1 / (self.gamma_tmp1 - 1)) * \
        (self.Pa / self.Pc_def) ** (1 / self.gamma_tmp1) * \
        np.sqrt((self.gamma_tmp1 + 1) / (self.gamma_tmp1 - 1) * (1 - (self.Pa/ self.Pc_def) ** ((self.gamma_tmp1 - 1) / self.gamma_tmp1)))
        self.epsilon_new = 1/self.epsilon_new
        print("calcrated epsilon = ", self.epsilon_new, "[-]")

        # 初期CEA計算
        (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1,
         self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1,
         self.Pe_tmp1, self.Mach_tmp1) = CEAInterface.compute(self.Pc_def, self.OF_def, self.epsilon_new)

        self.Pe_old = self.Pe_tmp1
        self.diff_exit = 2
        self.iter_logger = IterationLogger()
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
            self.Me_new = np.sqrt(2 * self.R_tmp1 * self.T_c_tmp1 * (self.gamma_tmp1/(self.gamma_tmp1-1)) * (1- (self.Pe_tmp1/self.Pc_def)**((self.gamma_tmp1-1)/self.gamma_tmp1))) / \
            np.sqrt(self.gamma_tmp1 * self.R_tmp1 * self.T_e_tmp1)

            # 開口比、出口面積
            # self.epsilon_new = ((1 + (self.gamma_tmp1 - 1) * (self.Me_new**2) / 2)/(1 + (self.gamma_tmp1 - 1)/2)) ** ((self.gamma_tmp1+1) / 2*(self.gamma_tmp1-1)) / self.Me_new
            # self.epsilon_new = 1/self.epsilon_new
            self.Ae_new = (self.Pthroat_tmp1/self.Pa) ** (1/self.gamma_tmp1) * (1/self.Me_new) * self.At_new
            self.epsilon_new = self.Ae_new / self.At_new
            print("calc.epsilon  = ", self.epsilon_new, "[-]")

            # 出口面積
            self.Ae_new = self.At_new * self.epsilon_new
            print("exit pres. = ",self.Pe_tmp1)

            #推力計算
            self.CF_tmp1 = self.CF_tmp1 + (self.Pe_tmp1 - self.Pa) * self.epsilon_new / self.Pc_def
            self.F = self.CF_tmp1 * self.Cstar_tmp1 * self.eta * self.mdot_new
            self.diff_F = self.F_req - self.F

            self.Dt = 2 * np.sqrt(self.At_new / math.pi)
            self.De = 2 * np.sqrt(self.Ae_new / math.pi)

            log.append(f"--- Iteration {self.j} ---")
            log.append(f"Thrust = {self.F:.3f} [N]")
            log.append(f"diff_F = {self.diff_F:.6f} [N]")
            log.append(f"mdot = {self.mdot_new:.6f} [kg/s]")
            log.append(f"Pe = {self.Pe_tmp1:.4f} [MPa]")
            log.append(f"epsilon_new = {self.epsilon_new:.4f}")
            log.append(f"Dt = {self.Dt:.4f} m, De = {self.De:.4f} m")
            self.iter_logger.append(self.j, self.F, self.mdot_new, self.Pe_tmp1, self.epsilon_new)
            self.j += 1

        # 初期状態の計算
        self.mdot_ox_init = (self.OF_def / (self.OF_def + 1)) * self.mdot_new  # 初期酸化剤流量[kg/s]
        self.mdot_f_init = (1 / (self.OF_def + 1)) * self.mdot_new  # 初期燃料流量[kg/s]

        # Discharge coef. * orifice cross section
        self.Kstar = self.mdot_ox_init / np.sqrt(2 * self.rho_ox_init * ((self.Ptank_init - self.Pc_init) * 1e6))

        # O/F, 燃料形状
        self.OF_tmp1 = self.mdot_ox_init / self.mdot_f_init
        self.Ap_req = self.mdot_f_init / (self.rho_f_start * self.a_ox * ((4 * self.mdot_ox_init) / (math.pi * self.Df_init ** 2)) ** self.n_ox)  # 定義したOFを実現するのに必要な燃焼面積
        self.Lf = self.Ap_req / (self.Df_init * math.pi)

        log.append("-------------")
        log.append(f"最終推力 = {self.F:.3f} [N]")
        log.append(f"最終mdot = {self.mdot_new:.6f} [kg/s]")
        log.append(f"最終Pe = {self.Pe_tmp1:.4f} [MPa]")
        log.append(f"最終epsilon = {self.epsilon_new:.4f}")
        log.append(f"Dt = {self.Dt:.4f} m")
        log.append(f"De = {self.De:.4f} m")
        log.append(f"K* = {self.Kstar}")
        log.append(f"初期酸化剤流量 = {self.mdot_ox_init:.6f}")
        log.append(f"初期燃料流量 = {self.mdot_f_init:.6f}")
        log.append(f"初期燃料内径(入力値) = {self.Df_init:.6f}")
        log.append(f"燃料長さ = {self.Lf:.6f}")

        log.append("-------------")

        print("-------------")
        print("Thrust = ", self.F, "[N]")
        print("mdot = ", self.mdot_new, "[kg/s]")
        print("Pe = ", self.Pe_tmp1, "[MPa]")
        print("epsilon_new = ", self.epsilon_new)
        print("Dt, De = ", self.Dt, self.De, "[m]")
        print("Kstar = ", self.Kstar)
        print("mdot_ox = ", self.mdot_ox_init, "[kg/s]")
        print("mdot_f = ", self.mdot_f_init, "[kg/s]")
        print("Df(input value) = ", self.Df_init, "[kg/s]")
        print("Lf = ", self.Lf, "[m]")
        print("END")
        print("-------------")

        # Kstarグラフ描画
        self.cd_values = np.linspace(0.1, 0.9, 50)
        for i in range(len(self.cd_values)):
            self.kstar_cd_list = np.append(self.kstar_cd_list, np.sqrt(4*self.Kstar / (self.cd_values[i] * math.pi)))

        print("END simulation")
        return "\n".join(log), self.kstar_cd_list, self.cd_values

    # resultのグラフ呼び出し関数
    def get_iteration_plot_base64(self, Dovalue, cdvalue):
        return self.iter_logger.get_base64_plot(Dovalue, cdvalue)

    def integration_simulation(self, Pc, Df, OF, eta_cstar, eta_nozzle, Kstar, epsilon,
                    Lf, mdot, V_tank, P_init, P_final, rho_ox, rho_fuel, a, n, F, Dt):
        # 入力の設定
        self.Pc_tmp1 = Pc
        self.Df = Df
        self.OF_tmp1 = OF
        self.eta_cstar = eta_cstar
        self.eta_nozzle = eta_nozzle
        self.eta = self.eta_cstar * eta_nozzle
        self.Kstar = Kstar
        self.epsilon_new = epsilon
        self.Lf = Lf
        self.mdot_start = mdot
        self.Vol_ox = V_tank
        self.Ptank_tmp1 = P_init
        self.Ptank_init = P_init
        self.Ptank_fin = P_final
        self.rho_ox_init = rho_ox
        self.rho_f_start = rho_fuel
        self.a_ox = a
        self.n_ox = n
        self.It = 0
        self.F = F
        self.Dt = Dt
        self.Ae_new = math.pi / 4 * self.Dt ** 2 * self.epsilon_new

        self.Pa_tmp1 = Pa
        self.mdot_ox_init = (self.OF_tmp1 / (self.OF_tmp1 + 1)) * self.mdot_start  # 初期酸化剤流量[kg/s]
        self.mdot_f_init = (1 / (self.OF_tmp1 + 1)) * self.mdot_start  # 初期燃料流量[kg/s]
        self.delta_t = 0.001  # 微小時間[s]
        self.k = 0
        self.iter_logger = IterationLogger()

        self.Ap_req  = self.Lf * self.Df * math.pi
        self.Ap = self.Ap_req
        print("---------------START INTEGRATION---------------")
        #====================#
        # 積分計算
        #====================#
        self.Mass_ox = self.Vol_ox * self.rho_ox_init * 1000  
        self.Mass_ox_remain = self.Mass_ox  

        # 初期状態CEAを回しなおす
        (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1, 
             self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1, 
             self.Pe_tmp1, self.Mach_tmp1) = CEAInterface.compute(self.Pc_tmp1, self.OF_tmp1, self.epsilon_new)

        # log配列
        self.OX = self.Mass_ox
        self.Pt_arr = np.array([self.Ptank_tmp1])
        self.Pc_int_arr = np.array([self.Pc_tmp1])
        self.F_arr = np.array([self.F])
        self.OF_arr = np.array([self.OF_tmp1])
        self.Ap_arr = np.array([self.Ap_req])
        self.mdot_arr = np.array([self.mdot_ox_init + self.mdot_f_init])
        self.Cstar_arr = np.array([self.Cstar_tmp1])
        self.CF_arr = np.array([self.CF_tmp1])
        self.F_fte_arr = np.array([self.F])
        self.M_ox_arr = np.array([self.Mass_ox_remain])
        self.mdot_ox_arr = np.array([self.mdot_ox_init])
        self.gamma_arr = np.array([self.gamma_tmp1])

        print("epsilon_new = ", self.epsilon_new)

        while self.Mass_ox_remain >= 1:  # 酸化剤残量が0になるまで計算を続ける
            self.delta_p = (self.Ptank_tmp1 - self.Pc_tmp1) * 1000000
            self.mdot_ox = (self.Kstar * np.sqrt(2 * self.rho_ox_init * self.delta_p))  # 微小時間における流量[g/ms]
            self.mdot_f = (self.Ap * self.rho_f_start * self.a_ox * ((4 * self.mdot_ox) / (math.pi * self.Df ** 2)) ** self.n_ox)  # 微小時間における燃料流量[g/ms]
            print("Df = ", self.Df)

            # rdotによる評価
            self.rdot = self.a_ox * ((4 * self.mdot_ox) / (math.pi * self.Df ** 2)) ** self.n_ox
            print(self.rdot)
            self.Df = self.Df + (2 * self.rdot / 1000)
            self.Ap = self.Df * math.pi * self.Lf

            print("mdot_ox = ", self.mdot_ox, "[g/ms]")
            print("mdot_f = ", self.mdot_f, "[g/ms]")

            self.OF_tmp1 = self.mdot_ox / self.mdot_f
            print("OF_tmp1", self.OF_tmp1)

            (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1, 
             self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1, 
             self.Pe_tmp1, self.Mach_tmp1) = CEAInterface.compute(self.Pc_tmp1, self.OF_tmp1, self.epsilon_new)

            self.R_tmp1 = self.R_univ / self.Mole_tmp1  # 気体定数
            self.a_tmp1 = np.sqrt(self.gamma_tmp1 * self.R_tmp1 * self.T_e_tmp1)  # 音速

            # 推力の計算
            self.F_fte = self.eta * ((self.mdot_ox + self.mdot_f) * self.a_tmp1 * self.Mach_tmp1) + (self.Pe_tmp1 - self.Pa_tmp1) * self.Ae_new
            # CFを実装する
            self.CF_tmp1 = self.CF_tmp1 + (self.Pe_tmp1 - self.Pa) * self.epsilon_new / self.Pc_tmp1
            self.F_new = self.eta * self.Cstar_tmp1 * (self.mdot_ox + self.mdot_f) * self.CF_tmp1

            # Pcの上書き
            print("F = ", self.F_new)
            print("Pe = ", self.Pe_tmp1)
            self.Mass_ox_remain = self.Mass_ox_remain - self.mdot_ox

            # Ptの計算
            self.Ptank_tmp1 = self.Ptank_fin + (self.Ptank_init - self.Ptank_fin) * (self.Mass_ox_remain / self.Mass_ox)
            self.Pc_tmp1 = 4 * self.Cstar_tmp1 * (self.mdot_ox + self.mdot_f) /(math.pi * self.Dt ** 2 ) / 1000000
            self.k = self.k + 1

            print("Pc_tmp1 = ", self.Pc_tmp1)
            print("Ptank_tmp1 = ", self.Ptank_tmp1)
            print("Pt = ", self.Ptank_tmp1)
            print("Mass_ox = ", self.Mass_ox)
            print("Remain ox = ", self.Mass_ox_remain)
            print("Lf = ", self.Lf)
            print("k = ", self.k)
            print("---------------")

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

            self.It = self.It + self.F_new * 0.001
        
        # print result
        print("----------RESULT----------")
        print("Kstar = ", self.Kstar)
        print("O/F_init = ", OF, "[-]")
        print("It = ", self.It, "[Ns]")
        print("Lf = ", Lf * 1000, "[mm]")
        print("Df_init = ", Df * 1000, "[mm]")
        print("Df_final = ", self.Df * 1000, "[mm]")
        print("F_ave =", self.It * 1000 / self.k, "[N]")
        time_ms = list(range(len(self.F_arr)))
        return time_ms, self.F_arr, self.F_fte_arr, self.OF_arr, self.Cstar_arr
    
    def get_evolution_plot_base64(self, time_ms, F_arr, F_fte_arr, OF_arr, Cstar_arr):
        return IterationLogger.plot_time_series(time_ms, F_arr, F_fte_arr, OF_arr, Cstar_arr)

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