import numpy as np
import math
from inputprograms.cea_interface import CEAInterface
from inputprograms.iteration_logger import IterationLogger
from inputprograms.interp_density import OxidizerDatabase

# ct-cea loading
from inputprograms import cea

# 定数定義
R_univ = 8314 # 一般気体定数 [mJ/mol-K]
Pa = 0.1013   # 大気圧 [MPa]

class RocketSimulation:
    def __init__(self):
        # 定数・初期パラメータのセットアップ
        self.R_univ = R_univ
        self.Pa = Pa
        self.ox_db = OxidizerDatabase()
        self.iter_logger = IterationLogger()
        self.diffuse_deg = 12 #ディフューザーの拡大角[deg]
        self.T_init = 290 # cantera計算時の初期温度[K]

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

    # Ptからrho_oxを計算
    def calc_rho_ox(self, pressure, phasedata):
        phase, rho = self.ox_db.get_density(pressure, phase = phasedata)
        return phase, rho
    
    # cantera計算からCEA結果調整
    def ctcea_compute(self, OF, T_init, P_init, epsilon):
        # ノズル内の凍結流定義．基本は2
        nfz = 2
        chamber_props, throat_props, exit_props, throat_perf, exit_perf = cea.CEA(OF, T_init, P_init, epsilon, nfz)

        #大結果シュート大会
        self.gamma_tmp1 = chamber_props['Gamma']
        self.Cstar_tmp1 = exit_perf['Cstar']
        self.CF_tmp1 = exit_perf['Cf']
        self.T_c_tmp1 = chamber_props["T"]
        self.T_t_tmp1 = throat_props["T"]
        self.T_e_tmp1 = exit_props["T"]
        self.Mole_tmp1 = chamber_props["M"]
        self.Pthroat_tmp1 = throat_props["P"] / 1e6 # MPaに直す
        self.Pe_tmp1 = exit_props["P"] / 1e6 # MPaに直す
        self.Mach_tmp1 = exit_props["Mach"]
        self.a_tmp1 = exit_props["a"]

        return self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1,\
        self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1,\
        self.Pe_tmp1, self.Mach_tmp1, self.a_tmp1

    # 初期値計算本体
    def initial_convergence(self, F_req, Pc_def, OF_def, mdot_new, Df_init, eta_cstar, eta_nozzle, Ptank_init, rho_ox_init, rho_f_start, a_ox, n_ox):
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
        self.rho_ox_init = rho_ox_init 
        self.Ptank_init = Ptank_init
        self.rho_f_start = rho_f_start
        self.a_ox = a_ox
        self.n_ox = n_ox

        # 最適epsilon調整
        # (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1,
        #  self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1,
        #  self.Pe_tmp1, self.Mach_tmp1, self.a_tmp1) = CEAInterface.compute(self.Pc_def, self.OF_def, epsilon=3)
        (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1,
        self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1,
        self.Pe_tmp1, self.Mach_tmp1, self.a_tmp1) = RocketSimulation.ctcea_compute(self, self.OF_def, self.T_init, self.Pc_def, epsilon=0) 
        
        
        # CEA入力明示
        print(self.gamma_tmp1)
        print(self.Pc_def)
        print(self.Pa)

        # epsilonの計算式
        self.epsilon_new = \
        ((self.gamma_tmp1 + 1) / 2) ** (1 / (self.gamma_tmp1 - 1)) * \
        (self.Pa / self.Pc_def) ** (1 / self.gamma_tmp1) * \
        np.sqrt((self.gamma_tmp1 + 1) / (self.gamma_tmp1 - 1) * (1 - (self.Pa/ self.Pc_def) ** ((self.gamma_tmp1 - 1) / self.gamma_tmp1)))
        self.epsilon_new = 1/self.epsilon_new
        print("calculated epsilon = ", self.epsilon_new, "[-]")

        # 初期CEA計算
        # (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1,
        #  self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1,
        #  self.Pe_tmp1, self.Mach_tmp1, self.a_tmp1) = CEAInterface.compute(self.Pc_def, self.OF_def, self.epsilon_new)
        
        (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1,
        self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1,
        self.Pe_tmp1, self.Mach_tmp1, self.a_tmp1) = RocketSimulation.ctcea_compute(self, self.OF_def, self.T_init, self.Pc_def, self.epsilon_new) 

        # iteration設定
        self.Pe_old = self.Pe_tmp1
        self.diff_exit = 2
        self.i = 0
        self.j = 1

        # ループ初期条件計算
        self.R_tmp1 = self.R_univ / self.Mole_tmp1
        self.Ve_tmp1 = self.a_tmp1 * self.Mach_tmp1
        self.F = self.mdot_new * self.Ve_tmp1
        self.diff_F = self.F_req - self.F

        # 目標推力に対して収束させる
        while abs(self.diff_F) > 0.1:
            # mdotを微小量動かして調整
            self.mdot_new = self.mdot_old + 0.0001 if self.diff_F >= 0.1 else self.mdot_old - 0.0001
            self.mdot_old = self.mdot_new

            # CEAによる計算
            # (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1,
            #  self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1,
            #  self.Pe_tmp1, self.Mach_tmp1, self.a_tmp1) = CEAInterface.compute(self.Pc_def, self.OF_def, self.epsilon_new)
            
            (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1,
            self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1,
            self.Pe_tmp1, self.Mach_tmp1, self.a_tmp1) = RocketSimulation.ctcea_compute(self, self.OF_def, self.T_init, self.Pc_def, self.epsilon_new) 

            # 出口速度計算
            self.R_tmp1 = self.R_univ / self.Mole_tmp1
            self.Ve_tmp1 = self.a_tmp1 * self.Mach_tmp1

            # スロート断面積計算
            self.At_new = self.eta_cstar * self.Cstar_tmp1 * self.mdot_new / (self.Pc_def * 10 ** 6)
            
            # 出口マッハ数
            # self.Me_new = np.sqrt(2 * self.R_tmp1 * self.T_c_tmp1 * (self.gamma_tmp1/(self.gamma_tmp1-1)) * (1- (self.Pe_tmp1/self.Pc_def)**((self.gamma_tmp1-1)/self.gamma_tmp1))) / \
            # np.sqrt(self.gamma_tmp1 * self.R_tmp1 * self.T_e_tmp1)

            # 開口比、出口面積
            self.Ae_new = (self.Pthroat_tmp1/self.Pa) ** (1/self.gamma_tmp1) * (1/self.Mach_tmp1) * self.At_new
            self.epsilon_new = self.Ae_new / self.At_new

            #推力計算
            self.CF_tmp1 = self.CF_tmp1 + (self.Pe_tmp1 - self.Pa) * self.epsilon_new / self.Pc_def
            self.F = self.CF_tmp1 * self.Cstar_tmp1 * self.eta * self.mdot_new
            self.diff_F = self.F_req - self.F

            self.Dt = 2 * np.sqrt(self.At_new / math.pi)
            self.De = 2 * np.sqrt(self.Ae_new / math.pi)

            # iteration log管理
            log.append(f"--- Iteration {self.j} ---")
            log.append(f"Thrust = {self.F:.3f} [N]")
            log.append(f"diff_F = {self.diff_F:.6f} [N]")
            log.append(f"mdot = {self.mdot_new:.6f} [kg/s]")
            log.append(f"Pe = {self.Pe_tmp1:.4f} [MPa]")
            log.append(f"epsilon_new = {self.epsilon_new:.4f}")
            log.append(f"Dt = {self.Dt:.4f} m, De = {self.De:.4f} m")
            self.iter_logger.append(self.j, self.F, self.mdot_new, self.Pe_tmp1, self.epsilon_new)

            # terminal出力管理
            print(f"--- Iteration {self.j} ---")
            print(f"Thrust = {self.F:.3f} [N]")
            print(f"diff_F = {self.diff_F:.6f} [N]")
            print(f"mdot = {self.mdot_new:.6f} [kg/s]")
            print(f"Pe = {self.Pe_tmp1:.4f} [MPa]")
            print(f"Dt = {self.Dt:.4f} m, De = {self.De:.4f} m")
            print(f"epsilon_new = {self.epsilon_new:.4f}")
            self.j += 1

        # 収束した初期状態の計算
        self.mdot_ox_init = (self.OF_def / (self.OF_def + 1)) * self.mdot_new  # 初期酸化剤流量[kg/s]
        self.mdot_f_init = (1 / (self.OF_def + 1)) * self.mdot_new             # 初期燃料流量[kg/s]

        # Kstar = Discharge coef. * orifice cross section
        self.Kstar = self.mdot_ox_init / np.sqrt(2 * self.rho_ox_init * ((self.Ptank_init - self.Pc_def) * 1e6))

        # O/F, 燃料形状
        self.OF_tmp1 = self.mdot_ox_init / self.mdot_f_init

        # 定義したOFを実現するのに必要な燃焼面積
        # ここで算出されるのは有効長さ
        self.Ap_req = self.mdot_f_init / (self.rho_f_start * self.a_ox * ((4 * self.mdot_ox_init) / (math.pi * self.Df_init ** 2)) ** self.n_ox)
        
        # Dfから燃料長さを計算
        self.Lf = self.Ap_req / (self.Df_init * math.pi)

        # 上の有効長さから全長を計算
        self.Lf_total = self.Lf + self.Df_init / 2 / math.tan(math.radians(self.diffuse_deg))

        # 最終結果log保存
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

        # 最終結果terminal出力
        print("-------------")
        print("Thrust(input value) = ", self.F, "[N]")
        print("chamber pressure(input value) = ", self.Pc_def, "[MPa]")
        print("O/F(input value) = ", self.OF_def, "[-]")
        print("mdot = ", self.mdot_new, "[kg/s]")
        print("Df(input value) = ", self.Df_init, "[m]")
        print("eta_cstar(input value) = ", self.eta_cstar, "[-]")
        print("eta_nozzle(input value) = ", self.eta_nozzle, "[-]")
        print("Kstar = ", self.Kstar)
        print("epsilon_new = ", self.epsilon_new)
        print("Lf(effective length) = ", self.Lf, "[m]")
        print("Lf(total length) = ", self.Lf_total, "[m]")
        print("Dt, De = ", self.Dt, self.De, "[m]")
        print(f"\n")
        print("Pe = ", self.Pe_tmp1, "[MPa]")
        print("mdot_ox = ", self.mdot_ox_init, "[kg/s]")
        print("mdot_f = ", self.mdot_f_init, "[kg/s]")
        print("END")
        print("-------------")

        # Kstarグラフ描画
        self.cd_values = np.linspace(0.1, 0.9, 50)
        for i in range(len(self.cd_values)):
            self.kstar_cd_list = np.append(self.kstar_cd_list, np.sqrt(4*self.Kstar / (self.cd_values[i] * math.pi)))

        print("END initial condition simulation")

        return "\n".join(log), self.kstar_cd_list, self.cd_values

    # 時間発展計算
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

        # 必要値の定義
        self.Pa_tmp1 = Pa
        self.mdot_ox_init = (self.OF_tmp1 / (self.OF_tmp1 + 1)) * self.mdot_start  # 初期酸化剤流量[kg/s]
        self.mdot_f_init = (1 / (self.OF_tmp1 + 1)) * self.mdot_start              # 初期燃料流量[kg/s]
        
        # 時間管理
        # 微小時間[s] 現状kg/sとg/msが数値上同じ値になることを利用した残念仕様があるので，変えると壊れる
        self.delta_t = 0.001
        # iteration管理  
        self.k = 0

        self.Ap_req  = self.Lf * self.Df * math.pi
        self.Ap = self.Ap_req
        print("---------------START INTEGRATION---------------")
        #====================#
        # 積分計算
        #====================#
        # タンク内部定義
        self.Mass_ox = self.Vol_ox * self.rho_ox_init * 1000  
        self.Mass_ox_remain = self.Mass_ox

        # 有効長さを定義
        diseffect_length = self.Df / 2 / math.tan(math.radians(self.diffuse_deg))
        self.Lf = self.Lf - diseffect_length

        # 初期状態CEA
        # (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1, 
        #      self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1, 
        #      self.Pe_tmp1, self.Mach_tmp1, self.a_tmp1) = CEAInterface.compute(self.Pc_tmp1, self.OF_tmp1, self.epsilon_new)
        
        (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1,
        self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1,
        self.Pe_tmp1, self.Mach_tmp1, self.a_tmp1) = RocketSimulation.ctcea_compute(self, self.OF_tmp1, self.T_init, self.Pc_tmp1, self.epsilon_new) 

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
        self.Df_arr = np.array([self.Df])

        print("epsilon_new = ", self.epsilon_new)

        # 終了時酸化剤質量を定義
        _, rho_ox_end =  RocketSimulation.calc_rho_ox(self, self.Ptank_fin, "gas")
        Mass_ox_end = self.Vol_ox * rho_ox_end * 1000
        print(Mass_ox_end)
        print(self.Mass_ox_remain)

        # 酸化剤が規定量になるまで時間発展
        while self.Mass_ox_remain >= Mass_ox_end:
            # 流量計算
            self.delta_p = (self.Ptank_tmp1 - self.Pc_tmp1) * 1000000
            # 微小時間における酸化剤流量[g/ms]
            self.mdot_ox = (self.Kstar * np.sqrt(2 * self.rho_ox_init * self.delta_p))  
            # 微小時間における燃料流量[g/ms]
            self.mdot_f = (self.Ap * self.rho_f_start * self.a_ox * ((4 * self.mdot_ox) / (math.pi * self.Df ** 2)) ** self.n_ox) 
            print("Df = ", self.Df)

            # rdot計算
            self.rdot = self.a_ox * ((4 * self.mdot_ox) / (math.pi * self.Df ** 2)) ** self.n_ox
            print(self.rdot)
            #燃料後退，反応表面積計算
            self.Df = self.Df + (2 * self.rdot / 1000)
            self.Ap = self.Df * math.pi * self.Lf

            print("mdot_ox = ", self.mdot_ox, "[g/ms]")
            print("mdot_f = ", self.mdot_f, "[g/ms]")

            # OF算出
            self.OF_tmp1 = self.mdot_ox / self.mdot_f
            print("OF_tmp1", self.OF_tmp1)

            # CEA計算
            # (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1, 
            #  self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1, 
            #  self.Pe_tmp1, self.Mach_tmp1, self.a_tmp1) = CEAInterface.compute(self.Pc_tmp1, self.OF_tmp1, self.epsilon_new)
            
            (self.gamma_tmp1, self.Cstar_tmp1, self.CF_tmp1, self.T_c_tmp1,
            self.T_t_tmp1, self.T_e_tmp1, self.Mole_tmp1, self.Pthroat_tmp1,
            self.Pe_tmp1, self.Mach_tmp1, self.a_tmp1) = RocketSimulation.ctcea_compute(self, self.OF_tmp1, self.T_init, self.Pc_tmp1, self.epsilon_new) 

            # 気体物性値評価
            self.R_tmp1 = self.R_univ / self.Mole_tmp1  # 気体定数

            # 推力の計算
            self.F_fte = self.eta * ((self.mdot_ox + self.mdot_f) * self.a_tmp1 * self.Mach_tmp1 + (self.Pe_tmp1 - self.Pa_tmp1) * self.Ae_new * 1e6)
            
            # CFの圧力補正
            self.CF_tmp1 = self.CF_tmp1 + (self.Pe_tmp1 - self.Pa) * self.epsilon_new / self.Pc_tmp1
            self.F_new = self.eta * self.Cstar_tmp1 * (self.mdot_ox + self.mdot_f) * self.CF_tmp1

            # 酸化剤残量の更新
            print("F = ", self.F_new)
            print("Pe = ", self.Pe_tmp1)
            self.Mass_ox_remain = self.Mass_ox_remain - self.mdot_ox

            # 次iterationへ投げる圧力の計算
            self.Ptank_tmp1 = self.Ptank_init + (self.Ptank_fin - self.Ptank_init) / (Mass_ox_end - self.Mass_ox) * (self.Mass_ox_remain - self.Mass_ox)
            self.Pc_tmp1 = 4 * self.eta_cstar * self.Cstar_tmp1 * (self.mdot_ox + self.mdot_f) /(math.pi * self.Dt ** 2 ) / 1000000
            self.k = self.k + 1

            # iteration log terminal管理
            print("Pc_tmp1 = ", self.Pc_tmp1)
            print("Ptank_tmp1 = ", self.Ptank_tmp1)
            print("Pt = ", self.Ptank_tmp1)
            print("Mass_ox = ", self.Mass_ox)
            print("Remain ox = ", self.Mass_ox_remain)
            print("Lf = ", self.Lf)
            print("k = ", self.k)
            print("---------------")

            # 配列管理
            self.Pt_arr = np.append(self.Pt_arr, self.Ptank_tmp1)
            self.Pc_int_arr = np.append(self.Pc_int_arr, self.Pc_tmp1)
            self.F_arr = np.append(self.F_arr, self.F_new)
            self.OF_arr = np.append(self.OF_arr, self.OF_tmp1)
            self.Ap_arr = np.append(self.Ap_arr, self.Ap)
            self.Df_arr = np.append(self.Df_arr, self.Df)
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
        print("end time evolution simulation")
        time_ms = list(range(len(self.F_arr)))
        evolution_result = np.stack([self.F_arr, self.F_fte_arr, self.Pt_arr, self.Pc_int_arr, self.OF_arr, self.mdot_arr, self.Df_arr, self.Cstar_arr, self.CF_arr, self.M_ox_arr, self.mdot_ox_arr, self.gamma_arr]).T
        return time_ms, self.F_arr, self.F_fte_arr, self.OF_arr, self.Cstar_arr, self.Pc_int_arr, self.Pt_arr, evolution_result, self.It
    
    # GUIでグラフを書くためだけに存在する関数たち

    # 初期条件のresultでの収束例歴表示関数
    def get_iteration_plot_base64(self, Dovalue, cdvalue):
        return self.iter_logger.get_base64_plot(Dovalue, cdvalue)

    # 時間発展のresultでの結果表示関数
    def get_evolution_plot_base64(self, time_ms, F_arr, F_fte_arr, OF_arr, Cstar_arr, Pc_arr, Pt_arr):
        return IterationLogger.plot_time_series(time_ms, F_arr, F_fte_arr, OF_arr, Cstar_arr, Pc_arr, Pt_arr)