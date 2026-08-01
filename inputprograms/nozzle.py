import numpy as np
import cantera as ct
#from numba import jit
from inputprograms import kyleniemeyer as k

def nozzle_flow(chamber_props, gas, P_chamber, P_exit, gas_origin, const, stoich_coeffs, mode=0, nfz = 2):
    """
    等エントロピー展開を仮定したCDノズルのスロート・出口状態を計算し、
    CanteraでTP再計算して物性を出力する

    Parameters
    ----------
    chamber_props : Dict
        燃焼室平衡後の値を格納している子
    gas : ct.Solution
        燃焼室平衡後のgasオブジェクト
    P_chamber : float
        燃焼室圧力 [Pa]
    P_exit : float
        出口圧力 [Pa]（mode=0のとき使用）
    gas_origin : ct.Solution
        元のgasオブジェクト（組成保持用）
    const : dict
        {"gamma": γ固定値, "Cstar": 特性速度}
    mode : int
        0 → 出口圧を大気圧に設定
        ≠0 → 開口比 Ae/At を指定値として出口計算
    nfz : int
        凍結流の判定
        1~3で決定し，1は燃焼室，2はスロート，3はノズル出口で凍結する
        初期設定は2
    """

    gamma = const["gamma"]
    R = ct.gas_constant / gas.mean_molecular_weight

    # スロート条件（M=1）
    T_throat = chamber_props["T"] * (2 / (gamma + 1))
    P_throat = P_chamber * (2 / (gamma + 1)) ** (gamma / (gamma - 1))

    # throat状態をCanteraで再計算
    gas_throat = gas_origin
    gas_throat.TP = T_throat, P_throat
    gas_throat.X = gas.X
    gas_throat.equilibrate('HP')
    gas_throat.TP = T_throat, P_throat
    gas_throat.equilibrate('SP')
    X_throat = gas_throat.X
    derivs = k.get_thermo_derivatives(gas_throat, stoich_coeffs)
    _, _, cp, gamma_s = k.get_thermo_properties(
        gas_throat, stoich_coeffs, derivs[0], derivs[1], derivs[2]
    )
    throat_props = {
        "T": T_throat,
        "P": P_throat,
        "rho": gas_throat.density,
        "H": gas_throat.enthalpy_mass/1000,
        "U": gas_throat.int_energy_mass/1000,
        "G": gas_throat.gibbs_mass/1000,
        "S": gas_throat.entropy_mass/1000,
        "M": gas_throat.mean_molecular_weight,
        "Cp": cp/1000,
        "Gamma": gamma_s,
        "a": gas_throat.sound_speed,
        "Mach": 1.0
    }
    throat_perf = nozzle_performance(throat_props["Gamma"], R, throat_props["T"],
                                P_chamber, throat_props["P"], ct.one_atm, throat_props["Mach"], const["Cstar"])

    # nfz = 2でこの先を凍結させる．
    if nfz == 2:
        # 各パラメータの計算
        # canteraに物性値を計算してもらう
        gas_exit_frozen = gas_origin
        gas_exit_frozen.X = X_throat
        entropy_throat = throat_props["S"] * 1000 # 出力用に1000で割ったのを戻す

        # 出口マッハ数の確定
        # mode != 0で開口比から計算
        if mode != 0:   # mode = epsilon（開口比指定）
            Ae_At_target = mode  # 指定開口比

            # 初期推定
            P_exit = ct.one_atm

            def compute_Ae_At_from_pressure(P_trial):
                # 仮の出口圧 P_trial を与えたときの Ae/At を計算する
                def exit_temperture_calc(T):
                    gas_exit_frozen.TP = T, P_trial
                    return gas_exit_frozen.entropy_mass - entropy_throat

                def df_T(T):
                    eps = 1e-2
                    return (exit_temperture_calc(T+eps) - exit_temperture_calc(T-eps)) / (2*eps)

                T_exit = 1800
                for _ in range(50):
                    T_exit -= exit_temperture_calc(T_exit) / df_T(T_exit)
                    if abs(exit_temperture_calc(T_exit)) < 1e-4:
                        break

                # gamma を計算
                gamma = gas_exit_frozen.cp / gas_exit_frozen.cv

                # --- 出口マッハ数を計算（圧力比から） ---
                P_ratio = P_trial / P_throat

                def f_M(M):
                    return ((gamma + 1) / (2 + (gamma - 1)*M**2))**(gamma/(gamma - 1)) - P_ratio

                def df_M(M):
                    eps = 1e-6
                    return (f_M(M+eps) - f_M(M-eps)) / (2*eps)

                M_exit = 2.0
                for _ in range(50):
                    M_exit -= f_M(M_exit) / df_M(M_exit)
                    if abs(f_M(M_exit)) < 1e-8:
                        break

                # --- Ae/At を計算 ---
                Ae_At_calc = (1/M_exit) * ((2/(gamma+1))*(1+(gamma-1)/2*M_exit**2))**((gamma+1)/(2*(gamma-1)))
                return Ae_At_calc, M_exit, T_exit, gamma

            # --- Newton-Raphson で P_exit を収束させる ---
            def F(P):
                Ae_at_calc, _, _, _ = compute_Ae_At_from_pressure(P)
                return Ae_at_calc - Ae_At_target

            def dF(P):
                eps = 0.01 * P
                return (F(P+eps) - F(P-eps)) / (2*eps)

            for _ in range(50):
                P_exit -= F(P_exit) / dF(P_exit)
                if abs(F(P_exit)) < 1e-6:
                    break

            # --- 収束した P_exit で最終状態を計算 ---
            _, M_exit, T_exit, gamma = compute_Ae_At_from_pressure(P_exit)

        # mode = 0で，出口大気圧
        elif mode == 0:
            # Newton-Raphsonでentropyを収束させる
            def exit_temperture_calc(T):
                gas_exit_frozen.TP = T, P_exit
                Entropy = gas_exit_frozen.entropy_mass
                return Entropy - entropy_throat
            def df(T):
                eps = 1e-2
                return (exit_temperture_calc(T+eps) - exit_temperture_calc(T-eps)) / (2*eps)
            T_exit = 1800 # 初期値
            for _ in range(50):
                T_exit -= exit_temperture_calc(T_exit)/df(T_exit)
                if abs(exit_temperture_calc(T_exit)) < 1e-4:
                    break

            # gammaを計算する
            gamma = gas_exit_frozen.cp / gas_exit_frozen.cv
        
            # 出口圧を大気圧に設定
            P_ratio = P_exit / P_throat
            # Newton-RaphsonでM_exitを解く
            def P_ratio_from_M(M):
                return ((gamma + 1) / (2 + (gamma - 1)*M**2))**((gamma / (gamma - 1))) - P_ratio
            def df(M):
                eps = 1e-6
                return (P_ratio_from_M(M+eps)-P_ratio_from_M(M-eps))/(2*eps)
            M_exit = 2.0  # 初期値
            for _ in range(50):
                M_exit -= P_ratio_from_M(M_exit)/df(M_exit)
                if abs(P_ratio_from_M(M_exit)) < 1e-8:
                    break

        # 等エントロピー流れ計算
        P_exit = throat_props["P"] * ((gamma + 1) / (2 + (gamma - 1) * M_exit ** 2)) ** (gamma / (gamma - 1))
        rho_exit = throat_props["rho"] * ((gamma + 1) / (2 + (gamma - 1) * M_exit ** 2)) ** (1 / (gamma - 1))
        S_exit = throat_props["S"] # 等エントロピー
        Mw_exit = throat_props["M"] # 反応しないので変化なし

        # 全パラメータ代入
        exit_props = {
        "T": T_exit,
        "P": P_exit,
        "rho": rho_exit,
        "H": gas_exit_frozen.enthalpy_mass / 1000,
        "U": gas_exit_frozen.int_energy_mass / 1000,
        "G": gas_exit_frozen.gibbs_mass / 1000,
        "S": S_exit,
        "M": Mw_exit,
        "Cp": gas_exit_frozen.cp/1000,
        "Gamma": gamma, # 凍結しているので組成変化なし
        "a": gas_exit_frozen.sound_speed,
        "Mach": M_exit
        }
        exit_perf = nozzle_performance(exit_props["Gamma"], R, exit_props["T"],
                                    P_chamber, exit_props["P"], ct.one_atm, exit_props["Mach"], const["Cstar"])
    # ここからnfz =3
    # 出口条件
    elif nfz == 3:
        if mode == 0:
            # 出口圧を大気圧に設定
            P_ratio = P_exit / P_chamber
            M_exit = np.sqrt((2 / (gamma - 1)) * (P_ratio ** (-(gamma - 1) / gamma) - 1))
        else:
            # 開口比 Ae/At から出口マッハ数を逆算
            Ae_At = mode
            # Newton-RaphsonでM_exitを解く
            def f(M):
                return (1/M) * ((2/(gamma+1))*(1+(gamma-1)/2*M**2))**((gamma+1)/(2*(gamma-1))) - Ae_At
            def df(M):
                eps = 1e-6
                return (f(M+eps)-f(M-eps))/(2*eps)
            M_exit = 2.0  # 初期値
            for _ in range(50):
                M_exit -= f(M_exit)/df(M_exit)
                if abs(f(M_exit)) < 1e-8:
                    break

        T_exit = chamber_props["T"] / (1 + (gamma - 1) / 2 * M_exit ** 2)
        P_exit_calc = P_chamber * (1 + (gamma - 1) / 2 * M_exit ** 2) ** (-gamma / (gamma - 1))

        # exit状態をCanteraで再計算
        gas_exit = gas_origin
        gas_exit.TP = T_exit, P_exit_calc
        gas_exit.X = gas.X
        gas_exit.equilibrate('HP')
        gas_exit.TP = T_exit, P_exit_calc
        gas_exit.equilibrate('SP')
        derivs = k.get_thermo_derivatives(gas_exit, stoich_coeffs)
        _, _, cp, gamma_s = k.get_thermo_properties(
            gas_exit, stoich_coeffs, derivs[0], derivs[1], derivs[2]
        )

        # calculate from EOS
        T_exit = chamber_props["T"] / (1 + ((gamma_s - 1) / 2 * M_exit ** 2))
        P_exit_calc =  gas_exit.density * ct.gas_constant * T_exit / gas_exit.mean_molecular_weight

        exit_props = {
            "T": T_exit,
            "P": P_exit_calc,
            "rho": gas_exit.density,
            "H": gas_exit.enthalpy_mass/1000,
            "U": gas_exit.int_energy_mass/1000,
            "G": gas_exit.gibbs_mass/1000,
            "S": gas_exit.entropy_mass/1000,
            "M": gas_exit.mean_molecular_weight,
            "Cp": cp/1000,
            "Gamma": gamma,
            "a": gas_exit.sound_speed,
            "Mach": M_exit
        }
        exit_perf = nozzle_performance(exit_props["Gamma"], R, exit_props["T"],
                                    P_chamber, exit_props["P"], ct.one_atm, exit_props["Mach"], const["Cstar"])

    return throat_props, exit_props, throat_perf, exit_perf

#@jit(nopython=True, cache=True)
def nozzle_performance(gamma, R, T, Pc, Pe, Pa, M, Cstar):
    # 膨張比 Ae/At
    Ae_At = (1/M) * ((2/(gamma+1))*(1+(gamma-1)/2*M**2))**((gamma+1)/(2*(gamma-1)))

    # 推力係数 Cf
    term1 = np.sqrt((2*gamma**2/(gamma-1)) * (2/(gamma+1))**((gamma+1)/(gamma-1)) * (1-(Pe/Pc)**((gamma-1)/gamma)))
    term2 = (Pe-Pa)/Pc * Ae_At
    Cf = term1 + term2

    # 真空比推力 Ivac
    Ivac = Cstar * (Cf + (Pa/Pc)*Ae_At)

    # 比推力 Isp
    g0 = 9.80665
    Isp = Cstar * Cf / g0

    return {"Ae/At": Ae_At, "Cstar": Cstar, "Cf": Cf, "Ivac": Ivac, "Isp": Isp}