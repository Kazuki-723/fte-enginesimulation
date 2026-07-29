import numpy as np
import cantera as ct
from inputprograms import kyleniemeyer as k

def nozzle_flow(gas, P_chamber, P_exit, gas_origin, const, mode=0):
    """
    等エントロピー展開を仮定したCDノズルのスロート・出口状態を計算し、
    CanteraでTP再計算して物性を出力する

    Parameters
    ----------
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
    """

    gamma = const["gamma"]
    R = ct.gas_constant / gas.mean_molecular_weight

    # スロート条件（M=1）
    T_throat = gas.T * (2 / (gamma + 1))
    P_throat = P_chamber * (2 / (gamma + 1)) ** (gamma / (gamma - 1))

    # 出口条件
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

    T_exit = gas.T / (1 + (gamma - 1) / 2 * M_exit ** 2)
    P_exit_calc = P_chamber * (1 + (gamma - 1) / 2 * M_exit ** 2) ** (-gamma / (gamma - 1))

    # throat状態をCanteraで再計算
    gas_throat = gas_origin
    gas_throat.TP = T_throat, P_throat
    gas_throat.X = gas.X
    gas_throat.equilibrate('HP')
    gas_throat.TP = T_throat, P_throat
    gas_throat.equilibrate('SP')
    derivs = k.get_thermo_derivatives(gas)
    dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = k.get_thermo_properties(
        gas, derivs[0], derivs[1], derivs[2]
    )
    throat_props = {
        "T": gas_throat.T,
        "P": gas_throat.P,
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

    # exit状態をCanteraで再計算
    gas_exit = gas_origin
    gas_exit.TP = T_exit, P_exit_calc
    gas_exit.X = gas.X
    gas_exit.equilibrate('HP')
    gas_exit.TP = T_exit, P_exit_calc
    gas_exit.equilibrate('SP')
    derivs = k.get_thermo_derivatives(gas)
    dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = k.get_thermo_properties(
        gas, derivs[0], derivs[1], derivs[2]
    )
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
        "Gamma": gamma_s,
        "a": gas_exit.sound_speed,
        "Mach": M_exit
    }
    exit_perf = nozzle_performance(exit_props["Gamma"], R, exit_props["T"],
                                P_chamber, exit_props["P"], ct.one_atm, exit_props["Mach"], const["Cstar"])

    return throat_props, exit_props, throat_perf, exit_perf

def nozzle_performance(gamma, R, T, Pc, Pe, Pa, M, Cstar):
    # 膨張比 Ae/At
    Ae_At = (1/M) * ((2/(gamma+1))*(1+(gamma-1)/2*M**2))**((gamma+1)/(2*(gamma-1)))

    # 特性速度 C*
    Cstar = Cstar

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