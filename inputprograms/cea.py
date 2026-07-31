import cantera as ct
import numpy as np
from inputprograms import kyleniemeyer as k
from inputprograms.nozzle import nozzle_flow

def CEA(O_F, T_init, P_init, epsilon, nfz = 2):
    # O/F指定（酸化剤/燃料の質量比）
    O_F = O_F

    # 質量比計算
    mass_PMMA = 1.0
    mass_N2O = mass_PMMA * O_F
    total_mass = mass_PMMA + mass_N2O

    Y_PMMA = mass_PMMA / total_mass
    Y_N2O = mass_N2O / total_mass

    # GRI-Mech + PMMA拡張の統合
    gri = ct.Solution('gri30.yaml')
    pmma = ct.Solution('inputprograms/pmma_extensionV3.yaml')

    gri_species_names = {sp.name for sp in gri.species()}
    pmma_species_filtered = [sp for sp in pmma.species() if sp.name not in gri_species_names]

    gas_origin = ct.Solution(thermo='ideal-gas',
                    species=gri.species() + pmma_species_filtered,
                    reactions=gri.reactions() + pmma.reactions())
    gas = gas_origin
    # 初期条件（燃焼室）
    T_init = T_init
    P_init = P_init * 10 * ct.one_atm # P_initはMPa想定

    gas.TP = T_init, P_init
    gas.Y = {'N2O': Y_N2O, 'MMA': Y_PMMA}

    # 平衡計算
    gas.equilibrate('HP')
    T_eq = gas.T
    gas.TP = T_eq, P_init
    gas.equilibrate('SP')

    derivs = k.get_thermo_derivatives(gas)

    dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = k.get_thermo_properties(
        gas, derivs[0], derivs[1], derivs[2]
        )

    Cstar = np.sqrt((ct.gas_constant * gas.T / (gas.mean_molecular_weight * gamma_s)) * ((gamma_s + 1) / 2) ** ((gamma_s + 1) / (gamma_s - 1)))
    rho_calc = P_init * gas.mean_molecular_weight / (ct.gas_constant * gas.T)
    speed_sound = np.sqrt(ct.gas_constant * gas.T * gamma_s / gas.mean_molecular_weight)

    chamber_props = {
        "T": gas.T,
        "P": P_init,
        "rho": gas.density,
        "H": gas.enthalpy_mass/1000,
        "U": gas.int_energy_mass/1000,
        "G": gas.gibbs_mass/1000,
        "S": gas.entropy_mass/1000,
        "M": gas.mean_molecular_weight,
        "Cp": cp/1000,
        "Gamma": gamma_s,
        "a": gas.sound_speed,
        "Mach": 0
    }

    P_chamber = P_init
    P_exit = ct.one_atm  # 大気圧

    # 定数固定する値の指定
    const = {"Cstar": Cstar, "gamma": gamma_s}

    # 反応流凍結制御
    """
    1：燃焼室凍結(未実装)
    2：スロート凍結
    3：全反応流
    """
    nfz_nozzle = nfz

    throat_props, exit_props, throat_perf, exit_perf = nozzle_flow(chamber_props, gas, P_chamber, P_exit, gas_origin, const, mode = epsilon, nfz = nfz_nozzle)
    return chamber_props, throat_props, exit_props, throat_perf, exit_perf