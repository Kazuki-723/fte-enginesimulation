import cantera as ct
import numpy as np
from input import kyleniemeyer as k
from input.nozzle import nozzle_flow

def CEA(O_F, T_init, P_init):
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
    pmma = ct.Solution('input/pmma_extensionV3.yaml')

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


    throat_props, exit_props, throat_perf, exit_perf = nozzle_flow(gas, P_chamber, P_exit, gas_origin, const, mode = 0)

    # print("---input---")
    # print("OF = ", O_F,"[-]")
    # print("Pc = ",P_init / (10 * ct.one_atm),"[MPa]")
    # print("initial Temp.",T_init,"[K]")

    # print("------------------------------------------------------------")
    # print("Chamber vs Throat vs Exit (CEA-style)")
    # print("------------------------------------------------------------")
    # print(f"{'Property':<20}{'Chamber':>15}{'Throat':>15}{'Exit':>15}")
    # print(f"{'P [bar]':<20}{chamber_props['P']/ct.one_atm:>15.4f}{results['throat']['P']/ct.one_atm:>15.4f}{results['exit']['P']/ct.one_atm:>15.4f}")
    # print(f"{'T [K]':<20}{chamber_props['T']:>15.4f}{results['throat']['T']:>15.4f}{results['exit']['T']:>15.4f}")
    # print(f"{'rho [kg/m^3]':<20}{chamber_props['rho']:>15.4f}{results['throat']['rho']:>15.4f}{results['exit']['rho']:>15.4f}")
    # print(f"{'H [kJ/kg]':<20}{chamber_props['H']:>15.4f}{results['throat']['H']:>15.4f}{results['exit']['H']:>15.4f}")
    # print(f"{'U [kJ/kg]':<20}{chamber_props['U']:>15.4f}{results['throat']['U']:>15.4f}{results['exit']['U']:>15.4f}")
    # print(f"{'G [kJ/kg]':<20}{chamber_props['G']:>15.4f}{results['throat']['G']:>15.4f}{results['exit']['G']:>15.4f}")
    # print(f"{'S [kJ/kg/K]':<20}{chamber_props['S']:>15.4f}{results['throat']['S']:>15.4f}{results['exit']['S']:>15.4f}")
    # print(f"{'M [1/n]':<20}{chamber_props['M']:>15.4f}{results['throat']['M']:>15.4f}{results['exit']['M']:>15.4f}")
    # print(f"{'Cp [kJ/kg/K]':<20}{chamber_props['Cp']:>15.4f}{results['throat']['Cp']:>15.4f}{results['exit']['Cp']:>15.4f}")
    # print(f"{'Gamma [-]':<20}{chamber_props['Gamma']:>15.4f}{results['throat']['Gamma']:>15.4f}{results['exit']['Gamma']:>15.4f}")
    # print(f"{'Sonic Vel. [m/s]':<20}{chamber_props['a']:>15.4f}{results['throat']['a']:>15.4f}{results['exit']['a']:>15.4f}")
    # print(f"{'Mach number':<20}{chamber_props['Mach']:>15.4f}{results['throat']['Mach']:>15.4f}{results['exit']['Mach']:>15.4f}")
    # print("------------------------------------------------------------")

    # print("------------------------------------------------------------")
    # print("Performance Parameters")
    # print("------------------------------------------------------------")
    # print(f"{'Property':<15}{'Throat':>15}{'Exit':>15}")
    # print(f"{'Ae/At':<15}{throat_perf['Ae/At']:>15.4f}{exit_perf['Ae/At']:>15.4f}")
    # print(f"{'Cstar [m/s]':<15}{throat_perf['Cstar']:>15.2f}{exit_perf['Cstar']:>15.2f}")
    # print(f"{'Cf':<15}{throat_perf['Cf']:>15.4f}{exit_perf['Cf']:>15.4f}")
    # print(f"{'Ivac [m/s]':<15}{throat_perf['Ivac']:>15.2f}{exit_perf['Ivac']:>15.2f}")
    # print(f"{'Isp [s]':<15}{throat_perf['Isp']:>15.2f}{exit_perf['Isp']:>15.2f}")
    # print("------------------------------------------------------------")

    # 組成表示（濃度が高い種のみ）
    # for sp in gas.species_names:
    #     X = gas.X[gas.species_index(sp)]
    #     if X > 5e-6:
    #         print(f"{sp}: {X:.6f}")
    
    return chamber_props, throat_props, exit_props, throat_perf, exit_perf