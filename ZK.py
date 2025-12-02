import numpy as np
import csv
from dataclasses import dataclass
from typing import List
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import root

# ----------------------------
# ユーティリティ
# ----------------------------
def clamp(x, xmin, xmax):
    return max(min(x, xmax), xmin)

# ----------------------------
# 物性テーブル補間クラス
# ----------------------------
M_N2O = 44.0128e-3  # kg/mol
R_univ = 8.314

@dataclass
class PhaseTable:
    T: np.ndarray
    P_MPa: np.ndarray
    rho: np.ndarray
    u_kJ_per_mol: np.ndarray
    h_kJ_per_mol: np.ndarray
    cp_J_per_molK: np.ndarray
    viscosity_uPa_s: np.ndarray
    k_W_mK: np.ndarray
    surf_tension_N_m: np.ndarray = None

    # スプライン補間器
    f_P_MPa: CubicSpline = None
    f_rho: CubicSpline = None
    f_u_kJ_per_mol: CubicSpline = None
    f_h_kJ_per_mol: CubicSpline = None
    f_cp_J_per_molK: CubicSpline = None
    f_visc_uPa_s: CubicSpline = None
    f_k_W_mK: CubicSpline = None
    f_sigma: CubicSpline = None

    def build(self, is_liquid: bool):
        self.f_P_MPa = CubicSpline(self.T, self.P_MPa, bc_type="natural")
        self.f_rho = CubicSpline(self.T, self.rho, bc_type="natural")
        self.f_u_kJ_per_mol = CubicSpline(self.T, self.u_kJ_per_mol, bc_type="natural")
        self.f_h_kJ_per_mol = CubicSpline(self.T, self.h_kJ_per_mol, bc_type="natural")
        self.f_cp_J_per_molK = CubicSpline(self.T, self.cp_J_per_molK, bc_type="natural")
        self.f_visc_uPa_s = CubicSpline(self.T, self.viscosity_uPa_s, bc_type="natural")
        self.f_k_W_mK = CubicSpline(self.T, self.k_W_mK, bc_type="natural")
        if is_liquid and self.surf_tension_N_m is not None:
            self.f_sigma = CubicSpline(self.T, self.surf_tension_N_m, bc_type="natural")

    def _clip(self, Tq):
        return np.clip(Tq, self.T[0], self.T[-1])

    def Psat_Pa(self, Tq):
        return float(self.f_P_MPa(self._clip(Tq))) * 1e6

    def rho_kg_m3(self, Tq):
        return float(self.f_rho(self._clip(Tq)))

    def u_J_kg(self, Tq):
        return float(self.f_u_kJ_per_mol(self._clip(Tq))) * 1e3 / M_N2O

    def h_J_kg(self, Tq):
        return float(self.f_h_kJ_per_mol(self._clip(Tq))) * 1e3 / M_N2O

    def cp_J_kgK(self, Tq):
        return float(self.f_cp_J_per_molK(self._clip(Tq))) / M_N2O

    def mu_Pa_s(self, Tq):
        return float(self.f_visc_uPa_s(self._clip(Tq))) * 1e-6

    def k_W_mK(self, Tq):
        return float(self.f_k_W_mK(self._clip(Tq)))

    def sigma_N_m(self, Tq):
        if self.f_sigma is None:
            return None
        return float(self.f_sigma(self._clip(Tq)))

def load_phase_csv(path: str, is_liquid: bool) -> PhaseTable:
    T, P_MPa, rho = [], [], []
    u_kJ_mol, h_kJ_mol, cp_J_molK = [], [], []
    visc_uPa_s, k_W_mK, sigma_N_m = [], [], []

    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            T.append(float(row["Temperature (K)"]))
            P_MPa.append(float(row["Pressure (MPa)"]))
            rho.append(float(row["Density (kg/m3)"]))
            u_kJ_mol.append(float(row["Internal Energy (kJ/mol)"]))
            h_kJ_mol.append(float(row["Enthalpy (kJ/mol)"]))
            cp_J_molK.append(float(row["Cp (J/mol*K)"]))
            visc_uPa_s.append(float(row["Viscosity (uPa*s)"]))
            k_W_mK.append(float(row["Therm. Cond. (W/m*K)"]))
            if is_liquid and "Surf. Tension (N/m)" in row:
                sigma_N_m.append(float(row["Surf. Tension (N/m)"]))

    table = PhaseTable(
        T=np.array(T), P_MPa=np.array(P_MPa), rho=np.array(rho),
        u_kJ_per_mol=np.array(u_kJ_mol), h_kJ_per_mol=np.array(h_kJ_mol),
        cp_J_per_molK=np.array(cp_J_molK), viscosity_uPa_s=np.array(visc_uPa_s),
        k_W_mK=np.array(k_W_mK),
        surf_tension_N_m=(np.array(sigma_N_m) if sigma_N_m else None)
    )
    table.build(is_liquid=is_liquid)
    return table
# ----------------------------
# Peng–Robinson EOS utilities for N2O
# ----------------------------
R_univ = 8.314462618  # J/mol/K

# N2O constants
Tc_N2O = 309.57       # K
Pc_N2O = 7.245e6      # Pa
omega_N2O = 0.274     # acentric factor

def PR_params_N2O(T: float):
    # PR parameters (a, b, alpha) for N2O at temperature T
    kappa = 0.37464 + 1.54226 * omega_N2O - 0.26992 * omega_N2O**2
    Tr = T / Tc_N2O
    alpha = (1 + kappa * (1 - np.sqrt(Tr)))**2
    a = 0.45724 * (R_univ**2 * Tc_N2O**2 / Pc_N2O) * alpha
    b = 0.07780 * (R_univ * Tc_N2O / Pc_N2O)
    return a, b, alpha

def PR_pressure_from_molar_volume(T: float, v_m: float) -> float:
    """
    P(T, v_m) by Peng–Robinson EoS; v_m: molar volume [m3/mol]
    """
    a, b, _ = PR_params_N2O(T)
    term1 = R_univ * T / (v_m - b)
    term2 = a / (v_m * (v_m + b) + b * (v_m - b))
    return term1 - term2

def PR_Z_roots(T: float, P: float):
    """
    Solve PR cubic for compressibility factor Z at given T, P.
    Returns roots (real). Gas root is the largest Z.
    """
    a, b, _ = PR_params_N2O(T)
    A = a * P / (R_univ**2 * T**2)
    B = b * P / (R_univ * T)
    # Coefficients for Z^3 + c2 Z^2 + c1 Z + c0 = 0
    c2 = -(1 - B)
    c1 = A - 3*B**2 - 2*B
    c0 = -(A*B - B**2 - B**3)
    # Solve cubic
    coeffs = [1.0, c2, c1, c0]
    roots = np.roots(coeffs)
    real_roots = np.real(roots[np.isreal(roots)])
    return np.sort(real_roots)

def PR_rho_gas_from_PT(T: float, P: float) -> float:
    """
    Gas density from PR at given T,P using gas root Z.
    rho = (P * M) / (Z * R * T)
    """
    Zs = PR_Z_roots(T, P)
    Z = Zs[-1] if len(Zs) else 1.0
    return (P * M_N2O) / (Z * R_univ * T)

def PR_pressure_from_state(m_gas: float, T_gas: float, V_gas: float) -> float:
    """
    Compute gas pressure from state using PR via molar volume.
    """
    n_mol = m_gas / M_N2O
    v_m = V_gas / max(n_mol, 1e-12)
    return PR_pressure_from_molar_volume(T_gas, v_m)


# ----------------------------
# タンクモデル
# ----------------------------
@dataclass
class TankConfig:
    V_tank: float
    A_orifice: float
    Cd: float
    T_wall: float
    P_down: float
    R_g: float
    gamma_g: float
    P_down_gas: float

@dataclass
class TankState:
    m_liq: float
    m_gas: float
    T_liq: float
    T_gas: float

def ullage_volume(cfg: TankConfig, st: TankState, liquid: PhaseTable):
    V_liq = st.m_liq / liquid.rho_kg_m3(st.T_liq)
    return cfg.V_tank - V_liq

def tank_pressure(st: TankState, liquid: PhaseTable):
    return liquid.Psat_Pa(st.T_liq)

# ----------------------------
# Orifice mass flow models
# ----------------------------
def mass_flow_orifice_liquid(P_tank: float, T_liq: float, cfg: TankConfig, liquid: PhaseTable, k_flash: float = 1.0) -> float:
    """
    Liquid orifice flow with simple flashing clamp:
    """
    rho_l = liquid.rho_kg_m3(T_liq)
    dP = max(P_tank - cfg.P_down, 0.0)
    return cfg.Cd * cfg.A_orifice * np.sqrt(2.0 * rho_l * dP)

def mass_flow_orifice_gas_PR(P_tank: float, T_gas: float, cfg: TankConfig) -> float:
    """
    Gas orifice flow using upstream gas density from PR(T, P_tank).
    Non-choked Bernoulli-like formula for robustness.
    """
    rho_g = PR_rho_gas_from_PT(T_gas, P_tank)
    dP = max(P_tank - cfg.P_down_gas, 0.0)
    return cfg.Cd * cfg.A_orifice * np.sqrt(2.0 * rho_g * dP)

# 内部エネルギーから温度を逆算する関数
def invert_u_to_T(table: PhaseTable, u_target: float, T_guess: float = 290.0):
    f = lambda T: table.u_J_kg(T) - u_target
    T_sol = fsolve(f, T_guess)
    return float(T_sol[0])

def step(cfg: TankConfig, st: TankState, dt: float, liquid: PhaseTable, vapor: PhaseTable,
         alpha_evap: float = 0.3, alpha_sat: float = 0.5, clip_frac: float = 0.01,
         Vliq_hysteresis: tuple = (0.03, 0.01), q_int_W: float = 0.0, k_flash: float = 1.0):
    """
    Switches to gas-only discharge and PR EOS when gas volume fraction >= 99%.
    """
    # Internal energies
    U_liq = st.m_liq * liquid.u_J_kg(st.T_liq)
    U_gas = st.m_gas * vapor.u_J_kg(st.T_gas)

    # Geometries with current temperatures
    rho_l_curr = liquid.rho_kg_m3(st.T_liq)
    V_liq_curr = st.m_liq / max(rho_l_curr, 1e-12)
    V_gas_curr = max(cfg.V_tank - V_liq_curr, 1e-12)
    gas_frac_vol = V_gas_curr / cfg.V_tank

    # Mode decision
    gas_mode = gas_frac_vol >= 0.99

    if not gas_mode:
        # -------------- LIQUID-ONLY MODE --------------
        # Pressure and liquid outflow
        P = liquid.Psat_Pa(st.T_liq)
        m_dot_out_l = mass_flow_orifice_liquid(P, st.T_liq, cfg, liquid, k_flash=k_flash)

        # Apply discharge before evaporation
        dm_out = min(m_dot_out_l * dt, st.m_liq * 0.5)
        m_liq_minus = max(st.m_liq - dm_out, 1e-9)
        U_liq_minus = U_liq - dm_out * liquid.h_J_kg(st.T_liq)

        # Provisional volumes
        rho_l_guess = liquid.rho_kg_m3(st.T_liq)
        V_liq_guess = m_liq_minus / max(rho_l_guess, 1e-12)
        V_gas_guess = max(cfg.V_tank - V_liq_guess, 1e-12)
        Vliq_frac = V_liq_guess / cfg.V_tank

        suppress_evap = Vliq_frac <= Vliq_hysteresis[0]

        # Latent heat at current T
        L_curr = vapor.h_J_kg(st.T_liq) - liquid.h_J_kg(st.T_liq)

        # Evap cap
        dm_evap_cap = clip_frac * m_liq_minus
        # Required gas mass to fill ullage at saturation (use T_liq as interface T)
        rho_g_sat = vapor.rho_kg_m3(st.T_liq)
        m_g_req = rho_g_sat * V_gas_guess
        dm_need = m_g_req - st.m_gas
        alpha_use = 0.5 * alpha_evap if suppress_evap else alpha_evap
        dm_evap_relaxed = np.clip(alpha_use * dm_need, -dm_evap_cap, dm_evap_cap)

        # Update masses
        m_liq_new = max(m_liq_minus - dm_evap_relaxed, 1e-9)
        m_gas_new = max(st.m_gas + dm_evap_relaxed, 1e-9)

        # Energy transfer (latent heat only)
        U_liq_new = U_liq_minus - dm_evap_relaxed * L_curr
        U_gas_new = U_gas + dm_evap_relaxed * L_curr

        # Invert to temperature
        T_liq_new = invert_u_to_T(liquid, U_liq_new / m_liq_new, T_guess=st.T_liq)
        T_gas_new = invert_u_to_T(vapor, U_gas_new / m_gas_new, T_guess=st.T_gas)

        # Soft saturation projection
        rho_l = liquid.rho_kg_m3(T_liq_new)
        rho_g = vapor.rho_kg_m3(T_liq_new)
        V_liq_new = m_liq_new / max(rho_l, 1e-12)
        V_gas_new = max(cfg.V_tank - V_liq_new, 1e-12)
        m_g_req_fin = rho_g * V_gas_new
        delta_m = m_g_req_fin - m_gas_new

        m_liq_new = max(m_liq_new - alpha_sat * delta_m, 1e-9)
        m_gas_new = max(m_gas_new + alpha_sat * delta_m, 1e-9)
        U_liq_new = U_liq_new - alpha_sat * delta_m * L_curr
        U_gas_new = U_gas_new + alpha_sat * delta_m * L_curr

        # Final temps
        T_liq_new = invert_u_to_T(liquid, U_liq_new / m_liq_new, T_guess=T_liq_new)
        T_gas_new = invert_u_to_T(vapor, U_gas_new / m_gas_new, T_guess=T_gas_new)

        # Pressure from saturation
        P_new = liquid.Psat_Pa(T_liq_new)
        m_dot_out = m_dot_out_l

    else:
        # -------------- GAS-ONLY MODE + PR EOS --------------
        # Pressure from PR EOS using current gas state
        P = PR_pressure_from_state(st.m_gas, st.T_gas, V_gas_curr)

        # Gas outflow using upstream PR density
        m_dot_out_g = mass_flow_orifice_gas_PR(P, st.T_gas, cfg)

        # Apply gas discharge (no evaporation enforced; liquid may be residual)
        dm_out = min(m_dot_out_g * dt, st.m_gas * 0.5)
        m_gas_new = max(st.m_gas - dm_out, 1e-9)
        U_gas_new = U_gas - dm_out * vapor.u_J_kg(st.T_gas)

        # Keep liquid unchanged (optional: small re-equilibration with latent heat can be added)
        m_liq_new = st.m_liq
        U_liq_new = U_liq

        # Temperatures from energy
        T_gas_new = invert_u_to_T(vapor, U_gas_new / m_gas_new, T_guess=st.T_gas)
        T_liq_new = invert_u_to_T(liquid, U_liq_new / m_liq_new, T_guess=st.T_liq)

        # New volumes and PR pressure
        rho_l = liquid.rho_kg_m3(T_liq_new)
        V_liq_new = m_liq_new / max(rho_l, 1e-12)
        V_gas_new = max(cfg.V_tank - V_liq_new, 1e-12)
        P_new = PR_pressure_from_state(m_gas_new, T_gas_new, V_gas_new)
        m_dot_out = m_dot_out_g

    # Update state
    st.m_liq, st.m_gas = m_liq_new, m_gas_new
    st.T_liq, st.T_gas = T_liq_new, T_gas_new

    return st, P_new, m_dot_out, V_gas_new

# ----------------------------
# 初期条件設定関数
# ----------------------------
def find_T_for_pressure(liquid: PhaseTable, P_target: float, T_guess: float = 290.0):
    """
    飽和圧 Psat(T) = P_target となる温度を数値的に探索
    """
    from scipy.optimize import fsolve
    f = lambda T: liquid.Psat_Pa(T) - P_target
    T_sol = fsolve(f, T_guess)
    return float(T_sol[0])

# ----------------------------
# 実行例（液相100%からブローダウン）
# ----------------------------
if __name__ == "__main__":
    # CSVファイルから物性テーブルを読み込む
    liquid_csv = "density_dependent/N2O_liquid.csv"
    vapor_csv = "density_dependent/N2O_vapor.csv"
    liquid = load_phase_csv(liquid_csv, is_liquid=True)
    vapor = load_phase_csv(vapor_csv, is_liquid=False)

    # ユーザ指定パラメータ
    V_tank = 0.002      # タンク容積 [m3]
    P_init = 4.2e6        # 初期内圧 [Pa]

    # 初期液温度を飽和圧から逆算
    T_init = find_T_for_pressure(liquid, P_init, T_guess=290.0)

    V_liq = 0.99 * V_tank
    V_gas = 0.01 * V_tank

    # タンク設定
    cfg = TankConfig(
        V_tank=V_tank,
        A_orifice=np.pi/4 * (0.004**2),  # オリフィス径4mm相当
        Cd=0.333,
        T_wall=273,
        P_down=2e6,                      # 下流圧 [Pa]
        R_g=R_univ / M_N2O,
        gamma_g=1.25,
        P_down_gas=1e5                   # 燃焼終了後下流圧を大気圧
    )

    # 初期状態：液相100%充填
    st = TankState(
        T_liq=T_init,
        T_gas=T_init,
        m_liq=V_liq * liquid.rho_kg_m3(T_init),
        m_gas=V_gas * vapor.rho_kg_m3(T_init)
    )


    # 時間発展パラメータ
    dt = 0.001
    t_end = 20.0
    t = 0.0

    # ログ
    log_t, log_P, log_Tl, log_Tg, log_mL, log_mG, log_mout = [], [], [], [], [], [], []

    while t < t_end:
        st, P, m_dot_out, V_gas = step(cfg, st, dt, liquid, vapor,
                                       alpha_evap=0.3, alpha_sat=0.5, clip_frac=0.01,
                                       Vliq_hysteresis=(0.03, 0.01), k_flash=1.0)

        t += dt
        log_t.append(t); log_P.append(P)
        log_Tl.append(st.T_liq); log_Tg.append(st.T_gas)
        log_mL.append(st.m_liq); log_mG.append(st.m_gas)
        log_mout.append(m_dot_out)

        gas_fraction = V_gas / V_tank
        mode = "GAS" if gas_fraction >= 0.99 else "LIQ"
        if int(t/dt) % int(0.1/dt) == 0:
            print(f"[{mode}] t={t:.2f}s, P={P/1e5:.2f} bar, T_liq={st.T_liq:.2f} K, "
                  f"T_gas={st.T_gas:.2f} K, m_liq={st.m_liq:.4f} kg, m_gas={st.m_gas:.4f} kg, "
                  f"gas_vol={gas_fraction*100:.1f}%, m_dot={m_dot_out:.4f} kg/s")

        # Optional: stop if almost fully gas AND gas mass very low
        if gas_fraction >= 0.99 and st.m_gas < 1e-5:
            print(">>> Gas-only phase nearly depleted, stopping.")
            break


    # 結果表示
    fig, axs = plt.subplots(1, 2, figsize=(8, 10))

    axs[0].plot(log_t, np.array(log_P)/1e5)
    axs[0].set_ylabel("Tank Pressure [bar]")

    axs[1].plot(log_t, log_Tl, label="Liquid T")
    axs[1].plot(log_t, log_Tg, label="Gas T")
    axs[1].set_ylabel("Temperature [K]")
    axs[1].legend()
    plt.xlabel("Time [s]")
    plt.tight_layout()
    plt.show()

    fig, axs = plt.subplots(1, 2, figsize=(8, 10))
    axs[0].plot(log_t, log_mout)
    axs[0].set_ylabel("Mass flow rate [kg/s]")

    axs[1].plot(log_t, log_mL, label="Liquid mass")
    axs[1].plot(log_t, log_mG, label="Gas mass")
    axs[1].set_ylabel("Mass [kg]")
    axs[1].legend()

    plt.xlabel("Time [s]")
    plt.tight_layout()
    plt.show()