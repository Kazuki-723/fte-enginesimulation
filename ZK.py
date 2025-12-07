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
M_N2O = 44.013e-3  # kg/mol
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
    muJT_K_per_MPa: np.ndarray = None

    # スプライン補間器
    f_P_MPa: CubicSpline = None
    f_rho: CubicSpline = None
    f_u_kJ_per_mol: CubicSpline = None
    f_h_kJ_per_mol: CubicSpline = None
    f_cp_J_per_molK: CubicSpline = None
    f_visc_uPa_s: CubicSpline = None
    f_k_W_mK: CubicSpline = None
    f_sigma: CubicSpline = None
    f_muJT: CubicSpline = None

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
        if self.muJT_K_per_MPa is not None:
            self.f_muJT = CubicSpline(self.T, self.muJT_K_per_MPa, bc_type="natural")


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
    
    def muJT_K_per_MPa_at(self, Tq):
        if self.f_muJT is None:
            return 0.0
        return float(self.f_muJT(self._clip(Tq)))

def load_phase_csv(path: str, is_liquid: bool) -> PhaseTable:
    T, P_MPa, rho = [], [], []
    u_kJ_mol, h_kJ_mol, cp_J_molK = [], [], []
    visc_uPa_s, k_W_mK, sigma_N_m = [], [], []
    muJT_K_MPa = []

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
            if "Joule-Thomson (K/MPa)" in row and row["Joule-Thomson (K/MPa)"] != "":
                muJT_K_MPa.append(float(row["Joule-Thomson (K/MPa)"]))
            else:
                muJT_K_MPa.append(0.0)

    table = PhaseTable(
        T=np.array(T), P_MPa=np.array(P_MPa), rho=np.array(rho),
        u_kJ_per_mol=np.array(u_kJ_mol), h_kJ_per_mol=np.array( h_kJ_mol),
        cp_J_per_molK=np.array(cp_J_molK), viscosity_uPa_s=np.array(visc_uPa_s),
        k_W_mK=np.array(k_W_mK),
        surf_tension_N_m=(np.array(sigma_N_m) if sigma_N_m else None),
        muJT_K_per_MPa=np.array(muJT_K_MPa)
    )
    table.build(is_liquid=is_liquid)
    return table
# ----------------------------
# Peng–Robinson EOS utilities for N2O
# ----------------------------

# N2O constants
# https://ntrs.nasa.gov/api/citations/20190001326/downloads/20190001326.pdf
Tc_N2O = 309.6       # K
Pc_N2O = 7.25e6      # Pa
omega_N2O = 0.160     # acentric factor
kappa = 0.37464 + 1.54226 * omega_N2O - 0.26992 * omega_N2O**2 # Peng-robinson EoS param.

def PR_params_N2O(T: float):
    # PR parameters (a, b, alpha) for N2O at temperature T
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
    C_down: float

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

# -------------------------------
# Joule-Thomson Cooling effect
# -------------------------------

def apply_joule_thomson_cooling(vapor: PhaseTable, T_gas: float, P_old: float, P_new: float):
    """
    μJT(T) [K/MPa] を用い、ΔT_JT = μJT(T_avg) * ΔP(MPa) を適用。
    圧力低下（ΔP<0）で T が低下（μJT>0の場合）。
    """
    T_avg = T_gas
    muJT = vapor.muJT_K_per_MPa_at(T_avg)
    dP_MPa = (P_old - P_new) / 1e6
    return T_gas + muJT * dP_MPa

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

def blended_pressure(T_liq, T_gas, m_gas, V_gas, Vliq_frac, liquid, cfg):
    P_sat = liquid.Psat_Pa(T_liq)
    P_pr  = PR_pressure_from_state(m_gas, T_gas, V_gas)
    # ブレンド開始・終了閾値
    f_start, f_end = 0.1, 0.01  # 液体体積比
    if Vliq_frac > f_start:
        return P_sat
    elif Vliq_frac < f_end:
        return P_pr
    else:
        # 線形補間
        w = (f_start - Vliq_frac) / (f_start - f_end)
        return (1-w)*P_sat + w*P_pr

# 内部エネルギーから温度を逆算する関数
def invert_u_to_T(table: PhaseTable, u_target: float, T_guess: float):
    f = lambda T: table.u_J_kg(T) - u_target
    T_sol = fsolve(f, T_guess)
    return float(T_sol[0])

def energy_init(cfg: TankConfig, st: TankState):
    U_liq = st.m_liq * liquid.u_J_kg(st.T_liq)
    U_gas = st.m_gas * vapor.u_J_kg(st.T_gas)
    return U_liq + U_gas

def step(cfg: TankConfig, st: TankState, dt: float, liquid: PhaseTable, vapor: PhaseTable,
         alpha_evap: float, alpha_sat: float, clip_frac: float,
         Vliq_hysteresis: tuple, k_flash: float, time: float, u_total: float, q_int_W: float = 0.0,
         h_int: float = 200.0, A_int: float = 0.003848, beta_mix: float = 0.4):
    """
    h_int: 界面伝熱係数 [W/m2/K]
    A_int: 界面面積 [m2]
    beta_mix: ガス温度の界面整合リラクゼーション係数 (0<β≤1)
    """

    # 現在の体積比
    rho_l_curr = liquid.rho_kg_m3(st.T_liq)
    V_liq_curr = st.m_liq / max(rho_l_curr, 1e-12)
    V_gas_curr = max(cfg.V_tank - V_liq_curr, 1e-12)
    gas_mode = V_gas_curr / cfg.V_tank >= 0.99

    # 内部エネルギー
    #U_liq = st.m_liq * liquid.u_J_kg(st.T_liq)
    U_gas = st.m_gas * vapor.u_J_kg(st.T_gas)
    U_liq = u_total - U_gas

    if not gas_mode:
        # ---------- 液相排出モード ----------
        P_old = liquid.Psat_Pa(st.T_liq)
        P = liquid.Psat_Pa(st.T_liq)
        cfg.P_down = C_down * np.sqrt(P)
        m_dot_out_l = mass_flow_orifice_liquid(P, st.T_liq, cfg, liquid, k_flash=k_flash)

        dm_out = min(m_dot_out_l * dt, st.m_liq * 0.5)
        m_liq_minus = max(st.m_liq - dm_out, 1e-9)
        U_liq_minus = U_liq - dm_out * liquid.h_J_kg(st.T_liq)

        # 蒸発の緩和
        rho_l_guess = liquid.rho_kg_m3(st.T_liq)
        V_liq_guess = m_liq_minus / max(rho_l_guess, 1e-12)
        V_gas_guess = max(cfg.V_tank - V_liq_guess, 1e-12)
        suppress_evap = (V_liq_guess / cfg.V_tank) <= Vliq_hysteresis[0]

        #L_curr = (18.535262 - st.T_liq * 0.011019) * 1000 / M_N2O
        #L_curr = 16.5 * ((1 - st.T_liq/Tc_N2O) / (1 - 184.7/Tc_N2O)) ** 0.38 * 1000 / M_N2O
        L_curr = vapor.h_J_kg(st.T_liq) - liquid.h_J_kg(st.T_liq)
        #print(L_curr)
        dm_evap_cap = clip_frac * m_liq_minus
        rho_g_sat = vapor.rho_kg_m3(st.T_liq)
        m_g_req = rho_g_sat * V_gas_guess
        dm_need = m_g_req - st.m_gas
        alpha_use = 0.5 * alpha_evap if suppress_evap else alpha_evap
        dm_evap_relaxed = np.clip(alpha_use * dm_need, -dm_evap_cap, dm_evap_cap)

        # 質量
        m_liq_new = max(m_liq_minus - dm_evap_relaxed, 1e-9)
        m_gas_new = max(st.m_gas + dm_evap_relaxed, 1e-9)

        # 潜熱のやり取り
        Latent_heat =  dm_evap_relaxed * L_curr
        U_liq_new = U_liq_minus - Latent_heat
        U_gas_new = U_gas + Latent_heat

        # 温度逆写像（一次）
        T_liq_new = invert_u_to_T(liquid, U_liq_new / m_liq_new, T_guess=st.T_liq)
        T_gas_new = invert_u_to_T(vapor, U_gas_new / m_gas_new, T_guess=st.T_gas)

        # 界面の感熱整合（ガス側へ Q_int_g を加え、温度もリラクゼーション）
        Q_int_g = h_int * A_int * (T_liq_new - T_gas_new)
        U_liq_new -= Q_int_g * dt
        U_gas_new += Q_int_g * dt
        # ガス温度の再逆写像
        T_gas_new = invert_u_to_T(vapor, U_gas_new / m_gas_new, T_guess=T_gas_new)
        # 温度リラクゼーション（数値安定化）
        T_gas_new = T_gas_new + beta_mix * (T_liq_new - T_gas_new)

        # 飽和投影の緩和
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

        # 最終温度
        T_liq_new = invert_u_to_T(liquid, U_liq_new / m_liq_new, T_guess=T_liq_new)
        T_gas_new = invert_u_to_T(vapor, U_gas_new / m_gas_new, T_guess=T_gas_new)
        # 温度リラクゼーション
        T_gas_new = T_gas_new + beta_mix * (T_liq_new - T_gas_new)

        # 圧力
        #Vliq_frac_new = V_liq_new / cfg.V_tank
        #P_new = blended_pressure(T_liq_new, T_gas_new, m_gas_new, V_gas_new, Vliq_frac_new, liquid, cfg)
        P_new = PR_pressure_from_state(m_gas_new, T_gas_new, V_gas_new)
        m_dot_out = m_dot_out_l

        # μJTによる補正（圧力変化に応じて）
        T_gas_corr = apply_joule_thomson_cooling(vapor, T_gas_new, P_old, P_new)

        # 内部エネルギーを整合（u(T)へ合わせる）
        U_gas_new = m_gas_new * vapor.u_J_kg(T_gas_corr)
        T_gas_new = T_gas_corr
        U_out = dm_out * liquid.u_J_kg(T_liq_new)

    else:
        # ---------- 気相排出モード（PR EOS） ----------
        # PR圧力
        P = PR_pressure_from_state(st.m_gas, st.T_gas, V_gas_curr)
        # ガス流量
        m_dot_out_g = mass_flow_orifice_gas_PR(P, st.T_gas, cfg)

        dm_out = min(m_dot_out_g * dt, st.m_gas * 0.5)
        m_gas_new = max(st.m_gas - dm_out, 1e-9)
        U_gas_new = U_gas - dm_out * vapor.u_J_kg(st.T_gas)

        # 液はそのまま（微小なら無視）
        m_liq_new = st.m_liq
        U_liq_new = U_liq
        Latent_heat = 0

        # 温度
        T_gas_new = invert_u_to_T(vapor, U_gas_new / m_gas_new, T_guess=st.T_gas)
        T_liq_new = invert_u_to_T(liquid, U_liq_new / m_liq_new, T_guess=st.T_liq)

        # 体積・PR圧力
        rho_l = liquid.rho_kg_m3(T_liq_new)
        V_liq_new = m_liq_new / max(rho_l, 1e-12)
        V_gas_new = max(cfg.V_tank - V_liq_new, 1e-12)
        P_new = PR_pressure_from_state(m_gas_new, T_gas_new, V_gas_new)
        m_dot_out = m_dot_out_g
        U_out = dm_out * vapor.u_J_kg(T_gas_new)

    # 状態更新
    st.m_liq, st.m_gas = m_liq_new, m_gas_new
    st.T_liq, st.T_gas = T_liq_new, T_gas_new
    u_total -= U_out
    properties ={'U_liq': U_liq_new, 
                 'U_gas': U_gas_new, 
                 'U_out': U_out,
                 'rho_liq': liquid.rho_kg_m3(T_liq_new), 
                 'rho_gas': vapor.rho_kg_m3(T_gas_new),
                 'Latent_Heat': Latent_heat
                }
    return st, P_new, m_dot_out, V_gas_new, properties, u_total

# ----------------------------
# 初期条件設定関数
# ----------------------------
def find_T_for_pressure(liquid: PhaseTable, P_target: float, T_guess: float):
    """
    飽和圧 Psat(T) = P_target となる温度を数値的に探索
    """
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
    V_tank = 0.00044      # タンク容積 [m3]
    P_init = 4.2e6        # 初期内圧 [Pa]
    P_down_init = 2e6     # 初期燃焼室圧[Pa]
    C_down = P_down_init / np.sqrt(P_init)

    # 初期液温度を飽和圧から逆算
    T_init = find_T_for_pressure(liquid, P_init, T_guess=250.0)

    V_liq = 0.99 * V_tank
    V_gas = 0.01 * V_tank

    # タンク設定
    cfg = TankConfig(
        V_tank=V_tank,
        A_orifice=np.pi/4 * (0.003**2),  # オリフィス径4mm相当
        Cd=0.333,
        T_wall=200,
        P_down=P_down_init,                      # 下流圧 [Pa]
        R_g=R_univ / M_N2O,
        gamma_g=1.30,
        P_down_gas=1e5,                   # 燃焼終了後下流圧を大気圧
        C_down = C_down
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
    i = 0

    # ログ
    log_t, log_P, log_P_down, log_Tl, log_Tg, log_mL, log_mG, log_mout = [], [], [], [], [], [], [], []
    log_u_liq, log_u_gas, log_u_out, log_latent_heat = [], [], [], []
    u_out_tot = 0
    latent_heat_tot = 0
    log_u_out_tot, log_latent_heat_tot, log_u_tot = [], [], []

    u_total = energy_init(cfg, st)
    while t < t_end:
        st, P, m_dot_out, V_gas, properties, u_total = step(cfg, st, dt, liquid, vapor,
                                       alpha_evap=0.2, alpha_sat=0.5, clip_frac=0.01, u_total = u_total,
                                       Vliq_hysteresis=(0.05, 0.01), k_flash=1,time = t)

        t += dt
        gas_fraction = V_gas / V_tank
        log_t.append(t); log_P.append(P)
        log_Tl.append(st.T_liq); log_Tg.append(st.T_gas)
        log_mL.append(st.m_liq); log_mG.append(st.m_gas)
        log_mout.append(m_dot_out)
        log_u_liq.append(properties["U_liq"]); log_u_gas.append(properties["U_gas"]); log_u_out.append(properties["U_out"]); log_latent_heat.append(properties["Latent_Heat"])
        u_out_tot += properties["U_out"]
        latent_heat_tot += properties["Latent_Heat"]
        log_u_out_tot.append(u_out_tot); log_latent_heat_tot.append(latent_heat_tot)
        log_u_tot.append(properties["U_liq"] + properties["U_gas"] + u_out_tot)
        if gas_fraction >= 0.99:
            log_P_down.append(cfg.P_down_gas)
        else:
            log_P_down.append(cfg.P_down)

        mode = "GAS" if gas_fraction >= 0.99 else "LIQ"
        if int(t/dt) % int(0.05/dt) == 0:
            print(f"[{mode}] t={t:.2f}s, P={P/1e5:.2f} bar, T_liq={st.T_liq:.2f} K, "
                  f"T_gas={st.T_gas:.2f} K, m_liq={st.m_liq:.4f} kg, m_gas={st.m_gas:.4f} kg, "
                  f"U_liq={properties["U_liq"]/1e3:.4f} kJ, U_gas={properties["U_gas"]/1e3:.4f} kJ, U_out={properties["U_out"]/1e3:.4f} kJ,"
                  f"Latent Heat={properties["Latent_Heat"]/1e3:.4f} kJ, "
                  f"gas_vol={gas_fraction*100:.1f}%, m_dot={m_dot_out:.4f} kg/s")

        # liquid only mode
        # if gas_fraction >= 0.99:
        #     print(">>> Liquid-only phase nearly depleted.")
        #     break

        # Optional: stop if almost fully gas AND gas mass very low
        if gas_fraction >= 0.99 and m_dot_out < 1e-5:
            print(">>> Gas-only phase nearly depleted, stopping.")
            break


    # 結果表示
    fig, axs = plt.subplots(1, 2, figsize=(8, 10))

    axs[0].plot(log_t, np.array(log_P)/1e5)
    axs[0].plot(log_t, np.array(log_P_down)/1e5)
    axs[0].set_ylabel("Tank Pressure [bar]")
    axs[0].set_ylim(0, 50)

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

    plt.figure(figsize=(10, 5))
    plt.plot(log_t, log_u_liq, label="liquid")
    plt.plot(log_t, log_u_gas, label="gas")
    plt.plot(log_t, log_u_out_tot, label="out")
    plt.plot(log_t, log_latent_heat_tot, label="latent heat")
    plt.plot(log_t, log_u_tot, label="total")
    plt.xlabel("Time [s]")
    plt.ylabel("Internal Energy [J]")
    plt.title("Internal Energy Over Time")
    #plt.ylim(65000,75000)

    plt.legend()
    plt.tight_layout()
    plt.show()

    # plt.figure(figsize=(10, 5))
    # plt.plot(log_t, log_latent_heat_tot, label="latent heat")
    # plt.xlabel("Time [s]")
    # plt.ylabel("Internal Energy [J]")
    # plt.title("Internal Energy Over Time")

    # plt.legend()
    # plt.tight_layout()
    # plt.show()