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

def choked_mass_flow(P_tank, T_liquid, cfg: TankConfig):
    # R = cfg.R_g
    # rho = P_tank / (R * T_liquid)
    rho = liquid.rho_kg_m3(T_liquid)
    dP = max(P_tank - cfg.P_down, 0.0)
    return cfg.Cd * cfg.A_orifice * np.sqrt(2 * dP * rho)

# 内部エネルギーから温度を逆算する関数
def invert_u_to_T(table: PhaseTable, u_target: float, T_guess: float = 290.0):
    f = lambda T: table.u_J_kg(T) - u_target
    T_sol = fsolve(f, T_guess)
    return float(T_sol[0])

def step(cfg: TankConfig, st: TankState, dt: float, liquid: PhaseTable, vapor: PhaseTable,
         alpha_evap: float = 0.3, alpha_sat: float = 0.5, clip_frac: float = 0.01,
         Vliq_hysteresis: tuple = (0.03, 0.01), q_int_W: float = 0.0):
    """
    alpha_evap: 蒸発量のアンダーリラクゼーション係数 (0<α≤1)
    alpha_sat : 飽和投影の差分緩和係数
    clip_frac : ステップあたりの蒸発上限 (液質量に対する比率)
    Vliq_hysteresis: (抑制開始液体体積比, ガスモード切替液体体積比)
    q_int_W: 界面熱流束 [W]（既知なら指定、未知なら0で無視）
    """
    # 内部エネルギー
    U_liq = st.m_liq * liquid.u_J_kg(st.T_liq)
    U_gas = st.m_gas * vapor.u_J_kg(st.T_gas)

    # 圧力と液用流量
    P = liquid.Psat_Pa(st.T_liq)
    m_dot_out_l = choked_mass_flow(P, st.T_liq, cfg)

    # 排出（蒸発前）
    dm_out = m_dot_out_l * dt
    dm_out = min(dm_out, st.m_liq * 0.5)  # セーフティ
    m_liq_minus = max(st.m_liq - dm_out, 1e-9)
    U_liq_minus = U_liq - dm_out * liquid.h_J_kg(st.T_liq)

    # 暫定液体体積比（旧Tで概算）
    V_liq_tmp = m_liq_minus / max(liquid.rho_kg_m3(st.T_liq), 1e-9)
    V_liq_frac_tmp = V_liq_tmp / cfg.V_tank

    # ヒステリシス判定
    suppress_evap = V_liq_frac_tmp <= Vliq_hysteresis[0]
    gas_mode = V_liq_frac_tmp <= Vliq_hysteresis[1]

    # 潜熱
    L_curr = vapor.h_J_kg(st.T_liq) - liquid.h_J_kg(st.T_liq)

    # 蒸発上限（物理＋数値）
    dm_evap_max_clip = clip_frac * m_liq_minus
    dm_evap_max_flux = (q_int_W / max(L_curr, 1e-6)) * dt if q_int_W > 0.0 else dm_evap_max_clip
    dm_evap_cap = min(dm_evap_max_clip, dm_evap_max_flux)

    if gas_mode:
        # ガス単相モード：蒸発を止め、圧力はPsatではなくガスEoSで更新すべきだが、
        # ここでは簡易に蒸発ゼロで更新（液がほぼ消失）
        dm_evap_relaxed = 0.0
    else:
        # 飽和要件から必要ガス質量を、"新しい質量"で一貫再計算するためにT_intを仮置き
        T_int_guess = st.T_liq
        rho_l_guess = liquid.rho_kg_m3(T_int_guess)
        rho_g_guess = vapor.rho_kg_m3(T_int_guess)
        V_liq_new_guess = m_liq_minus / max(rho_l_guess, 1e-9)
        V_gas_new_guess = max(cfg.V_tank - V_liq_new_guess, 1e-9)
        m_g_required = rho_g_guess * V_gas_new_guess

        dm_evap_needed = m_g_required - st.m_gas
        if suppress_evap:
            # 抑制域では過補正を避けるためさらに緩和
            alpha_use = 0.5 * alpha_evap
        else:
            alpha_use = alpha_evap

        dm_evap_relaxed = np.clip(alpha_use * dm_evap_needed, -dm_evap_cap, dm_evap_cap)

    # 新質量
    m_liq_new = max(m_liq_minus - dm_evap_relaxed, 1e-9)
    m_gas_new = max(st.m_gas + dm_evap_relaxed, 1e-9)

    # 新しい界面温度を内部エネルギーで決定（潜熱のみ移動）
    U_liq_new = U_liq_minus - dm_evap_relaxed * L_curr
    U_gas_new = U_gas + dm_evap_relaxed * L_curr

    # 逆写像で新温度
    T_liq_new = invert_u_to_T(liquid, U_liq_new / m_liq_new, T_guess=st.T_liq)
    T_gas_new = invert_u_to_T(vapor, U_gas_new / m_gas_new, T_guess=st.T_gas)

    # 飽和投影（穏やかに近づける）
    rho_l = liquid.rho_kg_m3(T_liq_new)
    rho_g = vapor.rho_kg_m3(T_liq_new)
    V_liq_new = m_liq_new / max(rho_l, 1e-9)
    V_gas_new = max(cfg.V_tank - V_liq_new, 1e-9)
    m_g_required_fin = rho_g * V_gas_new
    delta_m = m_g_required_fin - m_gas_new

    # 過補正防止の緩和
    m_liq_new = max(m_liq_new - alpha_sat * delta_m, 1e-9)
    m_gas_new = max(m_gas_new + alpha_sat * delta_m, 1e-9)

    # エネルギーも対応
    U_liq_new = U_liq_new - alpha_sat * delta_m * L_curr
    U_gas_new = U_gas_new + alpha_sat * delta_m * L_curr

    # 最終温度
    T_liq_new = invert_u_to_T(liquid, U_liq_new / m_liq_new, T_guess=T_liq_new)
    T_gas_new = invert_u_to_T(vapor, U_gas_new / m_gas_new, T_guess=T_gas_new)

    # 圧力更新
    P_new = liquid.Psat_Pa(T_liq_new)

    # 状態更新
    st.m_liq, st.m_gas = m_liq_new, m_gas_new
    st.T_liq, st.T_gas = T_liq_new, T_gas_new

    return st, P_new, m_dot_out_l, V_gas_new

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
    V_tank = 0.00044      # タンク容積 [m3]
    P_init = 4.2e6        # 初期内圧 [Pa]

    # 初期液温度を飽和圧から逆算
    T_init = find_T_for_pressure(liquid, P_init, T_guess=290.0)

    V_liq = 0.99 * V_tank
    V_gas = 0.01 * V_tank

    # タンク設定
    cfg = TankConfig(
        V_tank=V_tank,
        A_orifice=np.pi/4 * (0.003**2),  # オリフィス径4mm相当
        Cd=0.333,
        T_wall=273,
        P_down=2e6,                      # 下流圧 [Pa]
        R_g=R_univ / M_N2O,
        gamma_g=1.25
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
    t_end = 10.0
    t = 0.0

    # ログ
    log_t, log_P, log_Tl, log_Tg, log_mL, log_mG, log_mout = [], [], [], [], [], [], []

    while t < t_end:
        st, P, m_dot_out, V_gas = step(cfg, st, dt, liquid, vapor)
        t += dt
        log_t.append(t)
        log_P.append(P)
        log_Tl.append(st.T_liq)
        log_Tg.append(st.T_gas)
        log_mL.append(st.m_liq)
        log_mG.append(st.m_gas)
        log_mout.append(m_dot_out)

        total_mass = st.m_liq + st.m_gas
        gas_fraction = V_gas / V_tank

        # 一定時間ごとにprint
        if int(t/dt) % int(0.1/dt) == 0:
            print(f"t={t:.2f}s, P={P/1e5:.2f} bar, T_liq={st.T_liq:.2f} K, "
                  f"T_gas={st.T_gas:.2f} K, m_liq={st.m_liq:.3f} kg, "
                  f"m_gas={st.m_gas:.3f} kg, gas_frac={gas_fraction*100:.1f}%,"
                  f"m_dot={m_dot_out:.3f} kg/s")

        # 気相が99%以上になったら停止
        if gas_fraction >= 0.99:
            print(">>> no liquid")
            print(f"end time{t:.3f} s, end pressure{np.min(P)/1e6:.3f} MPa")
            break


    # 結果表示
    fig, axs = plt.subplots(1, 3, figsize=(8, 10))

    axs[0].plot(log_t, np.array(log_P)/1e5)
    axs[0].set_ylabel("Tank Pressure [bar]")

    axs[1].plot(log_t, log_Tl, label="Liquid T")
    axs[1].plot(log_t, log_Tg, label="Gas T")
    axs[1].set_ylabel("Temperature [K]")
    axs[1].legend()

    axs[2].plot(log_t, log_mL, label="Liquid mass")
    axs[2].plot(log_t, log_mG, label="Gas mass")
    axs[2].set_ylabel("Mass [kg]")
    axs[2].legend()

    plt.xlabel("Time [s]")
    plt.tight_layout()
    plt.show()