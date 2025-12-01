import numpy as np
import csv
from dataclasses import dataclass
from typing import List
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

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

def choked_mass_flow(P_tank, T_gas, cfg: TankConfig):
    R = cfg.R_g
    rho = P_tank / (R * T_gas)
    dP = max(P_tank - cfg.P_down, 0.0)
    return cfg.Cd * cfg.A_orifice * np.sqrt(2 * dP * rho)

# 内部エネルギーから温度を逆算する関数
def invert_u_to_T(table: PhaseTable, u_target: float, T_guess: float = 290.0):
    f = lambda T: table.u_J_kg(T) - u_target
    T_sol = fsolve(f, T_guess)
    return float(T_sol[0])

def step(cfg: TankConfig, st: TankState, dt: float, liquid: PhaseTable, vapor: PhaseTable):
    V_gas = ullage_volume(cfg, st, liquid)
    P = tank_pressure(st, liquid)

    # 吐出流量
    m_dot_out_g = choked_mass_flow(P, st.T_gas, cfg)

    # 飽和条件に合わせて蒸発量を調整
    rho_g_sat = vapor.rho_kg_m3(st.T_liq)
    m_g_required = rho_g_sat * V_gas
    dm_g_target = m_g_required - st.m_gas
    # 蒸発量を制限（液質量の数％/s）
    m_dot_evap = np.clip(dm_g_target/dt, -0.02*st.m_liq/dt, 0.02*st.m_liq/dt)

    # 内部エネルギー
    U_liq = st.m_liq * liquid.u_J_kg(st.T_liq)
    U_gas = st.m_gas * vapor.u_J_kg(st.T_gas)

    # エネルギー収支更新
    L = vapor.h_J_kg(st.T_liq) - liquid.h_J_kg(st.T_liq)  # 潜熱近似
    U_liq_new = U_liq - m_dot_evap * L * dt
    U_gas_new = U_gas + m_dot_evap * L * dt - m_dot_out_g * vapor.h_J_kg(st.T_gas) * dt

    # 質量更新
    m_liq_new = max(st.m_liq - m_dot_evap * dt, 1e-6)
    m_gas_new = max(st.m_gas + m_dot_evap * dt - m_dot_out_g * dt, 1e-9)

    # 温度更新（内部エネルギーから逆算）
    T_liq_new = invert_u_to_T(liquid, U_liq_new / m_liq_new, T_guess=st.T_liq)
    T_gas_new = invert_u_to_T(vapor, U_gas_new / m_gas_new, T_guess=st.T_gas)

    st.m_liq, st.m_gas = m_liq_new, m_gas_new
    st.T_liq, st.T_gas = T_liq_new, T_gas_new

    return st, P, m_dot_out_g


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
        Cd=0.4,
        T_wall=293.15,
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
    dt = 0.01
    t_end = 50.0
    t = 0.0

    # ログ
    log_t, log_P, log_Tl, log_Tg, log_mL, log_mG, log_mout = [], [], [], [], [], [], []

    while t < t_end and st.m_liq > 0.05:
        st, P, m_dot_out = step(cfg, st, dt, liquid, vapor)
        t += dt
        log_t.append(t)
        log_P.append(P)
        log_Tl.append(st.T_liq)
        log_Tg.append(st.T_gas)
        log_mL.append(st.m_liq)
        log_mG.append(st.m_gas)
        log_mout.append(m_dot_out)

        total_mass = st.m_liq + st.m_gas
        gas_fraction = st.m_gas / total_mass

        # 一定時間ごとにprint（例：1秒ごと）
        if int(t/dt) % int(1.0/dt) == 0:
            print(f"t={t:.2f}s, P={P/1e5:.2f} bar, T_liq={st.T_liq:.2f} K, "
                  f"T_gas={st.T_gas:.2f} K, m_liq={st.m_liq:.3f} kg, "
                  f"m_gas={st.m_gas:.3f} kg, gas_frac={gas_fraction*100:.1f}%")

        # 気相が99%以上になったら停止
        if gas_fraction >= 0.99:
            print(">>> no liquid")
            break


    # 結果表示
    fig, axs = plt.subplots(3, 1, figsize=(8, 10))

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