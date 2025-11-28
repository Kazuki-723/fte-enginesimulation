import csv
import math
import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Dict
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

# ----------------------------
# Utilities: linear and monotone interpolation with clamping
# ----------------------------

def clamp(x, xmin, xmax):
    return max(min(x, xmax), xmin)

class TableInterp1D:
    """
    1D interpolator keyed by Temperature (K).
    Assumes strictly increasing Temperature values.
    Clamps to endpoints outside the range.
    """
    def __init__(self, T: List[float], Y: List[float]):
        assert len(T) == len(Y) and len(T) >= 2
        # sort by temperature
        pairs = sorted(zip(T, Y), key=lambda p: p[0])
        self.T = [p[0] for p in pairs]
        self.Y = [p[1] for p in pairs]

    def __call__(self, Tq: float) -> float:
        Tq = clamp(Tq, self.T[0], self.T[-1])
        # binary search
        lo, hi = 0, len(self.T) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if self.T[mid] <= Tq:
                lo = mid
            else:
                hi = mid
        # linear interpolate
        t0, t1 = self.T[lo], self.T[hi]
        y0, y1 = self.Y[lo], self.Y[hi]
        if t1 == t0:
            return y0
        w = (Tq - t0) / (t1 - t0)
        return y0 * (1.0 - w) + y1 * w


# ----------------------------
# CSV loaders and property model from tables
# ----------------------------

M_N2O = 44.0128e-3  # kg/mol

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

    # 温度範囲外は端点でクリップ
    def _clip(self, Tq):
        return np.clip(Tq, self.T[0], self.T[-1])

    # 物性取得関数
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
# EQ model for liquid-only discharge using table properties
# ----------------------------

class EQTankModel:
    def __init__(self, V_tank_m3, liquid, vapor,
                 A_orifice_m2, Cd, C_down,
                 T_wall_K=None, h_wall_Wm2K=0.0, A_wall_m2=0.0):
        self.V = V_tank_m3
        self.liq = liquid
        self.vap = vapor
        self.A_or = A_orifice_m2
        self.Cd = Cd
        self.C_down = C_down
        self.T_wall = T_wall_K
        self.h_wall = h_wall_Wm2K
        self.A_wall = A_wall_m2

    def Psat(self, T):
        return self.liq.Psat_Pa(T)

    def Q_wall(self, T_l):
        if self.h_wall <= 0.0 or self.A_wall <= 0.0 or self.T_wall is None:
            return 0.0
        return self.h_wall * self.A_wall * (self.T_wall - T_l)

    def algebraic_state(self, m_l, T):
        rho_l = self.liq.rho_kg_m3(T)
        V_v = max(self.V - m_l / rho_l, 0.0)
        P = self.Psat(T)
        rho_v = self.vap.rho_kg_m3(T)
        m_v = rho_v * V_v
        return V_v, P, m_v

    def flow_liquid(self, T, Pt):
        rho_l = self.liq.rho_kg_m3(T)
        P_down = max(self.C_down * math.sqrt(Pt), 1.0e5)
        dp = max(Pt - P_down, 0.0)
        return self.Cd * self.A_or * rho_l * math.sqrt(2.0 * dp / rho_l)

    def step(self, m_l, T, dt, max_iter=10, tol=1e-6):
        """
        EQモデル逐次収束ステップ
        """
        # 初期推定
        Vv0, Pt, mv0 = self.algebraic_state(m_l, T)

        # 液相流出量
        m_out = self.flow_liquid(T, Pt)
        m_l_pred = max(m_l - m_out * dt, 0.0)

        # 逐次収束ループ
        T_guess = T
        for _ in range(max_iter):
            # 平衡拘束から蒸気量を計算
            Vv, P, m_v = self.algebraic_state(m_l_pred, T_guess)

            # エネルギー収支
            Qw = self.Q_wall(T_guess)
            h_l = self.liq.h_J_kg(T_guess)
            u_l = self.liq.u_J_kg(T_guess)
            h_v = self.vap.h_J_kg(T_guess)
            cp_l = self.liq.cp_J_kgK(T_guess)

            # 蒸発量を体積拘束から算出
            m_ph = (m_v - mv0) / dt

            rhs = Qw - m_out * (h_l - u_l) - m_ph * (h_v - h_l)
            dT = rhs / (max(m_l_pred * cp_l, 1e-9)) * dt

            # 更新
            T_new = np.clip(T_guess + dT, self.liq.T[0], self.liq.T[-1])

            if abs(T_new - T_guess) < tol:
                T_guess = T_new
                break
            T_guess = T_new

        # 最終状態
        VvN, PtN, mvN = self.algebraic_state(m_l_pred, T_guess)

        diag = {
            "P_tank_Pa": PtN,
            "P_down_Pa": self.C_down * math.sqrt(PtN),
            "T_K": T_guess,
            "m_l_kg": m_l_pred,
            "m_v_kg": mvN,
            "m_out_kg_s": m_out,
            "m_ph_kg_s": m_ph,
            "V_v": VvN,
            "rho_l_kg_m3": self.liq.rho_kg_m3(T)
        }
        return m_l_pred, T_guess, diag

# ----------------------------
# Initialization: find T from target Psat using table
# ----------------------------

def find_T_from_Psat_table(target_P_Pa: float, table: PhaseTable) -> float:
    # Invert Psat(T) by scanning and interpolating
    Tvals = table.T
    Pvals = [table.f_P_MPa(T) * 1e6 for T in Tvals]
    # clamp target within range
    target_P_Pa = clamp(target_P_Pa, min(Pvals), max(Pvals))
    # find bracket
    for i in range(len(Tvals) - 1):
        P0, P1 = Pvals[i], Pvals[i + 1]
        if (P0 - target_P_Pa) * (P1 - target_P_Pa) <= 0.0:
            # linear inverse
            w = (target_P_Pa - P0) / (P1 - P0) if P1 != P0 else 0.0
            return Tvals[i] * (1 - w) + Tvals[i + 1] * w
    # fallback: nearest
    j = min(range(len(Pvals)), key=lambda k: abs(Pvals[k] - target_P_Pa))
    return Tvals[j]


# ----------------------------
# Example run configuration
# ----------------------------

if __name__ == "__main__":
    # Paths to your CSVs
    liquid_csv = "density_dependent\\N2O_liquid.csv"
    vapor_csv = "density_dependent\\N2O_vapor.csv"

    # Load tables
    liquid = load_phase_csv(liquid_csv, is_liquid=True)
    vapor = load_phase_csv(vapor_csv, is_liquid=False)

    # 初期条件
    T0 = 285.2  # K
    Pc = 2e6
    V_tank_m3 = 0.002
    rho_l0 = liquid.rho_kg_m3(T0)
    m_l0 = rho_l0 * V_tank_m3
    Pt = liquid.Psat_Pa(T0)
    C_down = Pc / math.sqrt(Pt)
    

    model = EQTankModel(
        V_tank_m3= V_tank_m3,
        liquid=liquid,
        vapor=vapor,
        A_orifice_m2= math.pi / 4 * 0.004 ** 2,
        Cd=0.4,
        C_down=C_down
    )

    t_hist, P_hist, T_hist, m_out_hist = [], [], [], []
    m_l, T = m_l0, T0
    t, dt, T_end = 0.0, 1e-3, 10.0

    print(f"Initial: T={T0:.2f} K, Psat={liquid.Psat_Pa(T0)/1e6:.3f} MPa, "
            f"rho_l={rho_l0:.1f} kg/m^3, m_l={m_l0:.3f} kg")

    while t < T_end:
        m_l, T, diag = model.step(m_l, T, dt)
        t += dt
        t_hist.append(t)
        P_hist.append(diag["P_tank_Pa"]/1e6)
        T_hist.append(diag["T_K"])
        m_out_hist.append(diag["m_out_kg_s"])

        if diag["V_v"] / V_tank_m3 > 0.99:
            print("end simulation")
            break

        if int(t*1000) % 50 == 0:
            print(f"t={t:.2f} s | Pt={diag['P_tank_Pa']/1e6:.3f} MPa | Pc={diag['P_down_Pa']/1e6:.3f} MPa | T={T:.2f} K | "
                  f"m_l={m_l:.3f} kg | m_out={diag['m_out_kg_s']:.3f} kg/s | "
                  f"Vv={diag['V_v']*1e6:.1f} cm^3 | rho_l={diag['rho_l_kg_m3']:.1f}")
        
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 3, 1)
    plt.plot(t_hist, P_hist)
    plt.xlabel("Time [s]")
    plt.ylabel("Pressure [MPa]")
    plt.title("Tank, Chamber Pressure Over Time")

    plt.subplot(1, 3, 2)
    plt.plot(t_hist, m_out_hist)
    plt.xlabel("Time [s]")
    plt.ylabel("Mass flow rate")
    plt.title("mass flow rate Over Time")

    plt.subplot(1, 3, 3)
    plt.plot(t_hist,T_hist)
    plt.xlabel("Time [s]")
    plt.ylabel("Temperture[kg/s]")
    plt.title("temperture Over Time")


    plt.tight_layout()
    plt.show()