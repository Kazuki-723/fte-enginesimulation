import re
import csv
from inputprograms.rocket_simulation import RocketSimulation
from inputprograms.interp_density import OxidizerDatabase
from inputprograms.importjson import JsoncLoader
sim = RocketSimulation()
ox_db = OxidizerDatabase()

def parse_initial_results(text: str) -> dict:
    result = {}
    # K*
    match_k = re.search(r"K\* *= *([\d\.Ee+-]+)", text)
    if match_k:
        result["Kstar"] = float(match_k.group(1))
    # epsilon
    match_eps = re.search(r"最終epsilon *= *([\d\.Ee+-]+)", text)
    if match_eps:
        result["epsilon"] = float(match_eps.group(1))
    # Lf（燃料長さ）
    match_lf = re.search(r"燃料長さ *= *([\d\.Ee+-]+)", text)
    if match_lf:
        result["Lf"] = float(match_lf.group(1))
    # mdot
    match_mdot = re.search(r"最終mdot *= *([\d\.Ee+-]+)", text)
    if match_mdot:
        result["mdot"] = float(match_mdot.group(1))
    # 初期推力F
    match_F = re.search(r"最終推力 *= *([\d\.Ee+-]+)", text)
    if match_F:
        result["F"] = float(match_F.group(1))
    # Dt
    match_Dt = re.search(r"Dt *= *([\d\.Ee+-]+)", text)
    if match_Dt:
        result["Dt"] = float(match_Dt.group(1))

    return result

# Ptからrho_oxを計算
def calc_rho_ox(pressure):
    ox_calc_result = ox_db.get_density(pressure)
    rho_ox = float(ox_calc_result.split(":")[-1].replace("kg/m³", "").strip())
    return rho_ox

# -------------------------
# 初期条件計算モード
# -------------------------
def run_initial_condition_mode():
    print("initial condition mode selected.")
    print("input jsonc filename(example.jsonc):")
    filename = input("> ").strip()

    try:
        loader = JsoncLoader(filename)
        inputvalues = loader.load()
    except Exception as e:
        print(f"loading error: {e}")
        exit(1)

    print("loading values:")
    for key, val in inputvalues.items():
        print(f"{key} = {val}")

    # Ptからrho_oxを計算
    rho_ox = calc_rho_ox(inputvalues["Pt_init"])

    # initial_convergence()に投げる部分
    F_req = inputvalues["F_req"]
    Pc_def = inputvalues["Pc_def"]
    OF_def = inputvalues["OF_def"]
    mdot_new = inputvalues["mdot_new"]
    Df_init = inputvalues["Df_init"]
    eta_cstar = inputvalues["eta_cstar"]
    eta_nozzle = inputvalues["eta_nozzle"]
    Ptank_init = inputvalues["Pt_init"]
    rho_ox_init = rho_ox
    rho_f_start = inputvalues["rho_f"]
    a_ox = inputvalues["a_ox"]
    n_ox = inputvalues["n_ox"]

    init_output, _, _ = sim.initial_convergence(
        F_req, Pc_def, OF_def, mdot_new, Df_init,
        eta_cstar, eta_nozzle, Ptank_init,
        rho_ox_init, rho_f_start, a_ox, n_ox
    )

# -------------------------
# 時間発展計算モード（これから作成）
# -------------------------
def run_time_evolution_mode():
    print("time_evolution mode selected")
    print("input jsonc filename(example.jsonc):")
    filename = input("> ").strip()

    try:
        loader = JsoncLoader(filename)
        inputvalues = loader.load()
    except Exception as e:
        print(f"loading error: {e}")
        exit(1)

    print("loading values:")
    for key, val in inputvalues.items():
        print(f"{key} = {val}")

    # Ptからrho_oxを計算
    rho_ox = calc_rho_ox(inputvalues["Pt_init"])

    # integration_simulation()に投げる部分
    F = inputvalues["F_init"]
    Pc = inputvalues["Pc_def"]
    OF = inputvalues["OF_def"]
    mdot = inputvalues["mdot_new"]
    Df = inputvalues["Df_init"]
    eta_cstar = inputvalues["eta_cstar"]
    eta_nozzle = inputvalues["eta_nozzle"]
    P_init = inputvalues["Pt_init"]
    rho_f = inputvalues["rho_f"]
    a_ox = inputvalues["a_ox"]
    n_ox = inputvalues["n_ox"]
    Kstar = inputvalues["Kstar"]
    epsilon = inputvalues["epsilon"]
    Lf = inputvalues["Lf"]
    V_tank = inputvalues["Vol_ox"]
    P_final = inputvalues["Pt_end"]
    Dt = inputvalues["Dt"]

    (_, _, _, _, _, _, _, evolution_result, It,) = sim.integration_simulation(
        Pc=Pc, Df=Df, OF=OF, eta_cstar=eta_cstar, eta_nozzle=eta_nozzle, Kstar=Kstar,
        epsilon=epsilon, Lf=Lf, mdot=mdot, V_tank=V_tank, P_init=P_init, P_final=P_final,
        rho_ox=rho_ox, rho_fuel=rho_f, a=a_ox, n=n_ox, F=F, Dt=Dt)
    
    # 結果出力
    print("input output csv filename(example.csv):")
    output_filename = input("> ").strip()

    try:
        # input記載用
        input_params = [
                ("Pc", Pc), ("Df", Df), ("OF", OF),
                ("eta_cstar", eta_cstar), ("eta_nozzle", eta_nozzle), ("Kstar", Kstar),
                ("epsilon", epsilon), ("Lf", Lf), ("mdot", mdot),
                ("V_tank", V_tank), ("P_init", P_init), ("P_final", P_final),
                ("rho_ox", rho_ox), ("rho_fuel", rho_f),
                ("a", a_ox), ("n", n_ox), ("F", F), ("Dt", Dt)
            ]
        
        # 時間発展記載用
        evolution_headers = [
                "F [N]",
                "F_fte [N]",
                "Ptank [MPa]",
                "Pc [MPa]",
                "O/F [-]",
                "mdot [kg/s]",
                "Df [m]",
                "C* [m/s]",
                "CF [-]",
                "tank mass [g]",
                "mdot_ox [g/ms]",
                "gamma [-]"
            ]

        # データ本体出力
        filename = output_filename
        with open(filename, "w", newline="", encoding="utf-8") as file:
            writer = csv.writer(file, quoting=csv.QUOTE_NONE)
            # 入力パラメータの書き出し
            writer.writerow(["# input params."])
            for i in range(0, len(input_params), 3):
                row = []
                for j in range(3):
                    if i + j < len(input_params):
                        key, val = input_params[i + j]
                        row.extend([key, val])
                writer.writerow(row)


            writer.writerow([])  # 空行
            writer.writerow(["# evolution params."])
            writer.writerow(evolution_headers)
            writer.writerows(evolution_result)
    except Exception as e:
        print(f"loading error: {e}")
        exit(1)

# -------------------------
# メイン処理
# -------------------------
if __name__ == "__main__":
    print("Select mode:")
    print("1: initial condition mode")
    print("2: time_evolution mode")
    mode = input("> ").strip()

    if mode == "1":
        run_initial_condition_mode()
    elif mode == "2":
        run_time_evolution_mode()
    else:
        print("Invalid mode selected.")