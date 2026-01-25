import re
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
    ox_calc_result = ox_db.get_density(inputvalues["Pt_init"])
    rho_ox = float(ox_calc_result.split(":")[-1].replace("kg/m³", "").strip())

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

    # 結果のパース
    init_parse_result = parse_initial_results(init_output)
    print(init_parse_result)


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

    # ここに時間発展コードを追加していく
    print("時間発展計算を開始します（まだ未実装）")
    # TODO: ここに時間発展ロジックを実装


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