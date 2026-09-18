import csv
from inputprograms.rocket_simulation import RocketSimulation
from inputprograms.importjson import JsoncLoader
sim = RocketSimulation()

# 入力のjsonエラー判定関数
def validate_inputs(required_keys: dict, inputvalues: dict):
    """
    required_keys: {"key_name": expected_type}
    inputvalues: JSONC から読み込んだ dict
    """

    errors = []

    # -------------------------
    # 1. 欠損キーのチェック
    #--------------------------
    missing_keys = [key for key in required_keys if key not in inputvalues]
    if missing_keys:
        errors.append("Missing required keys:")
        for key in missing_keys:
            errors.append(f"  - {key}")

    # -------------------------
    # 2. 型チェック
    #--------------------------
    type_errors = []
    for key, expected_type in required_keys.items():
        if key in inputvalues:
            val = inputvalues[key]
            if expected_type is float:
                if not isinstance(val, (int, float)):
                    type_errors.append(f"{key}: expected float, got {type(val).__name__}")
            elif expected_type is str:
                if not isinstance(val, str):
                    type_errors.append(f"{key}: expected str, got {type(val).__name__}")
            elif expected_type is int:
                if not isinstance(val, int):
                    type_errors.append(f"{key}: expected int, got {type(val).__name__}")
            else:
                type_errors.append(f"{key}: unknown expected type {expected_type}")

    if type_errors:
        errors.append("Invalid types:")
        for err in type_errors:
            errors.append(f"  - {err}")

    # -------------------------
    # 3. 不要キーのチェック
    #--------------------------
    extra_keys = [key for key in inputvalues if key not in required_keys]
    if extra_keys:
        errors.append("Extra keys found (not used by simulation):")
        for key in extra_keys:
            errors.append(f"  - {key}")

    # -------------------------
    # 4. エラーがあればまとめて出力して終了
    #--------------------------
    if errors:
        print("ERROR: Invalid JSON input detected:")
        for e in errors:
            print(e)
        exit(1)



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

    # 入力エラー判定用正解の提示
    required_keys_init = {
    "F_req": float,
    "Pc_def": float,
    "OF_def": float,
    "mdot_new": float,
    "Df_init": float,
    "eta_cstar": float,
    "eta_nozzle": float,
    "Pt_init": float,
    "rho_f": float,
    "a_ox": float,
    "n_ox": float,
    "fuel_material": str,
    }

    # キー欠損，余剰，型チェック
    validate_inputs(required_keys_init, inputvalues)

    # 入力チェック後にinitial_convergence()に投げる用
    F_req         = inputvalues["F_req"]
    Pc_def        = inputvalues["Pc_def"]
    OF_def        = inputvalues["OF_def"]
    mdot_new      = inputvalues["mdot_new"]
    Df_init       = inputvalues["Df_init"]
    eta_cstar     = inputvalues["eta_cstar"]
    eta_nozzle    = inputvalues["eta_nozzle"]
    Ptank_init    = inputvalues["Pt_init"]
    rho_f_start   = inputvalues["rho_f"]
    a_ox          = inputvalues["a_ox"]
    n_ox          = inputvalues["n_ox"]
    fuel_material = inputvalues["fuel_material"]

    # Ptからrho_oxを計算
    _, rho_ox = sim.calc_rho_ox(Ptank_init, "liquid")

    # 入力チェック後にinput提示
    print("loading values:")
    for key, val in inputvalues.items():
        print(f"{key} = {val}")

    _ = sim.initial_convergence(
        F_req, Pc_def, OF_def, mdot_new, Df_init,
        eta_cstar, eta_nozzle, Ptank_init,
        rho_ox, rho_f_start, a_ox, n_ox, fuel_material
    )

# -------------------------
# 時間発展計算モード
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

    # 入力エラー判定用正解の提示
    required_keys_time = {
        "F_init": float,
        "Pc_def": float,
        "OF_def": float,
        "mdot_new": float,
        "Df_init": float,
        "eta_cstar": float,
        "eta_nozzle": float,
        "Pt_init": float,
        "rho_f": float,
        "a_ox": float,
        "n_ox": float,
        "fuel_material": str,
        "Kstar": float,
        "epsilon": float,
        "Lf": float,
        "Vol_ox": float,
        "Pt_end": float,
        "Dt": float,
        "is_fast": int,
    }

    # キー欠損，余剰，型チェック
    validate_inputs(required_keys_time, inputvalues)

    # integration_simulation()に投げる部分
    F             = inputvalues["F_init"]
    Pc            = inputvalues["Pc_def"]
    OF            = inputvalues["OF_def"]
    mdot          = inputvalues["mdot_new"]
    Df            = inputvalues["Df_init"]
    eta_cstar     = inputvalues["eta_cstar"]
    eta_nozzle    = inputvalues["eta_nozzle"]
    Ptank_init    = inputvalues["Pt_init"]
    rho_f         = inputvalues["rho_f"]
    a_ox          = inputvalues["a_ox"]
    n_ox          = inputvalues["n_ox"]
    fuel_material = inputvalues["fuel_material"]
    Kstar         = inputvalues["Kstar"]
    epsilon       = inputvalues["epsilon"]
    Lf            = inputvalues["Lf"]
    V_tank        = inputvalues["Vol_ox"]
    P_final       = inputvalues["Pt_end"]
    Dt            = inputvalues["Dt"]
    is_fast       = inputvalues["is_fast"]

    # Ptからrho_oxを計算
    _, rho_ox = sim.calc_rho_ox(Ptank_init, "liquid")

    print("loading values:")
    for key, val in inputvalues.items():
        print(f"{key} = {val}")

    if is_fast == 1:
        print("normal mode")
    elif is_fast > 1:
        print("fast mode")
    else:
        print("invailed settings, set to normal mode")
        is_fast = 1

    cea_interval = is_fast
    # normal
    (_, _, _, _, _, _, _, evolution_result, _,) = sim.integration_simulation(
        Pc=Pc, Df=Df, OF=OF, eta_cstar=eta_cstar, eta_nozzle=eta_nozzle, Kstar=Kstar,
        epsilon=epsilon, Lf=Lf, mdot=mdot, V_tank=V_tank, P_init=Ptank_init, P_final=P_final,
        rho_ox=rho_ox, rho_fuel=rho_f, a=a_ox, n=n_ox, fuel_material = fuel_material,
        F=F, Dt=Dt, cea_interval = cea_interval)
    
    # 結果出力
    print("input output csv filename(example.csv):")
    output_filename = input("> ").strip()

    try:
        # input記載用
        input_params = [
                ("Pc", Pc), ("Df", Df), ("OF", OF),
                ("eta_cstar", eta_cstar), ("eta_nozzle", eta_nozzle), ("Kstar", Kstar),
                ("epsilon", epsilon), ("Lf", Lf), ("mdot", mdot),
                ("V_tank", V_tank), ("P_init", Ptank_init), ("P_final", P_final),
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