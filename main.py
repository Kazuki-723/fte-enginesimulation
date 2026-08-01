import csv
import time 
from inputprograms.rocket_simulation import RocketSimulation
from inputprograms.importjson import JsoncLoader
sim = RocketSimulation()

# import tracemalloc
# tracemalloc.start(15)

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

    # 実行時間計測(初期設定に依存するのでsampleに対する参考値)
    start = time.perf_counter()

    print("loading values:")
    for key, val in inputvalues.items():
        print(f"{key} = {val}")

    # Ptからrho_oxを計算
    _, rho_ox = sim.calc_rho_ox(inputvalues["Pt_init"], "liquid")

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

    _, _, _ = sim.initial_convergence(
        F_req, Pc_def, OF_def, mdot_new, Df_init,
        eta_cstar, eta_nozzle, Ptank_init,
        rho_ox_init, rho_f_start, a_ox, n_ox
    )

    end = time.perf_counter()
    print(f"Elapsed time: {end - start:.6f} seconds")

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

    # 実行時間計測
    start = time.perf_counter()

    print("loading values:")
    for key, val in inputvalues.items():
        print(f"{key} = {val}")

    # Ptからrho_oxを計算
    _, rho_ox = sim.calc_rho_ox(inputvalues["Pt_init"], "liquid")

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

    # normal
    (_, _, _, _, _, _, _, evolution_result, _,) = sim.integration_simulation(
        Pc=Pc, Df=Df, OF=OF, eta_cstar=eta_cstar, eta_nozzle=eta_nozzle, Kstar=Kstar,
        epsilon=epsilon, Lf=Lf, mdot=mdot, V_tank=V_tank, P_init=P_init, P_final=P_final,
        rho_ox=rho_ox, rho_fuel=rho_f, a=a_ox, n=n_ox, F=F, Dt=Dt)

    # 計算本体の時間をはかる
    end = time.perf_counter()
    print(f"Elapsed time: {end - start:.6f} seconds")

    # snapshot = tracemalloc.take_snapshot()
    # top_stats = snapshot.statistics('traceback')

    # print("[ Top 10 ]")
    # for stat in top_stats[:10]:
    #     print(stat)
    #     for line in stat.traceback.format():
    #         print(line)
    #     print("=====")

        
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