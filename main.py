from inputprograms.rocket_simulation import RocketSimulation
from inputprograms.interp_density import OxidizerDatabase
from inputprograms.importjson import JsoncLoader
sim = RocketSimulation()
ox_db = OxidizerDatabase()


if __name__ == "__main__":
    print("input jsonc filename(examplr.jsonc):")
    filename = input("> ").strip()

    try:
        loader = JsoncLoader(filename)
        inputvalues = loader.load()
    except Exception as e:
        print(f"loading error: {e}")
        exit(1)

    # inputvalues に JSONC の内容が dict として入っている
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
    output, _, _ = sim.initial_convergence(F_req, Pc_def, OF_def, mdot_new, Df_init, eta_cstar, eta_nozzle, Ptank_init, rho_ox_init, rho_f_start, a_ox, n_ox)