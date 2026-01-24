from inputprograms.rocket_simulation import RocketSimulation
from inputprograms.interp_density import OxidizerDatabase
from inputprograms.importjson import JsoncLoader


if __name__ == "__main__":
    print("input jsonc filename(examplr.jsonc):")
    filename = input("> ").strip()

    try:
        loader = JsoncLoader(filename)
        values = loader.load()
    except Exception as e:
        print(f"loading error: {e}")
        exit(1)

    # values に JSONC の内容が dict として入っている
    print("loading values:")
    for key, val in values.items():
        print(f"{key} = {val}")

    # Ptからrho_oxを計算
    ox_db = OxidizerDatabase()
    ox_calc_result = ox_db.get_density(values["Pt_init"])
    rho_ox = float(ox_calc_result.split(":")[-1].replace("kg/m³", "").strip())
    print(rho_ox)
