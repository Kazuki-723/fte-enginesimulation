import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 読み込むファイルの名前を入れる
filename = "2000ccV2.csv"
df = pd.read_csv(filename, skiprows=9, header=0)

columns = df.columns
data_arrays = {col: df[col].to_numpy() for col in columns}

# arrays
# F [N],F_fte [N],Ptank [MPa],Pc [MPa],O/F [-],mdot [kg/s],Df [m],C* [m/s],CF [-],tank mass [g],mdot_ox [g/ms],gamma [-]
F_arr = data_arrays["F [N]"]
F_fte_arr = data_arrays["F_fte [N]"]
Pt_arr = data_arrays["Ptank [MPa]"]
Pc_arr  = data_arrays["Pc [MPa]"]
OF_arr = data_arrays["O/F [-]"]
mdot_arr = data_arrays["mdot [kg/s]"]
Df_arr = data_arrays["Df [m]"]
Cstar_arr = data_arrays["C* [m/s]"]
CF_arr  = data_arrays["CF [-]"]
tank_mass_arr = data_arrays["tank mass [g]"]
mdot_ox_arr = data_arrays["mdot_ox [g/ms]"]
gamma_arr = data_arrays["gamma [-]"]
mdot_f_arr = mdot_arr - mdot_ox_arr

# make timearray
time = np.arange(0, len(F_arr) / 1000, 0.001)

# plot

#--
# 1. plot two methods of thrust
#--
plt.figure(figsize=(10, 5))
plt.plot(time, F_arr, label="thrust1")
plt.plot(time, F_fte_arr, label="thrust2", linestyle="--")
plt.xlabel("Time [s]")
plt.ylabel("Thrust [N]")
plt.title("compare two thrust methods")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#--
# 2.Cstar vs time
#--

plt.figure(figsize=(10, 5))
plt.plot(time, Cstar_arr, label="Cstar")
plt.xlabel("Time [s]")
plt.ylabel("Cstar [m/s]")
plt.title("Cstar vs Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#--
# 3. Pressure vs time
#--
plt.figure(figsize=(10, 5))
plt.plot(time, Pt_arr, label="tank")
plt.plot(time, Pc_arr, label="chamnber", linestyle="--")
plt.xlabel("Time [s]")
plt.ylabel("Pressure [MPa]")
plt.title("Pressure vs Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#--
# 4. mdot vs time(1)
#--
plt.figure(figsize=(10, 5))
plt.plot(time, mdot_ox_arr, label="ox")
plt.plot(time, mdot_arr, label="total", linestyle="--")
plt.xlabel("Time [s]")
plt.ylabel("Mass flow Rate [kg/s]")
plt.title("Mass Flow Rate vs Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#--
# 5. fuel flow rate vs time
#--
plt.figure(figsize=(10, 5))
plt.plot(time, mdot_f_arr, label="fuel")
plt.xlabel("Time [s]")
plt.ylabel("Mass flow Rate [kg/s]")
plt.title("Mass Flow Rate vs Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()