import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# liquid csv load
df = pd.read_csv("inputdatas\\N2O_liquid.csv")
li_pressure = df['Pres'].values
li_density = df['Dens'].values
li_temperture = df['Temp'].values

# 補完関数の作成（線形補完またはスプライン補完）
li_interp = interp1d(li_pressure, li_density, kind='cubic', fill_value='extrapolate')

# 補完用の圧力範囲を生成
li_pressure_new = np.linspace(min(li_pressure), max(li_pressure), 500)
li_density_new = li_interp(li_pressure_new)

# vapor csv load
df = pd.read_csv("inputdatas\\N2O_vapor.csv")
va_pressure = df['Pres'].values
va_density = df['Dens'].values
va_temperture = df['Temp'].values

# 補完関数の作成（線形補完またはスプライン補完）
va_interp = interp1d(va_pressure, va_density, kind='cubic', fill_value='extrapolate')

# 補完用の圧力範囲を生成
va_pressure_new = np.linspace(min(va_pressure), max(va_pressure), 500)
va_density_new = va_interp(va_pressure_new)

# グラフの描画
plt.figure(figsize=(8, 5))
plt.plot(li_pressure, li_density, 'o', label='liquid Original Data')
plt.plot(li_pressure_new, li_density_new, '-', label=';liquid Interpolated Curve')
plt.plot(va_pressure, va_density, 'o', label='vapor Original Data')
plt.plot(va_pressure_new, va_density_new, '-', label=';vapor Interpolated Curve')
plt.xlabel('Pressure')
plt.ylabel('Density')
plt.title('Pressure vs Density')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()