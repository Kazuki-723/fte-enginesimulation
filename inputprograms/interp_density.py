import pandas as pd
from scipy.interpolate import interp1d

class OxidizerDatabase:
    def __init__(self):
        df_liq = pd.read_csv("inputdatas\\N2O_liquid.csv")
        self.li_interp = interp1d(df_liq['Pres'].values, df_liq['Dens'].values, kind='cubic', fill_value='extrapolate')
        self.liq_range = (min(df_liq['Pres'].values), max(df_liq['Pres'].values))

        df_vap = pd.read_csv("inputdatas\\N2O_vapor.csv")
        self.va_interp = interp1d(df_vap['Pres'].values, df_vap['Dens'].values, kind='cubic', fill_value='extrapolate')
        self.vap_range = (min(df_vap['Pres'].values), max(df_vap['Pres'].values))

    def get_density(self, pressure: float) -> str:
        if self.liq_range[0] <= pressure <= self.liq_range[1]:
            rho = self.li_interp(pressure)
            phase = "液相"
        elif self.vap_range[0] <= pressure <= self.vap_range[1]:
            rho = self.va_interp(pressure)
            phase = "気相"
        else:
            return "圧力範囲外"
        return f"{phase}密度: {rho:.2f} kg/m³"