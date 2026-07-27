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

    def get_density(self, pressure: float, phase ="liquid") -> str:
        if phase =="liquid":
            rho = self.li_interp(pressure)
        elif phase == "gas":
            rho = self.va_interp(pressure)
        else:
            return "error"
        return phase, rho