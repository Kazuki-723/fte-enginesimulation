import pandas as pd
from scipy.interpolate import interp1d

class OxidizerDatabase:
    """
    NIST databaseより液相気相それぞれの気液平衡時における密度を圧力から補完するクラス
    databaseの作り次第だが，作成者は0.1 - 7MPaの範囲で補完
    圧力範囲外のエラーは出さない(原理上外挿するはず)が値の信頼性は担保されてない．そもそも7超えたら臨界圧すぐだし
    """

    def __init__(self):
        # liquid phase data loading
        df_liq = pd.read_csv("inputdatas\\N2O_liquid.csv")
        self.li_interp = interp1d(df_liq['Pres'].values, df_liq['Dens'].values, kind='cubic', fill_value='extrapolate')
        self.liq_range = (min(df_liq['Pres'].values), max(df_liq['Pres'].values))

        # gas phase data loading
        df_vap = pd.read_csv("inputdatas\\N2O_vapor.csv")
        self.va_interp = interp1d(df_vap['Pres'].values, df_vap['Dens'].values, kind='cubic', fill_value='extrapolate')
        self.vap_range = (min(df_vap['Pres'].values), max(df_vap['Pres'].values))

    # csvから線形補完する部分
    def get_density(self, pressure: float, phase ="liquid") -> str:
        if phase =="liquid":
            rho = self.li_interp(pressure)
        elif phase == "gas":
            rho = self.va_interp(pressure)
        else:
            return "error" # 相入力判定失敗時
        return phase, rho