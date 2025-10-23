import flet as ft
import pandas as pd
from scipy.interpolate import interp1d

# CSV読み込みと補間関数の準備
df = pd.read_csv("inputdatas\\N2O.csv")  # Temp, Pres, Dens の列がある前提

# 圧力に対する温度・密度の補間関数を構築
temp_interp = interp1d(df["Pres"], df["Temp"], kind="linear", fill_value="extrapolate")
dens_interp = interp1d(df["Pres"], df["Dens"], kind="linear", fill_value="extrapolate")

def main(page: ft.Page):
    page.title = "気液平衡 補間検索"
    page.scroll = ft.ScrollMode.AUTO

    pres_input = ft.TextField(label="圧力 [MPa]", width=200)
    result_text = ft.Text()

    def search_data(e):
        try:
            input_pres = float(pres_input.value)
        except ValueError:
            result_text.value = "⚠️ 数値で圧力を入力してください"
            page.update()
            return

        # 補間による温度・密度の推定
        temp = float(temp_interp(input_pres))
        dens = float(dens_interp(input_pres))

        result_text.value = f"🔍 入力圧力: {input_pres:.2f} MPa\n🌡️ 推定温度: {temp:.2f} K\n🧪 推定密度: {dens:.4f} kg/m³"
        page.update()

    page.add(
        ft.Row([pres_input, ft.ElevatedButton("補間検索", on_click=search_data)]),
        result_text
    )

ft.app(target=main)