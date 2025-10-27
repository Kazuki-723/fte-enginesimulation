import flet as ft
from inputprograms.rocket_simulation import RocketSimulation
from inputprograms.interp_density import OxidizerDatabase
import re


def main(page: ft.Page):
    page.title = "Rocket Simulation GUI"
    page.scroll = ft.ScrollMode.AUTO

    # ページ切り替え処理
    def route_change(route):
        page.views.clear()
        if page.route == "/":
            page.views.append(main_view())
        elif page.route == "/evolution":
            page.views.append(evolution_view())
        page.update()

    page.on_route_change = route_change
    page.go("/")

    # メインビュー（初期条件＋収束）
    def main_view():
        inputs = {
            "F_req": ft.TextField(label="要求推力 [N]", width=150, value=650),
            "Pc_def": ft.TextField(label="初期燃焼室圧力 [MPa]", width=150, value=2),
            "OF_def": ft.TextField(label="初期O/F比", width=150, value=6.5),
            "mdot_new": ft.TextField(label="初期流量 [kg/s]", width=150, value=0.33),
            "Df_init":ft.TextField(label="初期燃料内径 [m]", width=150, value=0.034),
            "eta_cstar": ft.TextField(label="C*効率", width=150, value=0.8),
            "eta_nozzle": ft.TextField(label="ノズル効率", width=150,value=0.98),
        }

        result_text = ft.Text()
        graph_image = ft.Image(visible=False, width=page.window_width - 200)

        def run_simulation(e):
            try:
                values = {k: float(inputs[k].value) for k in inputs}
            except ValueError:
                result_text.value = "⚠️ 全ての値を数値で入力してください"
                page.update()
                return

            sim = RocketSimulation()
            output, Dovalue, cdvalue = sim.initial_convergence(**values)
            result_text.value = output

            graph_image.src_base64 = sim.get_iteration_plot_base64(Dovalue, cdvalue)
            graph_image.visible = True

            page.session.set("initial_conditions", values)  # 初期条件保存
            page.session.set("initial_results", output)  # 出力保存
            page.update()
        
        # 実行ボタンと遷移ボタンを並べる
        action_row = ft.Row([
            ft.ElevatedButton("収束計算", on_click=run_simulation),
            ft.TextButton("▶ 時間発展ページへ", on_click=lambda _: page.go("/evolution"))
        ])

        # 左側：入力群＋結果＋ボタン群＋ログ
        input_column = ft.Column(
            controls=[
                *list(inputs.values()),
                action_row,
                result_text
            ],
            spacing=10,
            expand=True,
            height=page.window_height + 100,
            scroll=ft.ScrollMode.AUTO
        )

        # 右側：収束グラフと K* グラフを縦に並べる
        graph_image = ft.Image(visible=False, expand=True)

        graph_column = ft.Column(
            controls=[graph_image],
            spacing=10,
            alignment=ft.MainAxisAlignment.START
        )

        return ft.View(
            route="/",
            controls=[
                ft.Row(
                    controls=[input_column, graph_column],
                    alignment=ft.MainAxisAlignment.START,
                    vertical_alignment=ft.CrossAxisAlignment.START
                )
            ]
        )
    # ちょっとパース
    def parse_initial_results(text: str) -> dict:
        result = {}

        # K*
        match_k = re.search(r"K\* *= *([\d\.Ee+-]+)", text)
        if match_k:
            result["Kstar"] = float(match_k.group(1))

        # epsilon
        match_eps = re.search(r"最終epsilon *= *([\d\.Ee+-]+)", text)
        if match_eps:
            result["epsilon"] = float(match_eps.group(1))

        # Lf（燃料長さ）
        match_lf = re.search(r"燃料長さ *= *([\d\.Ee+-]+)", text)
        if match_lf:
            result["Lf"] = float(match_lf.group(1))

        return result


    # 時間発展ビュー（別ページ）
    def evolution_view():
        initial = page.session.get("initial_conditions")
        results = page.session.get("initial_results")  # K*, epsilon, Lf を含む
        if results != None:
            results_parsed = parse_initial_results(results)

        print(results)
        # 初期値がある場合は値を埋める、なければ空欄
        Pc_def = str(initial["Pc_def"]) if initial else ""
        Df_init = str(initial["Df_init"]) if initial else ""
        eta_cstar = str(initial["eta_cstar"]) if initial else ""
        eta_nozzle = str(initial["eta_nozzle"]) if initial else ""

        Kstar = str(results_parsed["Kstar"]) if initial else ""
        epsilon = str(results_parsed["epsilon"]) if initial else ""
        Lf = str(results_parsed["Lf"]) if initial else ""

        # 入力欄の定義
        Pc_box = ft.TextField(label="燃焼室圧力 Pc [MPa]", value=Pc_def, width=150)
        Df_box = ft.TextField(label="初期ポート径 Df [m]", value=Df_init, width=150)
        eta_cstar_box = ft.TextField(label="C*効率", value=eta_cstar, width=150)
        eta_nozzle_box = ft.TextField(label="ノズル効率", value=eta_nozzle, width=150)

        Kstar_box = ft.TextField(label="K*", value=Kstar, width=150)
        epsilon_box = ft.TextField(label="膨張比 ε", value=epsilon, width=150)
        Lf_box = ft.TextField(label="燃焼長 Lf [m]", value=Lf, width=150)

        # 登録物質と物性値（a, n は仮値）
        materials = {
            "PMMA": {"密度": 1180, "a": 0.85, "n": 1.2},
            "ABS":  {"密度": 1040, "a": 0.90, "n": 1.1}
        }

        # 表示用テキスト群
        density_text = ft.Text(value="密度: -", size=16)
        a_text = ft.Text(value="a: -", size=16)
        n_text = ft.Text(value="n: -", size=16)

        def on_material_change(e):
            name = e.control.value
            props = materials.get(name, {})
            density_text.value = f"密度: {props.get('密度', '-')} kg/m³"
            a_text.value = f"a: {props.get('a', '-')}"
            n_text.value = f"n: {props.get('n', '-')}"
            page.session.set("material_properties", props)  # RocketSimulation側に渡す準備
            page.update()

        material_dropdown = ft.Dropdown(
            label="固体燃料を選択",
            options=[ft.dropdown.Option(name) for name in materials.keys()],
            on_change=on_material_change,
            width=250
        )

        property_column = ft.Column(
            controls=[density_text, a_text, n_text],
            spacing=5
        )

        # 酸化剤補完データベース
        ox_db = OxidizerDatabase()

        pressure_input = ft.TextField(label="酸化剤圧力 [MPa]", width=150)
        density_output = ft.Text(value="酸化剤密度: -", size=16)

        def on_pressure_change(e):
            try:
                p = float(pressure_input.value)
                result = ox_db.get_density(p)
                density_output.value = result
            except ValueError:
                density_output.value = "⚠️ 数値で入力してください"
            page.update()

        pressure_input.on_change = on_pressure_change

        evolution_output = ft.Text("🕒 時間発展シミュレーション（仮表示）")

        return ft.View(
            route="/evolution",
            controls=[
                ft.Text("時間発展ページ", size=20, weight=ft.FontWeight.BOLD),
                ft.Row(
                    controls=[
                        ft.Column([
                            ft.Text("初期状態パラメータ："),
                            Pc_box, Df_box, eta_cstar_box, eta_nozzle_box,
                            Kstar_box, epsilon_box, Lf_box,
                            material_dropdown,
                            property_column,
                            pressure_input,
                            density_output
                        ],)
                    ],
                    alignment=ft.MainAxisAlignment.START
                ),
                evolution_output,
                ft.TextButton("◀ 戻る", on_click=lambda _: page.go("/"))
            ]
        )


ft.app(target=main)