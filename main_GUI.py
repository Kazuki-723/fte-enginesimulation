import flet as ft
from inputprograms.rocket_simulation import RocketSimulation

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

    # 時間発展ビュー（別ページ）
    def evolution_view():
        initial = page.session.get("initial_conditions")
        if not initial:
            return ft.View(route="/evolution", controls=[
                ft.Text("⚠️ 初期条件が未設定です"),
                ft.TextButton("◀ 戻る", on_click=lambda _: page.go("/"))
            ])

        # 物質と密度の辞書
        materials = {
            "液体酸素 (LOX)": 1141,
            "液体水素 (LH2)": 71,
            "RP-1 (ケロシン)": 810,
            "メタン (CH4)": 422,
            "N2O4 (四酸化二窒素)": 1440,
            "UDMH (ジメチルヒドラジン)": 791
        }

        # 密度表示用テキスト
        selected_density = ft.Text(value="密度: -", size=16)

        # プルダウン選択イベント
        def on_material_change(e):
            name = e.control.value
            rho = materials.get(name, "-")
            selected_density.value = f"密度: {rho} kg/m³" if rho != "-" else "密度: -"
            page.update()

        # プルダウンコンポーネント
        material_dropdown = ft.Dropdown(
            label="推進剤を選択",
            options=[ft.dropdown.Option(name) for name in materials.keys()],
            on_change=on_material_change,
            width=250,
            value="液体酸素 (LOX)"
        )

        # 仮の時間発展出力
        evolution_output = ft.Text("🕒 時間発展シミュレーション（仮表示）")

        return ft.View(
            route="/evolution",
            controls=[
                ft.Text("時間発展ページ", size=20, weight=ft.FontWeight.BOLD),
                ft.Row(
                    controls=[material_dropdown, selected_density],
                    alignment=ft.MainAxisAlignment.START
                ),
                evolution_output,
                ft.TextButton("◀ 戻る", on_click=lambda _: page.go("/"))
            ]
        )
ft.app(target=main)