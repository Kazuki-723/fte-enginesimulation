import flet as ft
import csv
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
            "Df_init": ft.TextField(label="初期燃料内径 [m]", width=150, value=0.034),
            "eta_cstar": ft.TextField(label="C*効率", width=150, value=0.8),
            "eta_nozzle": ft.TextField(label="ノズル効率", width=150, value=0.98),
        }

        result_text = ft.Text()
        graph_image = ft.Image(visible=False, width=page.window.width - 200)

                # 登録物質と物性値（a, n は仮値）
        materials = {
            "PMMA": {"密度": 1190, "a": 0.000131, "n": 0.34},
            "ABS": {"密度": 1040, "a": 0.90, "n": 1.1},
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
            page.session.set(
                "material_properties", props
            )  # RocketSimulation側に渡す準備
            page.update()

        material_dropdown = ft.Dropdown(
            label="固体燃料を選択",
            options=[ft.dropdown.Option(name) for name in materials.keys()],
            on_change=on_material_change,
            width=250,
        )

        property_column = ft.Column(controls=[density_text, a_text, n_text], spacing=5)

        # 酸化剤補完データベース
        ox_db = OxidizerDatabase()

        pressure_input = ft.TextField(label="初期酸化剤圧力 [MPa]", width=150)
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

        def run_simulation(e):
            try:
                values = {k: float(inputs[k].value) for k in inputs}
                print(values)
                # 追加項目の取得と格納
                # 酸化剤密度（補完済みテキストから抽出）
                rho_ox = float(
                    density_output.value.split(":")[-1].replace("kg/m³", "").strip()
                )
                # 燃料密度・定数a,n（Dropdown選択から取得）
                material_props = page.session.get("material_properties")
                rho_fuel = float(material_props["密度"])
                a = float(material_props["a"])
                n = float(material_props["n"])

                values["Ptank_init"] = float(pressure_input.value)
                values["rho_ox_init"] = rho_ox
                values["rho_f_start"] = rho_fuel

                material_props = page.session.get("material_properties")
                values["a_ox"] = a
                values["n_ox"] = n
                print(values)

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
        action_row = ft.Row(
            [
                ft.ElevatedButton("収束計算", on_click=run_simulation),
                ft.TextButton(
                    "▶ 時間発展ページへ", on_click=lambda _: page.go("/evolution")
                ),
            ]
        )

        # 左側：入力群＋結果＋ボタン群＋ログ
        input_column = ft.Column(
            controls=[*list(inputs.values()),
                    pressure_input,
                    density_output,
                    material_dropdown,
                    property_column,
                    action_row, 
                    result_text],
            spacing=10,
            expand=True,
            height=page.window.height + 100,
            scroll=ft.ScrollMode.AUTO,
        )

        # 右側：収束グラフと K* グラフを縦に並べる
        graph_image = ft.Image(visible=False)

        graph_column = ft.Column(
            controls=[graph_image],
            spacing=10,
            expand=True,
            height=page.window.height + 100,
            scroll=ft.ScrollMode.AUTO,
            alignment=ft.MainAxisAlignment.START,
        )

        return ft.View(
            route="/",
            controls=[
                ft.Row(
                    controls=[input_column, graph_column],
                    alignment=ft.MainAxisAlignment.START,
                    vertical_alignment=ft.CrossAxisAlignment.START,
                )
            ],
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
        # mdot
        match_mdot = re.search(r"最終mdot *= *([\d\.Ee+-]+)", text)
        if match_mdot:
            result["mdot"] = float(match_mdot.group(1))
        # 初期推力F
        match_F = re.search(r"最終推力 *= *([\d\.Ee+-]+)", text)
        if match_F:
            result["F"] = float(match_F.group(1))
        # Dt
        match_Dt = re.search(r"Dt *= *([\d\.Ee+-]+)", text)
        if match_Dt:
            result["Dt"] = float(match_Dt.group(1))

        return result

    # 時間発展ビュー（別ページ）
    def evolution_view():
        initial = page.session.get("initial_conditions")
        results = page.session.get("initial_results")  # K*, epsilon, Lf を含む
        if results != None:
            results_parsed = parse_initial_results(results)

        results_graph_image = ft.Image(visible=False, width=600)

        # 初期値がある場合は値を埋める、なければ空欄
        Pc_def = str(initial["Pc_def"]) if initial else ""
        Df_init = str(initial["Df_init"]) if initial else ""
        eta_cstar = str(initial["eta_cstar"]) if initial else ""
        eta_nozzle = str(initial["eta_nozzle"]) if initial else ""
        OF_def = str(initial["OF_def"]) if initial else ""
        Pt_init = str(initial["Ptank_init"]) if initial else ""
        rho_ox = str(initial["rho_ox_init"]) if initial else ""
        rho_f = str(initial["rho_f_start"]) if initial else ""
        a_ox = str(initial["a_ox"]) if initial else ""
        n_ox = str(initial["n_ox"]) if initial else ""

        Kstar = str(results_parsed["Kstar"]) if results_parsed else ""
        epsilon = str(results_parsed["epsilon"]) if results_parsed else ""
        Lf = str(results_parsed["Lf"]) if results_parsed else ""
        mdot = str(results_parsed["mdot"]) if results_parsed else ""
        F = str(results_parsed["F"]) if results_parsed else ""
        Dt = str(results_parsed["Dt"]) if results_parsed else ""

        # 入力欄の定義
        Pc_box = ft.TextField(label="燃焼室圧力 Pc [MPa]", value=Pc_def, width=150)
        Df_box = ft.TextField(label="初期ポート径 Df [m]", value=Df_init, width=150)
        OF_box = ft.TextField(label="初期OF比", value=OF_def, width=150)
        eta_cstar_box = ft.TextField(label="C*効率", value=eta_cstar, width=150)
        eta_nozzle_box = ft.TextField(label="ノズル効率", value=eta_nozzle, width=150)


        Kstar_box = ft.TextField(label="K*", value=Kstar, width=150)
        epsilon_box = ft.TextField(label="膨張比 ε", value=epsilon, width=150)
        Lf_box = ft.TextField(label="燃焼長 Lf [m]", value=Lf, width=150)
        mdot_box = ft.TextField(label="推進剤流量 mdot [kg/s]", value=mdot, width=150)
        F_box = ft.TextField(label="初期推力F [N]", value=F, width=150)
        Dt_box = ft.TextField(label="スロート径 [m]", value=Dt, width=150)

        # タンク容積と最終酸化剤圧力の入力欄
        tank_volume_input = ft.TextField(label="タンク容積 [m³]", width=150)
        initial_pressure_input = ft.TextField(label="初期酸化剤圧力 [MPa]", value=Pt_init,width=150)
        rho_ox_input = ft.TextField(label="初期酸化剤密度(DBから自分で算出してください) [kg/s]", value=rho_ox,width=150)
        final_pressure_input = ft.TextField(label="最終酸化剤圧力 [MPa]", width=150)

        csv_download_button = ft.ElevatedButton(
            text="CSV出力 ⬇",
            icon=ft.Icons.DOWNLOAD,
            visible=False,
            on_click=lambda _: None,
        )

        def get_csv_download_link(evolution_result):
            print("output")
            filename = f"result.csv"
            with open(filename, "w", newline="", encoding="utf-8") as file:
                writer = csv.writer(file, quoting=csv.QUOTE_NONE)
                writer.writerows(evolution_result)

        # 関数に放り込む部分
        sim = RocketSimulation()

        def on_run_simulation(e):
            try:
                # 各入力欄から値を取得
                Pc = float(Pc_box.value)
                Df = float(Df_box.value)
                OF = float(OF_box.value)
                eta_cstar = float(eta_cstar_box.value)
                eta_nozzle = float(eta_nozzle_box.value)
                Kstar = float(Kstar_box.value)
                epsilon = float(epsilon_box.value)
                Lf = float(Lf_box.value)
                mdot = float(mdot_box.value)
                V_tank = float(tank_volume_input.value)
                P_init = float(initial_pressure_input.value)
                P_final = float(final_pressure_input.value)
                F_init = float(F_box.value)
                Dt = float(Dt_box.value)
                rho_ox = float(rho_ox_input.value)
                rho_f = initial["rho_f_start"]
                a_ox = initial["a_ox"]
                n_ox = initial["n_ox"]

                # RocketSimulation呼び出し
                (
                    time_ms,
                    F_arr,
                    F_fte_arr,
                    OF_arr,
                    Cstar_arr,
                    Pc_arr,
                    Pt_arr,
                    evolution_result,
                ) = sim.integration_simulation(
                    Pc=Pc,
                    Df=Df,
                    OF=OF,
                    eta_cstar=eta_cstar,
                    eta_nozzle=eta_nozzle,
                    Kstar=Kstar,
                    epsilon=epsilon,
                    Lf=Lf,
                    mdot=mdot,
                    V_tank=V_tank,
                    P_init=P_init,
                    P_final=P_final,
                    rho_ox=rho_ox,
                    rho_fuel=rho_f,
                    a=a_ox,
                    n=n_ox,
                    F=F_init,
                    Dt=Dt,
                )

                # 結果表示（仮）
                evolution_output.value = f"✅ 計算完了"
            except Exception as ex:
                evolution_output.value = f"⚠️ 計算エラー: {ex}"
                print(ex)

            def on_csv_download_click(e):
                csv_data_url = get_csv_download_link(evolution_result)
                page.launch_url(csv_data_url)

            csv_download_button.on_click = on_csv_download_click
            csv_download_button.visible = True
            results_graph_image.src_base64 = sim.get_evolution_plot_base64(
                time_ms, F_arr, F_fte_arr, OF_arr, Cstar_arr, Pc_arr, Pt_arr
            )
            results_graph_image.visible = True
            page.update()

        run_button = ft.ElevatedButton(
            text="時間発展計算 ▶", on_click=on_run_simulation
        )
        evolution_output = ft.Text("🕒 時間発展シミュレーション（仮表示）")

        return ft.View(
            route="/evolution",
            controls=[
                ft.Text("時間発展ページ", size=20, weight=ft.FontWeight.BOLD),
                ft.Row(
                    controls=[
                        # 1列目
                        ft.Column(
                            [
                                ft.Text("初期状態パラメータ①："),
                                Pc_box,
                                Df_box,
                                OF_box,
                                eta_cstar_box,
                                eta_nozzle_box,
                            ],
                            spacing=10,
                        ),
                        # 2列目
                        ft.Column(
                            [
                                ft.Text("初期状態パラメータ②："),
                                Kstar_box,
                                epsilon_box,
                                Lf_box,
                                mdot_box,
                                F_box,
                            ],
                            spacing=10,
                        ),
                        # 3列目
                        ft.Column(
                            [
                                ft.Text("初期状態パラメータ③："),
                                Dt_box,
                                tank_volume_input,
                                initial_pressure_input,
                                rho_ox_input, 
                                final_pressure_input,
                            ],
                            spacing=10,
                        ),
                        # ✅ 4列目：グラフ表示
                        ft.Column(
                            [ft.Text("時間発展グラフ："), results_graph_image],
                            spacing=10,
                        ),
                    ],
                    alignment=ft.MainAxisAlignment.START,
                    vertical_alignment=ft.CrossAxisAlignment.START,
                ),
                ft.Row(
                    controls=[run_button, csv_download_button, evolution_output],
                    alignment=ft.MainAxisAlignment.START,
                    vertical_alignment=ft.CrossAxisAlignment.START,
                ),
                ft.TextButton("◀ 戻る", on_click=lambda _: page.go("/")),
            ],
        )


ft.app(target=main)
