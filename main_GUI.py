import flet as ft
from inputprograms.rocket_simulation import RocketSimulation

def main(page: ft.Page):
    page.title = "Rocket Simulation GUI"
    page.scroll = ft.ScrollMode.AUTO

    # 入力フィールド群
    inputs = {
        "F_req": ft.TextField(label="要求推力 [N]", width=150),
        "Pc_def": ft.TextField(label="初期燃焼室圧力 [MPa]", width=150),
        "OF_def": ft.TextField(label="初期O/F比", width=150),
        "mdot_new": ft.TextField(label="初期流量 [kg/s]", width=150),
        "eta_cstar": ft.TextField(label="C*効率", width=150),
        "eta_nozzle": ft.TextField(label="ノズル効率", width=150),
    }

    result_text = ft.Text()
    graph_image = ft.Image(visible=False, expand=True)

    def run_simulation(e):
        try:
            values = {k: float(inputs[k].value) for k in inputs}
        except ValueError:
            result_text.value = "⚠️ 全ての値を数値で入力してください"
            page.update()
            return

        sim = RocketSimulation()
        output = sim.initial_convergence(**values)
        result_text.value = output

        # グラフ表示
        graph_image.src_base64 = sim.get_iteration_plot_base64()
        graph_image.visible = True
        page.update()

    # 左側：入力群
    input_column = ft.Column(
        controls=list(inputs.values()) + [ft.ElevatedButton("収束計算", on_click=run_simulation), result_text],
        spacing=10
    )

    # 右側：グラフ
    graph_column = ft.Column(
        controls=[graph_image],
        spacing=10,
        alignment=ft.MainAxisAlignment.START
    )   


    # 横並びに配置
    page.add(
        ft.Row(
            controls=[input_column, graph_column],
            alignment=ft.MainAxisAlignment.START,  # 横並びの左寄せ
            vertical_alignment=ft.CrossAxisAlignment.START
        )
    )


ft.app(target=main)