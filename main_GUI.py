import flet as ft
import matplotlib.pyplot as plt
import io
import base64

def main(page: ft.Page):
    page.title = "y = ax + b グラフ描画"
    page.scroll = ft.ScrollMode.AUTO

    # 入力フィールド
    a_input = ft.TextField(label="a の値", width=150)
    b_input = ft.TextField(label="b の値", width=150)
    graph_image = ft.Image(visible=False)  # 初期状態では非表示(base64エンコードエラー回避)

    def plot_graph(e):
        try:
            a = float(a_input.value)
            b = float(b_input.value)
        except ValueError:
            page.snack_bar = ft.SnackBar(ft.Text("a, b は数値で入力してください"))
            page.snack_bar.open = True
            page.update()
            return

        # グラフ描画
        x = [i for i in range(-10, 11)]
        y = [a * xi + b for xi in x]

        plt.figure()
        plt.plot(x, y, label=f"y = {a}x + {b}")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.grid(True)
        plt.legend()

        # 画像として保存してbase64変換
        buf = io.BytesIO()
        plt.savefig(buf, format="png")
        buf.seek(0)
        encoded = base64.b64encode(buf.read()).decode("utf-8")
        buf.close()
        plt.close()

        # Fletで表示
        graph_image.src_base64 = encoded
        graph_image.visible = True
        page.update()

    # UI構成
    page.add(
        ft.Row([a_input, b_input, ft.ElevatedButton("グラフ描画", on_click=plot_graph)]),
        graph_image
    )

ft.app(target=main)