import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import io, base64

class IterationLogger:
    """
    GUI専用クラス
    GUI側で実行後に結果をmatplotlibで表示する機能のためのもの
    エンコーダーはfletとの相性の問題かAgg
    """
    def __init__(self):
        self.log = []

    def append(self, j, F, mdot, Pe, epsilon):
        self.log.append({
            "j": j,
            "F": F,
            "mdot": mdot,
            "Pe": Pe,
            "epsilon": epsilon
        })

    def get_lists(self):
        return {
            "j": [entry["j"] for entry in self.log],
            "F": [entry["F"] for entry in self.log],
            "mdot": [entry["mdot"] for entry in self.log],
            "Pe": [entry["Pe"] for entry in self.log],
            "epsilon": [entry["epsilon"] for entry in self.log],
        }

    # 初期値側のグラフ表示
    def get_base64_plot(self, Dovalue, cdvalue):
        data = self.get_lists()
        print("plot start")
        fig, axs = plt.subplots(3, 2, figsize=(12, 10))

        # 1段目　
        # Thrustの収束履歴
        axs[0, 0].plot(data["j"], data["F"], marker='o')
        axs[0, 0].set_title("Thrust [N]")

        # 総流量の収束履歴
        axs[0, 1].plot(data["j"], data["mdot"], marker='o')
        axs[0, 1].set_title("Mass Flow Rate [kg/s]")

        # 2段目
        # 出口圧の収束履歴
        axs[1, 0].plot(data["j"], data["Pe"], marker='o')
        axs[1, 0].set_title("Exit Pressure [MPa]")

        # 開口比の収束履歴
        axs[1, 1].plot(data["j"], data["epsilon"], marker='o')
        axs[1, 1].set_title("Expansion Ratio [-]")

        # 3段目：Do vs Cd
        # K* =const.として，オリフィス径と流量係数の関係計算
        axs[2, 0].plot(Dovalue, cdvalue, color='blue')
        axs[2, 0].set_title("Do vs Cd")
        axs[2, 0].set_xlabel("Do [m]")
        axs[2, 0].set_ylabel("Cd [-]")

        # 3段目右側は非表示にせず、空白として残す（警告回避）
        axs[2, 1].axis("off")

        # 軸ラベルとグリッド
        for i in range(2):
            for j in range(2):
                axs[i, j].set_xlabel("Iteration")
                axs[i, j].grid(True)
        axs[2, 0].grid(True)

        # 基本設定，データエンコード
        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format="png")
        buf.seek(0)
        encoded = base64.b64encode(buf.read()).decode("utf-8")
        buf.close()
        plt.close(fig)
        return encoded

    # 時間発展側のグラフ表示
    def plot_time_series(time_ms, F_arr, F_fte_arr, OF_arr, Cstar_arr, Pc_arr, Pt_arr):
        fig, axs = plt.subplots(2, 2, figsize=(10, 8))

        # 1段目
        # 推力は二通りの計算方法が存在しているため両方表示
        axs[0, 0].plot(time_ms, F_arr, marker='o', label = "F_CF")
        axs[0, 0].plot(time_ms, F_fte_arr, marker='x', label = "F_mdot")
        axs[0, 0].set_xlabel("Time [ms]")
        axs[0, 0].set_ylabel("Thrust [N]")

        # OFの時間履歴
        axs[0, 1].plot(time_ms, OF_arr, marker='o')
        axs[0, 1].set_xlabel("Time [ms]")
        axs[0, 1].set_ylabel("O/F [-]")

        # 2段目
        # Cstarの時間履歴
        axs[1, 0].plot(time_ms, Cstar_arr, marker='o')
        axs[1, 0].set_xlabel("Time [ms]")
        axs[1, 0].set_ylabel("Charactaristic velosity [m/s]")

        # タンク圧，燃焼室圧の時間履歴
        axs[1, 1].plot(time_ms, Pc_arr, marker='o')
        axs[1, 1].plot(time_ms, Pt_arr, marker='o')
        axs[1, 1].set_xlabel("Time [ms]")
        axs[1, 1].set_ylabel("Pressure [MPa]")

        # 基本設定，データエンコード
        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format="png")
        buf.seek(0)
        encoded = base64.b64encode(buf.read()).decode("utf-8")
        buf.close()
        plt.close(fig)
        return encoded
