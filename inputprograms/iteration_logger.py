class IterationLogger:
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

    def get_base64_plot(self):
        import matplotlib.pyplot as plt
        import io, base64

        data = self.get_lists()
        fig, axs = plt.subplots(2, 2, figsize=(10, 8))

        axs[0, 0].plot(data["j"], data["F"], marker='o')
        axs[0, 0].set_title("Thrust [N]")

        axs[0, 1].plot(data["j"], data["mdot"], marker='o')
        axs[0, 1].set_title("Mass Flow Rate [kg/s]")

        axs[1, 0].plot(data["j"], data["Pe"], marker='o')
        axs[1, 0].set_title("Exit Pressure [MPa]")

        axs[1, 1].plot(data["j"], data["epsilon"], marker='o')
        axs[1, 1].set_title("Expansion Ratio [-]")

        for ax in axs.flat:
            ax.set_xlabel("Iteration")
            ax.grid(True)

        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format="png")
        buf.seek(0)
        encoded = base64.b64encode(buf.read()).decode("utf-8")
        buf.close()
        plt.close()
        return encoded