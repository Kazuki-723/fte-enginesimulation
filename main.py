from pythonfiles.rocket_simulation import RocketSimulation

if __name__ == '__main__':
    sim = RocketSimulation()

    # 初期状態のみ検証
    #sim.run_initialstate()

    # 全部回す
    sim.run_simulation()
