import numpy as np
import cea

ATM = 1.01325  # bar

class RocketCEA:
    def __init__(self, fuel="MMA", oxidizer="N2O", Pc_MPa=3.5, OF=6.0, n_frz = 2, epsilon=None):
        """
        ABS fuel uses fixed mass fractions:
        ABS1 = C3.85 H4.85 N0.43 → A=0.399, B=0.474, S=0.128
        ref https://www.sciencedirect.com/science/article/pii/S1270963825012209
        ref2 https://www.mdpi.com/2076-3417/16/16/8063#sec3dot3-applsci-16-08063 section 3.3.1 
        """
        self.fuel = fuel
        self.oxidizer = oxidizer
        self.Pc_bar = Pc_MPa * 10.0  # MPa → bar
        self.OF = OF
        self.n_frz = n_frz
        self.epsilon = epsilon

        # ============================
        # Fuel definition
        # ============================
        if self.fuel == "ABS":
            # ABS fuel → 3 monomers + oxidizer
            self.reac_names = ["ACRYLONITRILE", "BUTADIENE13", "STYRENE", oxidizer]

            # ABS1 mass fractions
            wA, wB, wS = 0.399, 0.473, 0.128

            self.fuel_weights = np.array([wA, wB, wS, 0.0])
            self.oxid_weights = np.array([0.0, 0.0, 0.0, 1.0])

            # Reactant temperatures
            self.T_reactants = np.array([300.0, 300.0, 300.0, 300.0])

        else:
            # MMA or other single-species fuel
            # MMAはthermo databaseもMMAで登録しているので，そのままdbに渡す
            self.reac_names = [fuel, oxidizer]
            self.fuel_weights = np.array([1.0, 0.0])
            self.oxid_weights = np.array([0.0, 1.0])
            self.T_reactants = np.array([300.0, 300.0])

        # Mixtures
        self.reac = cea.Mixture(self.reac_names)
        self.prod = cea.Mixture(self.reac_names, products_from_reactants=True)

        # Solver
        self.solver = cea.RocketSolver(self.prod, reactants=self.reac)
        self.solution = cea.RocketSolution(self.solver)

    def compute_pressure_ratio_for_exit_atm(self):
        """
        Pc / Pe = pressure ratio
        Pe = 1 atm = 1.01325 bar
        """
        Pe_bar = ATM
        return self.Pc_bar / Pe_bar

    def run(self):
        # Convert OF → reactant weights
        weights = self.reac.of_ratio_to_weights(
            self.oxid_weights, self.fuel_weights, self.OF
        )

        # Chamber enthalpy
        hc = self.reac.calc_property(
            cea.ENTHALPY, weights, self.T_reactants
        ) / cea.R

        # Decide pressure ratio or area ratio
        if self.epsilon is None:
            # Exit pressure = 1 atm
            pi_p = [self.compute_pressure_ratio_for_exit_atm()]
            subar = [1.01]      # throat　set 1.01 to avoid diverge
            supar = None       # not used
        else:
            # User-specified area ratio
            pi_p = None       # not used
            subar = [1.01]
            supar = [self.epsilon]

        # Solve
        self.solver.solve(
            self.solution,
            weights,
            self.Pc_bar,
            pi_p,
            subar=subar,
            supar=supar,
            hc=hc,
            iac=True,
            n_frz=self.n_frz
        )

        return self.format_output(self.epsilon)

    def format_output(self, epsilon):
        sol = self.solution
        if epsilon is None:
            # Exit pressure = 1 atm
            self.exit = 2
        else:
            self.exit = -1
        return {
            "chamber": {
                "T": sol.T[0],
                "P": sol.P[0],
                "rho": sol.density[0],
                "gamma": sol.gamma_s[0],
                "MW": sol.MW[0],
                "Cp": sol.cp_eq[0],
                "H": sol.enthalpy[0],
                "S": sol.entropy[0],
            },
            "throat": {
                "T": sol.T[1],
                "P": sol.P[1],
                "rho": sol.density[1],
                "Mach": sol.Mach[1],
                "gamma": sol.gamma_s[1],
                "MW": sol.MW[1],
                "Cp": sol.cp_eq[1],
                "H": sol.enthalpy[1],
                "S": sol.entropy[1],
            },
            "exit": {
                "T": sol.T[self.exit],
                "P": sol.P[self.exit],
                "rho": sol.density[self.exit],
                "Mach": sol.Mach[self.exit],
                "epsilon": sol.ae_at[self.exit],
                "gamma": sol.gamma_s[self.exit],
                "MW": sol.MW[self.exit],
                "Cp": sol.cp_eq[self.exit],
                "H": sol.enthalpy[self.exit],
                "S": sol.entropy[self.exit],
            },
            "performance": {
                "Cstar": sol.c_star[self.exit],
                "Cf": sol.coefficient_of_thrust[self.exit],
                "Isp": sol.Isp[self.exit],
                "Isp_vac": sol.Isp_vacuum[self.exit],
            },
            "mole_fractions": sol.mole_fractions,
        }