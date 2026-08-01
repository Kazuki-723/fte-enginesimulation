import numpy as np
import cantera as ct

def build_stoich_coeffs(gas):
    """元素×種の係数行列を一度だけ構築してキャッシュする"""
    ne = gas.n_elements
    ns = gas.n_species
    stoich = np.zeros((ne, ns))
    for i, elem in enumerate(gas.element_names):
        for j, sp in enumerate(gas.species_names):
            stoich[i, j] = gas.n_atoms(sp, elem)
    return stoich

def get_thermo_derivatives(gas, stoich_coeffs):
    '''Gets thermo derivatives based on shifting equilibrium.
    '''
    # unknowns for system with no condensed species:
    # dpi_i_dlogT_P (# elements)
    # dlogn_dlogT_P
    # dpi_i_dlogP_T (# elements)
    # dlogn_dlogP_T
    # total unknowns: 2*n_elements + 2

    # cantera cashe
    h_RT = gas.standard_enthalpies_RT
    X = gas.X
    MW = gas.mean_molecular_weight
    ne = gas.n_elements
    num_var = 2 * ne + 2

    coeff_matrix = np.zeros((num_var, num_var))
    right_hand_side = np.zeros(num_var)

    tot_moles = 1.0 / MW
    moles = X * tot_moles

    condensed = False

    # indices
    idx_dpi_dlogT_P = 0
    idx_dlogn_dlogT_P = idx_dpi_dlogT_P + ne
    idx_dpi_dlogP_T = idx_dlogn_dlogT_P + 1
    idx_dlogn_dlogP_T = idx_dpi_dlogP_T + ne

    # equations for derivatives with respect to temperature
    # first n_elements equations
    for k in range(ne):
        for i in range(ne):
            coeff_matrix[k,i] = np.sum(stoich_coeffs[k,:] * stoich_coeffs[i,:] * moles)
        coeff_matrix[k, ne] = np.sum(stoich_coeffs[k,:] * moles)
        right_hand_side[k] = -np.sum(stoich_coeffs[k,:] * moles * h_RT)

    # skip equation relevant to condensed species

    for i in range(ne):
        coeff_matrix[ne, i] = np.sum(stoich_coeffs[i, :] * moles)
    right_hand_side[ne] = -np.sum(moles * h_RT)

    # equations for derivatives with respect to pressure

    for k in range(ne):
        for i in range(ne):
            coeff_matrix[ne+1+k,ne+1+i] = np.sum(stoich_coeffs[k,:] * stoich_coeffs[i,:] * moles)
        coeff_matrix[ne+1+k, 2*ne+1] = np.sum(stoich_coeffs[k,:] * moles)
        right_hand_side[ne+1+k] = np.sum(stoich_coeffs[k,:] * moles)

    for i in range(ne):
        coeff_matrix[2*ne+1, ne+1+i] = np.sum(stoich_coeffs[i, :] * moles)
    right_hand_side[2*ne+1] = np.sum(moles)
    
    try:
        derivs = np.linalg.solve(coeff_matrix, right_hand_side)
    except np.linalg.LinAlgError:
        # print("Singular matrix detected, switching to least-squares solution.")
        derivs, *_ = np.linalg.lstsq(coeff_matrix, right_hand_side, rcond=None)


    dpi_dlogT_P = derivs[idx_dpi_dlogT_P : idx_dpi_dlogT_P + ne]
    dlogn_dlogT_P = derivs[idx_dlogn_dlogT_P]
    dpi_dlogP_T = derivs[idx_dpi_dlogP_T]
    dlogn_dlogP_T = derivs[idx_dlogn_dlogP_T]

    # dpi_dlogP_T is not used
    
    return dpi_dlogT_P, dlogn_dlogT_P, dlogn_dlogP_T

def get_thermo_properties(gas, stoich_coeffs, dpi_dlogT_P, dlogn_dlogT_P, dlogn_dlogP_T):
    '''Calculates specific heats, volume derivatives, and specific heat ratio.
    
    Based on shifting equilibrium for mixtures.
    '''

    # cantera cashe
    h_RT = gas.standard_enthalpies_RT
    cp_R = gas.standard_cp_R
    X = gas.X
    MW = gas.mean_molecular_weight
    ne = gas.n_elements

    tot_moles = 1.0 / MW
    moles = X * tot_moles
    
    spec_heat_p = ct.gas_constant * (
        np.sum([dpi_dlogT_P[i] * 
                np.sum(stoich_coeffs[i,:] * moles * h_RT) 
                for i in range(ne)
                ]) +
        np.sum(moles * h_RT) * dlogn_dlogT_P +
        np.sum(moles * cp_R) +
        np.sum(moles * h_RT**2)
        )
    
    dlogV_dlogT_P = 1 + dlogn_dlogT_P
    dlogV_dlogP_T = -1 + dlogn_dlogP_T
    
    spec_heat_v = (
        spec_heat_p + gas.P * gas.v / gas.T * dlogV_dlogT_P**2 / dlogV_dlogP_T
        )

    gamma = spec_heat_p / spec_heat_v
    gamma_s = -gamma/dlogV_dlogP_T
    
    return dlogV_dlogT_P, dlogV_dlogP_T, spec_heat_p, gamma_s