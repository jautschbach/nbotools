

#Start import from here
import pandas as pd
import numpy as np
from itertools import product

from my_func import select_cols_only
from visualizer import save_loc_data


from sep_mat import ao_basis_list_dict


from sep_mat import is_open
from sep_mat import aoX_trans
from sep_mat import pnaoX_trans
from sep_mat import opera_matrix

from pop_analys import population_analy
import time



if is_open:

    SAO_alpha = pd.DataFrame(opera_matrix["OVERLAP"]).values
    FAO_alpha = pd.DataFrame(opera_matrix["FOCK_ALPHA"]).values
    DM_alpha = pd.DataFrame(opera_matrix["DENSITY_ALPHA"]).values
    AOMO_alpha = pd.DataFrame(aoX_trans["AOMOs_ALPHA"]).values
    PNAOMO_alpha = pd.DataFrame(pnaoX_trans["PNAOMOs_ALPHA"]).values
    SPNAO_alpha = pd.DataFrame(opera_matrix["PNAO OVERLAP"]).values
    FPNAO_alpha = pd.DataFrame(opera_matrix['PNAO Alpha FOCK']).values
   
    num_electron_alpha = ((np.trace(DM_alpha @ SAO_alpha)))
    num_electron_alpha = int(num_electron_alpha.round(2))
    num_occ_orb_alpha = int(num_electron_alpha)   
    # print(num_electron_alpha, num_occ_orb_alpha)
    X_alpha = pd.DataFrame(aoX_trans['AOPNAOs_ALPHA'])   
    X_alpha =X_alpha.values    

    
    
    
    SAO_beta = pd.DataFrame(opera_matrix["OVERLAP"]).values
    FAO_beta = pd.DataFrame(opera_matrix["FOCK_BETA"]).values
    DM_beta = pd.DataFrame(opera_matrix["DENSITY_BETA"]).values
    AOMO_beta = pd.DataFrame(aoX_trans["AOMOs_BETA"]).values
    PNAOMO_beta = pd.DataFrame(pnaoX_trans["PNAOMOs_BETA"]).values
    SPNAO_beta = pd.DataFrame(opera_matrix["PNAO OVERLAP"]).values
    FPNAO_beta = pd.DataFrame(opera_matrix['PNAO Beta FOCK']).values
    
    num_electron_beta = ((np.trace(DM_beta @ SAO_beta)))
    num_electron_beta = int(num_electron_beta.round(2))
    num_occ_orb_beta = int(num_electron_beta)   
    # print(num_electron_beta, num_occ_orb_beta)
    # X_beta = T.values
    X_beta = X_alpha
    
    U_mat = np.linalg.inv(X_alpha)



else:    
    FAO = pd.DataFrame(opera_matrix["FOCK"]).values
    density = pd.DataFrame(opera_matrix["DENSITY"]).values
    SAO = pd.DataFrame(opera_matrix["OVERLAP"]).values
    FPNAO = pd.DataFrame(opera_matrix["PNAO FOCK"]).values
    SPNAO = pd.DataFrame(opera_matrix['PNAO OVERLAP']).values
    AOMO = pd.DataFrame(aoX_trans["AOMOs"]).values
    PNAOMOs = pd.DataFrame(pnaoX_trans["PNAOMOs"]).values
    X = pd.DataFrame(aoX_trans['AOPNAOs'])
    U_mat = np.linalg.inv(X.values)
       
    num_electron = ((np.trace(density @ SAO)))
    num_electron = int(num_electron.round(2))
    num_occ_orb = int(num_electron / 2)
    
    

def loc_all_occ(del_cmo, over_mat, fock_mat, n_occ, bas):
    
    
    # overlap = S
    # c_del = cmo 
    overlap = over_mat
    c_del = del_cmo
    cmo = del_cmo
    fock = fock_mat
    n_occ_oribitals = n_occ
    
    basis = bas 
    n_basis_fn = len(pd.DataFrame(cmo).index)
    t_sc = overlap @ c_del[:, :n_occ_oribitals]
    c = c_del[:, :n_occ_oribitals]
    sc = t_sc[:, :n_occ_oribitals]
    tol = 1e-8
    t_range = n_occ_oribitals
    s_range = n_occ_oribitals
    natom = len(basis)
    gamma_tol = 1e-10
    c_initial = c_del[:, :n_occ_oribitals]
    sc_initial = t_sc[:, :n_occ_oribitals] 
    tmp = np.empty((n_basis_fn, n_occ_oribitals))
    tmp_sc = np.empty((n_basis_fn, n_occ_oribitals))
    nrot = 1
    convergence_threshold = 1e-6 
    
    
    nrot = 1
    prev_result = None  # Initialize prev_result
    
    for itera in range(1, 2000):
        random_s = np.random.choice(s_range, size=s_range, replace=False)
        
        current_result = np.diag(c.T @ fock @ c)
        
        # print("----------", current_result)
        
        # Check for convergence
        if prev_result is not None and np.allclose(current_result, prev_result, atol=1e-12):
            print("Localization converged! after ", itera, " iteration...")
            break
        if itera == 1999:
            print("Localization not converged after 1999 iterations")
            break
        
        prev_result = current_result
        
        # Generate pairs of s and t where t != s
        pairs = [(s, t) for s, t in product(random_s, range(t_range)) if t != s]
        
        for s, t in pairs:
            ast = 0
            bst = 0
            gamma_max = 0
        
            for a in range(natom):
                qast = 0
                qat = 0
                qas = 0
        
                bflo = basis[a]['bflo']
                bfhi = basis[a]['bfhi']
        
                for u in range(bflo, bfhi + 1):
                    qast = qast + (c[u, t] * sc[u, s] + c[u, s] * sc[u, t]) * 0.5
                    qas = qas + c[u, s] * sc[u, s]
                    qat = qat + c[u, t] * sc[u, t]
        
                ast = ast + (qast ** 2) - 0.25 * ((qas - qat) ** 2)
                bst = bst + qast * (qas - qat)
        
            gamma = 0.25 * np.arccos((-1 * ast) / (np.sqrt(ast ** 2 + bst ** 2)))
            gamma = abs(gamma) * (bst / abs(bst))
            gamma_max = max(gamma_max, abs(gamma))
        
            if abs(gamma) > gamma_tol:   
                nrot = nrot + 1
                cosg = np.cos(gamma)
                sing = np.sin(gamma)
                tmp[:, 0] = c[:, s] * cosg + c[:, t] * sing
                tmp[:, 1] = c[:, t] * cosg - c[:, s] * sing
                
                c[:, s] = tmp[:, 0]
                c[:, t] = tmp[:, 1]
                
                tmp_sc[:, 0] = sc[:, s] * cosg + sc[:, t] * sing
                tmp_sc[:, 1] = sc[:, t] * cosg - sc[:, s] * sing
                
                sc[:, s] = tmp_sc[:, 0]
                sc[:, t] = tmp_sc[:, 1]
        
    
    T = sc_initial.T @ c
    
    ###### Begin sorting of c and sc based on the index of a sorted energy calculation
    energy_c = np.diag(c.T @ fock @ c)
    sort_energ_ind = np.argsort(energy_c)
    # print(sort_energ_ind)
    sorted_c = c[: , sort_energ_ind]
    c = sorted_c
    
    sorted_sc = sc[: , sort_energ_ind]
    sc = sorted_sc

    loc_energy = np.diag(sorted_c.T @ fock @ sorted_c)
    print(f"""Energy of localized orbitals:
{loc_energy}\n""")

    # population_analy(sorted_c, overlap, bas,  fock)
          
 
    return sorted_c




# Vectorized calculation of qast, qas, and qat
def calculate_q_values(c, sc, bflo, bfhi, t, s):
    c_slice = c[bflo:bfhi+1, [t, s]]
    sc_slice = sc[bflo:bfhi+1, [t, s]]
    
    qast = np.sum((c_slice[:, 0] * sc_slice[:, 1] + c_slice[:, 1] * sc_slice[:, 0]) * 0.5)
    qas = np.sum(c_slice[:, 1] * sc_slice[:, 1])
    qat = np.sum(c_slice[:, 0] * sc_slice[:, 0])
    
    return qast, qas, qat



def loc_select_orb_arg(del_cmo, over_mat, fock_mat, bas): #def loc_select_orb_arg(del_cmo, over_mat, fock_mat, bas):
    
    sel_orb = input("""Select orbitals  1,2,6,12 - 19:
""")
    index_arguments = sel_orb.split(",")
    # print(index_arguments)
    

    overlap = over_mat
    fock = fock_mat
    basis = bas
    # print(basis)
        
    cmo = del_cmo
    # print(np.diag(cmo.T @ fock @ cmo))
    cmo = pd.DataFrame(cmo)
    cmo.index = cmo.index + 1
    cmo.columns = cmo.columns + 1
    
    sel_cmo = select_cols_only(cmo,*index_arguments)
   
    
    
    c = sel_cmo.values
    col_list = pd.DataFrame(sel_cmo).columns.to_list()
    
    orbind = [x - 1 for x in col_list]
    
    # print(orbind, "---------------------")
    sc = overlap @ c
    num_orb_to_loc = len(pd.DataFrame(c).columns)
    
    # print(num_orb_to_loc, n_basis_fn)
    
    t_range = num_orb_to_loc
    s_range = num_orb_to_loc
    natom = len(basis)
    
    
    gamma_tol = 1e-10    
    nrot = 1
    prev_result = None  # Initialize prev_result
    

  
    stime = time.time()
    for itera in range(1, 2000):
        random_s = np.random.choice(s_range, size=s_range, replace=False)
        pairs = [(s, t) for s, t in product(random_s, range(t_range)) if t != s]
        
        current_result = np.diag(c.T @ fock @ c)
        
        # Check for convergence
        if prev_result is not None and np.allclose(current_result, prev_result, atol=1e-9):
            print("\nLocalization converged! after ", itera, " iteration...\n")
            break
        
        if itera == 1999:
            print("\nLocalization not converged after 1999 iterations\n")
            break
        
        prev_result = current_result

#Start with pairwise generation
        # Generate pairs of s and t where t != s
        
        # Initialize rotation counter
        nrot = 0
        
        # Main loop through pairs
        for s, t in pairs:
            ast = 0
            bst = 0
            gamma_max = 0
        
            for a in range(natom):
                bflo = basis[a]['bflo']
                bfhi = basis[a]['bfhi']
        
                qast, qas, qat = calculate_q_values(c, sc, bflo, bfhi, t, s)
        
                diff_as_at = qas - qat
                ast += qast * qast - 0.25 * diff_as_at**2
                bst += qast * diff_as_at
        
            # Avoid numerical instability in arccos
            cos_arg = -ast / np.sqrt(ast**2 + bst**2)
            cos_arg = np.clip(cos_arg, -1.0, 1.0)
            gamma = 0.25 * np.arccos(cos_arg)
            gamma *= np.sign(bst)
        
            gamma_max = max(gamma_max, abs(gamma))
        
            if abs(gamma) > gamma_tol:
                nrot += 1
                cosg = np.cos(gamma)
                sing = np.sin(gamma)
        
                rotation_matrix = np.array([[cosg, -sing], [sing, cosg]])
        
                c[:, [s, t]] = c[:, [s, t]] @ rotation_matrix
                sc[:, [s, t]] = sc[:, [s, t]] @ rotation_matrix
                        
    
    etime = time.time()
    
    print(f"Converged in {etime - stime}")

    ###### Begin sorting of c and sc based on the index of a sorted energy calculation
    energy_c = np.diag(c.T @ fock @ c)
    # print(energy_c)
    sort_energ_ind = np.argsort(energy_c)
    # print(sort_energ_ind)
    sorted_c = c[: , sort_energ_ind]
    c = sorted_c
    
    sorted_sc = sc[: , sort_energ_ind]
    sc = sorted_sc
    
    loc_energy = np.diag(sorted_c.T @ fock @ sorted_c)
    # print("\n Printing bas: ", bas)
    print(f"""Energy of localized orbitals: 
{loc_energy}\n""")
    

    return sorted_c #del_cmo
    



def get_localized_orbitals(basis_type, loc_func, convert_to_ao=False):
    """
    Localizes orbitals based on specified basis type and localization function.
    """
    if is_open:
        print(f"Localizing Alpha orbitals in {basis_type} basis")
        AOMO_alpha_basis = AOMO_alpha if basis_type == "AO" else PNAOMO_alpha
        SAO_alpha_basis = SAO_alpha if basis_type == "AO" else SPNAO_alpha
        FAO_alpha_basis = FAO_alpha if basis_type == "AO" else FPNAO_alpha
        
        # Conditional argument passing for loc_func
        if loc_func == loc_all_occ:
            res_alpha = loc_func(AOMO_alpha_basis, SAO_alpha_basis, FAO_alpha_basis, num_occ_orb_alpha, ao_basis_list_dict)
        else:
            res_alpha = loc_func(AOMO_alpha_basis, SAO_alpha_basis, FAO_alpha_basis, ao_basis_list_dict)
        
        print(f"\nLocalized alpha orbitals printed in {basis_type} basis\n")
        print(pd.DataFrame(res_alpha))

        res_alpha_ao = np.linalg.inv(U_mat) @ res_alpha if convert_to_ao else res_alpha
        population_analy(res_alpha_ao, SAO_alpha, ao_basis_list_dict, FAO_alpha)

        print(f"\nLocalizing Beta orbitals in {basis_type} basis")
        AOMO_beta_basis = AOMO_beta if basis_type == "AO" else PNAOMO_beta
        SAO_beta_basis = SAO_beta if basis_type == "AO" else SPNAO_beta
        FAO_beta_basis = FAO_beta if basis_type == "AO" else FPNAO_beta
        
        # Conditional argument passing for loc_func
        if loc_func == loc_all_occ:
            res_beta = loc_func(AOMO_beta_basis, SAO_beta_basis, FAO_beta_basis, num_occ_orb_beta, ao_basis_list_dict)
        else:
            res_beta = loc_func(AOMO_beta_basis, SAO_beta_basis, FAO_beta_basis, ao_basis_list_dict)

        print(f"\nLocalized beta orbitals printed in {basis_type} basis\n")
        print(pd.DataFrame(res_beta))

        res_beta_ao = np.linalg.inv(U_mat) @ res_beta if convert_to_ao else res_beta
        population_analy(res_beta_ao, SAO_beta, ao_basis_list_dict, FAO_beta)

        return res_alpha_ao, res_beta_ao
    else:
        print(f"Localizing Orbitals in {basis_type} basis")
        AOMO_basis = AOMO if basis_type == "AO" else PNAOMOs
        SAO_basis = SAO if basis_type == "AO" else SPNAO
        FAO_basis = FAO if basis_type == "AO" else FPNAO
        
        # Conditional argument passing for loc_func
        if loc_func == loc_all_occ:
            res = loc_func(AOMO_basis, SAO_basis, FAO_basis, num_occ_orb, ao_basis_list_dict)
        else:
            res = loc_func(AOMO_basis, SAO_basis, FAO_basis, ao_basis_list_dict)

        print(f"\nLocalized orbitals printed in {basis_type} basis\n")
        print(pd.DataFrame(res))

        res_ao = np.linalg.inv(U_mat) @ res if convert_to_ao else res
        population_analy(res_ao, SAO, ao_basis_list_dict, FAO)

        return res_ao

def save_orbitals(localized_orbitals, basis_type):
    """
    Saves localized orbitals to a file.
    """
    file_name = input("file_name without extension: ")
    file_name = file_name + '.40'
    if is_open:
        text = f"Localized MOs in full {basis_type} Basis"
        print(f"\nSaving localized orbitals in {basis_type} basis...")
        save_loc_data(pd.DataFrame(localized_orbitals[0]), file_name, text, pd.DataFrame(localized_orbitals[1]))
    else:
        text = f"Localized MOs in full {basis_type} Basis"
        print(f"\nSaving localized orbitals in {basis_type} basis...")
        save_loc_data(pd.DataFrame(localized_orbitals), file_name, text)


def localize_and_save(basis_type, loc_func, convert_to_ao=False):
    """
    Localizes orbitals and saves them based on user input.
    """
    localized_orbitals = get_localized_orbitals(basis_type, loc_func, convert_to_ao)
    save_option = input("Save localized orbitals to file (y/n): ").lower()
    if save_option == 'y':
        basis_text = "AO" if basis_type == "AO" else "pNAO"
        save_orbitals(localized_orbitals, basis_text)


def pm_run():
    while True:
        user_input = int(input("""
1. Localize only occupied orbitals
2. Localize selected range of orbitals: 
-1. Back                     
"""))

        if user_input == 1:
            while True:
                user_input = int(input("""
1. Localize AOMOs in AO basis. Returns localized MOs in AO basis
2. Localize pNAOMOs in pNAO basis. Returns localized MOs in AO basis
-1. Back             
"""))

                if user_input == 1:
                    localize_and_save("AO", loc_all_occ)
                elif user_input == 2:
                    localize_and_save("pNAO", loc_all_occ, convert_to_ao=True)
                elif user_input == -1:
                    break

        elif user_input == 2:
            while True:
                user_input = int(input("""
1. Localize selected AOMOs in AO basis. Returns localized selected MOs in AO basis
2. Localize selected pNAOMOs in pNAO basis. Returns localized selected MOs in AO basis
-1. Back             
"""))

                if user_input == 1:
                    localize_and_save("AO", loc_select_orb_arg)
                elif user_input == 2:
                    localize_and_save("pNAO", loc_select_orb_arg, convert_to_ao=True)
                elif user_input == -1:
                    break

        elif user_input == -1:
            return
