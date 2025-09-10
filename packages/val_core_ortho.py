
import pandas as pd
import numpy as np
from scipy.linalg import eigh
from my_func import get_user_input
from my_func import submatrix
from sep_mat import is_open
from sep_mat import opera_matrix, aoX_trans
from sep_mat import ao_basis
from sep_mat import list_to_blo_bi, ao_basis_list_dict
from pm_localization import loc_select_orb_arg
from visualizer import sub_to_full
from pop_analys import population_analy
from sep_mat import filename
import os
import fnmatch
import re

def find_nbo_out_file():
    pattern = f"{filename}.47*.out"
    matched_files = [os.path.join("./", file) for file in os.listdir("./") if fnmatch.fnmatch(file, pattern)]
    if len(matched_files) == 1:
        return matched_files[0]
    if not matched_files:
        print("No NBO output file found")
    else:
        print("Multiple files found...")
    return input("Enter NBO output file: ")

def extract_nao_data(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            header_lines = [i for i, line in enumerate(lines) if " NAO Atom No lang   Type(AO)    Occupancy" in line]
            if len(header_lines) == 1:
                header_line = header_lines[0]
            else:
                alpha_beta = input("Enter A for Alpha or B for Beta: ")
                if alpha_beta.lower() == 'a':
                    header_line = header_lines[1]
                elif alpha_beta.lower() == 'b':
                    header_line = header_lines[2]
                else:
                    raise ValueError("Invalid input. Please enter 'A' or 'B'.")
            nao_list = []
            pattern1 = re.compile(r'\((.*?)\)')
            pattern2 = re.compile(r'(\w+)\)')
            for line in lines[header_line+2:]:
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        nao = int(parts[0])
                        atom_info = parts[1] + parts[2]
                        val_core_ryd = parts[4][:3]
                        match1 = pattern1.search(parts[4])
                        if match1:
                            key = match1.group(1)
                            nao_list.append((atom_info, key, nao))
                        match2 = pattern2.match(parts[5])
                        if match2:
                            key = match2.group(1)
                            nao_list.append((atom_info, key, nao, val_core_ryd))
                    except ValueError:
                        break
            return nao_list
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found")
        return []
    except Exception as e:
        print(f"Error: {str(e)}")
        return []

def valence_core_index():
    file_path = find_nbo_out_file()
    nao_data = extract_nao_data(file_path)
    core_idx = [entry[2] for entry in nao_data if len(entry) == 4 and entry[3].lower() == 'cor']
    core_idx = ", ".join(map(str, core_idx))
    val_idx = [entry[2] for entry in nao_data if len(entry) == 4 and entry[3].lower() == 'val']
    val_idx = ", ".join(map(str, val_idx))
    return core_idx, val_idx

def search_nao_data(nao_data, search_criteria):
    results = {}
    for atom_orbitals in search_criteria.split(','):
        atom_orbitals = atom_orbitals.strip()
        if ':' in atom_orbitals:
            atom, orbitals = atom_orbitals.split(':')
            atom = atom.strip().lower()
            orbitals = [orb.strip().lower() for orb in orbitals.split()]
            atom_results = []
            for orbital in orbitals:
                orbital_results = [value for ai, k, value in nao_data if ai.lower() == atom and k.lower() == orbital]
                atom_results.extend(orbital_results)
            if atom_results:
                results[atom] = atom_results
    return results

def atom_orbital():
    file_path = find_nbo_out_file()
    nao_data = extract_nao_data(file_path)
    nao_data = [(entry[0], entry[1], entry[2]) for entry in nao_data]
    user_input = input("Enter search criteria (e.g., U1: 3d 2s 3p, Cl2: 3p 4f 3d): ")
    results = search_nao_data(nao_data, user_input)
    if results:
        all_matching_naos = set()
        for values in results.values():
            all_matching_naos.update(values)
        ordered_results = [str(nao) for _, _, nao in nao_data if nao in all_matching_naos]
        print("\nCombined result:")
        orbital_list = (', '.join(ordered_results))
    else:
        print("No matching data found")
    return orbital_list

def get_non_zero_positions(matrix, tol=1e-10):
    positions = []
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if abs(matrix[i, j]) > tol:
                positions.append((i, j, matrix[i, j]))
    return positions

def is_unit_matrix(matrix, tolerance=1e-10):
    if matrix.shape[0] != matrix.shape[1]:
        return False
    identity = np.eye(matrix.shape[0])
    return np.allclose(matrix, identity, atol=tolerance)

def save_vis_file(matrix, file_name):
    with open(file_name, 'w') as file:
        file.write(f"""NBO Tools: {file_name.replace('.40', '')}
Eigen vectors in Reduced Basis transformed to AO Basis 
-------------------------------------------------------------------------------
""")
        for col in matrix.columns:
            column_data = matrix[col].values.tolist()
            for i, value in enumerate(column_data):
                file.write(str("{:.13f}".format(value)))
                if (i + 1) % 5 == 0:
                    file.write('\n')
                else:
                    file.write(' ')
            file.write('\n')
    print(f"{file_name} saved...")

def gram_schma(cmo, core_index, val_index, S):
    orb = cmo.copy()
    core_orbs = orb[:, core_index]
    G = core_orbs.T @ S @ core_orbs
    for j in val_index:
        v_orbital = orb[:, j].copy()
        b = core_orbs.T @ S @ v_orbital
        a = np.linalg.solve(G, b)
        v_orbital -= core_orbs @ a
        norm = np.sqrt(v_orbital.T @ S @ v_orbital)
        if norm > 1e-10:
            v_orbital /= norm
        else:
            v_orbital[:] = 0
        orb[:, j] = v_orbital
        overlaps = core_orbs.T @ S @ orb[:, j]
    return orb

def val_core():
    core_idx, val_idx = valence_core_index()
    ORBITAL_SETTINGS = {
        "core": core_idx,
        "valence": val_idx
    }
    def get_orbital_data(orb_type):
        indices = ORBITAL_SETTINGS[orb_type].split(",")
        ao_key = 'AOPNAOs_ALPHA' if is_open else 'AOPNAOs'
        aopnao = pd.DataFrame(aoX_trans[ao_key])
        overlap = pd.DataFrame(opera_matrix['OVERLAP'])
        overlap = overlap.rename(index=lambda x: x+1, columns=lambda x: x+1)
        submatrix_df = submatrix(overlap, *indices)
        column_indices = [i - 1 for i in submatrix_df.columns]
        return column_indices, aopnao, overlap
    results = {}
    for orb_type in ORBITAL_SETTINGS:
        col_indices, aopnao, overlap = get_orbital_data(orb_type)
        results[orb_type] = {
            "columns": col_indices,
            "aopnao": aopnao,
            "overlap": overlap
        }
    ortho_pnao = gram_schma(
        results["core"]["aopnao"].values,
        results["core"]["columns"],
        results["valence"]["columns"],
        results["core"]["overlap"].values
    )
    return ortho_pnao

def solve_rhe():
    while True:
        basis_mode = input("""Enter basis to solve reduced and Full RHE equation:
1. AO Basis
2. PNAO Basis
3. Valence-Core orthogonalized PNAO basis
Enter 1, 2, or 3: """)
        try:
            basis_mode = int(basis_mode)
            if basis_mode in [1, 2, 3]:
                break
            else:
                print("Invalid input. Please enter 1, 2, or 3.")
        except ValueError:
            print("Invalid input. Please enter a number.")
    print(f"You selected option {basis_mode}")
    select_format = int(input(""" Choose the format to input selected orbitals
1. You will be able to select subset like 1,2, 6, 12 - 19                        
2: Atom1: orbitals, Atom2:orbitals. Example -  U1: 7S 6P  5F, Cl2: 3S 3P, Cl3: 3S 3P
"""))
    if select_format == 1:
        sel_orb = input("""Enter Fock and Overlap ranges e.g  1,2, 6, 12 - 19: """)
        index_arguments = sel_orb.split(",")
    elif select_format == 2:
        sel_orb = atom_orbital()
        index_arguments = sel_orb.split(",")
    if is_open:
        alpha_or_beta = get_user_input("""Get GEE solution for Alpha or Beta spin
Enter A or B for Alpha or Beta: """)
        if alpha_or_beta.upper() == "A":
            FAO = opera_matrix["FOCK_ALPHA"]
            BODM = pd.DataFrame(opera_matrix["DENSITY_ALPHA"])
        elif alpha_or_beta.upper() == "B":
            FAO = opera_matrix["FOCK_BETA"]
            BODM = pd.DataFrame(opera_matrix["DENSITY_BETA"])
        AOPNAO = aoX_trans['AOPNAOs_ALPHA']
        AONAO = aoX_trans['AONAOs_ALPHA']
    else:
        FAO = opera_matrix["FOCK"]
        AOPNAO = aoX_trans['AOPNAOs']
        AONAO = aoX_trans['AONAOs']
        BODM = pd.DataFrame(opera_matrix["DENSITY"])
    SAO = opera_matrix['OVERLAP']
    if basis_mode == 1:
        print("All calculations will be done in the AO basis")
        S = SAO
        F = FAO
    elif basis_mode == 2:
        print("All calculations will be done in the PNAO basis")
        S = AOPNAO.T @ SAO @ AOPNAO
        F = AOPNAO.T @ FAO @ AOPNAO
    elif basis_mode == 3:
        print("All calculations will be done in the Valence-Core orthogonalized PNAO basis")
        VCOPNAO = val_core()
        S = VCOPNAO.T @ SAO @ VCOPNAO
        F = VCOPNAO.T @ FAO @ VCOPNAO
    F = pd.DataFrame(F)
    S = pd.DataFrame(S)
    F.index += 1
    F.columns += 1
    S.index += 1
    S.columns += 1
    F = submatrix(F, *index_arguments)
    S = submatrix(S, *index_arguments)
    col_index = F.columns.tolist()
    sel_eigenvalues, sel_eigenvectors = eigh(F.values, S.values)
    sel_eigenvectors = pd.DataFrame(sel_eigenvectors)
    sel_eigenvectors.index = col_index
    sel_eigenvectors.columns = col_index
    print(f"""
Eigenvalues:
{sel_eigenvalues}
          """)
    print(f"""Eigenvectors: 
{sel_eigenvectors}""")
    full_eigenvectors = sub_to_full(sel_eigenvectors)
    if basis_mode == 1:
        print("Solved eigenvectors already in the AO basis")
        full_eigenvectors_ao = full_eigenvectors
    elif basis_mode == 2:
        print("Solved eigenvectors in PNAO basis. Transforming to AO basis.")
        full_eigenvectors_ao = AOPNAO @ full_eigenvectors
    elif basis_mode == 3:
        print("Solved eigenvectors in Valence-Core orthogonalized PNAO basis. Transforming to AO basis.")
        full_eigenvectors_ao = VCOPNAO @ full_eigenvectors
    while True:
        save_visual_file = input("Create and save visualization file for eigenvectors in AO basis? (y/n): ").lower()
        if save_visual_file in ['y', 'n']:
            break
        else:
            print("Invalid input. Please enter 'y' for yes or 'n' for no.")
    if save_visual_file == 'y':
        print("\nEigenvectors already in AO basis\n")
        while True:
            fname = input("Enter file name with no extension: ")
            if len(fname) != 0:
                fname += ".40"
                break
            else:
                print("Filename cannot be empty. Please enter a valid filename.")
        save_vis_file(pd.DataFrame(full_eigenvectors_ao), fname)
        print(f"File saved as {fname}")
    else:
        print("Visualization file will not be created.")
    user_input = input("Localize eigenvectors? (y/n): ")
    if user_input.lower() == 'y':
        nao_occ = np.linalg.inv(AONAO) @ BODM @ np.linalg.inv(AONAO.T)
        val_index = [i - 1 for i in col_index]
        diag_ele = np.diag(nao_occ)
        sel_occ = np.sum(diag_ele[val_index])
        print(f"\nThere are approximately {sel_occ} electrons in the selected subset of NAOs")
        loc_c = loc_select_orb_arg(full_eigenvectors_ao, SAO, FAO, ao_basis_list_dict)
        population_analy(loc_c, SAO, ao_basis_list_dict, FAO, False)
        while True:
            save_visual_file = input("Create and save visualization file for Localized orbitals in AO basis? (y/n): ").lower()
            if save_visual_file in ['y', 'n']:
                break
            else:
                print("Invalid input. Please enter 'y' for yes or 'n' for no.")
        if save_visual_file == 'y':
            print("\nEigenvectors already in AO basis\n")
            while True:
                fname = input("Enter file name with no extension: ")
                if len(fname) != 0:
                    fname += ".40"
                    break
                else:
                    print("Filename cannot be empty. Please enter a valid filename.")
            save_vis_file(pd.DataFrame(loc_c), fname)
            print(f"File saved as {fname}")
        else:
            print("Visualization file will not be created.")