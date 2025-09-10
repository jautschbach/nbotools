import numpy as np
import re
import time
import extract_nbo_files
from basis_info import pre_orb_from_orb

np.set_printoptions(suppress=True, precision=8, floatmode='fixed')

start_time = time.time()


matrix_dict, filename, NBAS = extract_nbo_files.run_processing()
matrix_keys = matrix_dict.keys()
# print(matrix_keys)


# def is_unit_matrix(matrix, tolerance=1e-6):
    # return np.allclose(matrix, np.eye(matrix.shape[0]), atol=tolerance)


# def is_unit_matrix(matrix, tolerance=1e-7):
#     # Check if the matrix is square
#     size = len(matrix)
#     if any(len(row) != size for row in matrix):
#         return False

#     # Check unit matrix properties with tolerance
#     for i in range(size):
#         for j in range(size):
#             value = matrix[i][j]
#             if i == j:
#                 # Check diagonal elements are close to 1
#                 if abs(value - 1) > tolerance:
#                     return False
#             else:
#                 # Check off-diagonal elements are close to 0
#                 if abs(value) > tolerance:
#                     return False
                    
#     return True

def is_unit_matrix(matrix, tolerance=1e-7):
    # Check if the matrix is square
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        return False

    # Check diagonal and off-diagonal elements
    for i in range(size):
        # Check diagonal element
        if abs(matrix[i][i] - 1) > tolerance:
            return False
        
        # Check off-diagonal elements in the current row
        for j in range(size):
            if i != j and abs(matrix[i][j]) > tolerance:
                return False

    return True






def get_user_choice():
    while True:
        usr_inp = input("""AOPNAOs matrix options:
1: Use NBO generated AOPNAOs
2: Recalculate AOPNAOs from AONAOs. SPNAOs will be recalulated
and may slightly be different from NBO generated AOPNAOs.
Enter your choice (1 or 2): """)
        
        if usr_inp in ['1', '2']:
            return int(usr_inp)
        else:
            print("Invalid input. Please enter 1 or 2.")
            
def test_orthogonality(mat_dict, ove_mat):
    orthogonal_results = {}
    pre_orthogonal_results = {}

    for key, matrix in mat_dict.items():
        
        
        result = np.transpose(matrix) @ ove_mat @ matrix
        is_orthogonal = is_unit_matrix(result)
        
        if   (key[0] == "A" and key[2] == "P") or (key[0] == "P" and key[4] == "P") :
            pre_orthogonal_results[key] = is_orthogonal
        else:
            orthogonal_results[key] = is_orthogonal

    return orthogonal_results, pre_orthogonal_results

def print_results(orthogonal_results, pre_orthogonal_results):
    print("Orthogonal basis results:")
    for key, is_orthogonal in orthogonal_results.items():
        print(f"{is_orthogonal} for {key}.T @ S @ {key} = 1")
    
    print("\nPre-orthogonal basis results:")
    for key, is_orthogonal in pre_orthogonal_results.items():
        print(f"{is_orthogonal} for {key}.T @ S @ {key} = 1")

    print("\nSummary:")
    print(f"All orthogonal bases are unit matrices: {all(orthogonal_results.values())}")
    print(f"All pre-orthogonal bases are non-unit matrices: {not any(pre_orthogonal_results.values())}")
    



def recalculate_aopnaos(matrix_dict):
    SAO = matrix_dict.get("OVERLAP")
    AONAO = matrix_dict.get("AONAOs")
    if SAO is None or AONAO is None:
        print("Error: Required matrices not found in matrix_dict.")
        return None
    
    AOPNAO_from_AONAO = pre_orb_from_orb(np.array(AONAO), np.array(SAO))
    return AOPNAO_from_AONAO

def update_aopnaos(matrix_dict):
    choice = get_user_choice()
    
    if choice == 2:
        AOPNAO = recalculate_aopnaos(matrix_dict)
        if AOPNAO is not None:
            matrix_dict["AOPNAOs"] = AOPNAO
            print("AOPNAOs recalculated and updated in matrix_dict.")
    else:
        print("Using NBO generated AOPNAOs.")
    
    return matrix_dict

matrix_dict = update_aopnaos(matrix_dict)
variables = matrix_dict

# print(matrix_keys)
is_open = [key for key in matrix_keys if "ALPHA" and "BETA"  in key ]


def process_all_shell(matrix_dict, is_open):
    def categorize_matrices(matrix_dict):
        
        if is_open:
            aoX_trans = {k: v for k, v in matrix_dict.items() if k.startswith("AO") and "AODM" not in k }

        else:            
            aoX_trans = {k: v for k, v in matrix_dict.items() if k.startswith("AO") and "AODM" not in k and "_" not in k}
        opera_matrix = {k: v for k, v in matrix_dict.items() if not k.startswith("AO") or "AODM" in k}
        return aoX_trans, opera_matrix

    def calculate_energies(matrices, F, FA=None, FB=None):
        if is_open:
            return {
                f'{key.replace("_", " ")} energy/Ha': 
                    np.einsum('ij,ji->i', np.dot(matrix.T, FA), matrix).tolist() if 'ALPHA' in key
                    else np.einsum('ij,ji->i', np.dot(matrix.T, FB), matrix).tolist() if 'BETA' in key
                    else []
                for key, matrix in matrices.items()
            }
        else:
            return {f'{key} Energies': np.einsum('ij,ji->i', np.dot(matrix.T, F), matrix).tolist()
                    for key, matrix in matrices.items()}

    def transform_to_pnao(matrices, U):
        return {f'PN{key}': U @ matrix for key, matrix in matrices.items()}

    aoX_trans, opera_matrix = categorize_matrices(matrix_dict)
    
    SAO = matrix_dict.get("OVERLAP")
    AOPNAO = matrix_dict.get("AOPNAOs")
    
    if AOPNAO is None or SAO is None:
        print("Warning: AOPNAO or SAO matrix not found in matrix_dict")
        return None
    
    pS = AOPNAO.T @ SAO @ AOPNAO
    U = np.linalg.inv(AOPNAO)
    pnaoX_trans = transform_to_pnao(aoX_trans, U)

    if is_open:
        fock_alpha = matrix_dict.get('FOCK_ALPHA')
        fock_beta = matrix_dict.get('FOCK_BETA')
        
        if fock_alpha is None or fock_beta is None:
            print("Warning: FOCK_ALPHA or FOCK_BETA matrix not found in matrix_dict")
            return None
        
        pFa = AOPNAO.T @ fock_alpha @ AOPNAO
        pFb = AOPNAO.T @ fock_beta @ AOPNAO
        
        # Remove spin independent transformations and modify the dictionary
        aonao = aoX_trans.get("AONAOs")
        if aonao is not None:
            aonao = np.array(aonao)
            aoX_trans['AONAOs_ALPHA'] = aonao
            aoX_trans['AONAOs_BETA'] = aonao
            aoX_trans.pop('AONAOs')
        else:
            print("Warning: AONAOs not found in aoX_trans.")
        
        aopnao = aoX_trans.get("AOPNAOs")
        if aopnao is not None:
            aopnao = np.array(aopnao)
            aoX_trans['AOPNAOs_ALPHA'] = aopnao
            aoX_trans['AOPNAOs_BETA'] = aopnao
            aoX_trans.pop('AOPNAOs')
        else:
            print("Warning: AOPNAOs not found in aoX_trans.")
        
        # Calculate energies
        aoX_energies_dict = calculate_energies(aoX_trans, None, fock_alpha, fock_beta)
        
        # Process PNAO data
        pnaonao = pnaoX_trans.get("PNAONAOs")
        if pnaonao is not None:
            pnaonao = np.array(pnaonao)
            pnaoX_trans['PNAONAOs_ALPHA'] = pnaonao
            pnaoX_trans['PNAONAOs_BETA'] = pnaonao
            pnaoX_trans.pop("PNAONAOs")
        else:
            print("Warning: PNAONAOs not found in pnaoX_trans.")
        
        pnaopnao = pnaoX_trans.get("PNAOPNAOs")
        if pnaopnao is not None:
            pnaopnao = np.array(pnaopnao)
            pnaoX_trans['PNAOPNAOs_ALPHA'] = pnaopnao
            pnaoX_trans['PNAOPNAOs_BETA'] = pnaopnao
            pnaoX_trans.pop("PNAOPNAOs")
        else:
            print("Warning: PNAOPNAOs not found in pnaoX_trans.")   
            
        pnaoX_energies_dict = calculate_energies(pnaoX_trans, None, pFa, pFb)

    
   

    else:
        FAO = matrix_dict.get("FOCK")
        
        if FAO is None:
            print("Warning: FOCK matrix not found in matrix_dict")
            return None
        
        pF = AOPNAO.T @ FAO @ AOPNAO
        
        aoX_energies_dict = calculate_energies(aoX_trans, FAO)
        pnaoX_energies_dict = calculate_energies(pnaoX_trans, pF)
                        
    
    opera_matrix['PNAO OVERLAP'] = pS
    
    if is_open:
        opera_matrix['PNAO Alpha FOCK'] = pFa
        opera_matrix['PNAO Beta FOCK'] = pFb

    else:
        opera_matrix['PNAO FOCK'] = pF
        
    
    
    
    """TESTS STARTS HERE"""    
    

    def ortho_pass_fail(orthogonal_results: dict, pre_orthogonal_results: dict) -> None:
        def check_and_print(results: dict, condition: bool, message: str) -> None:
            failed_keys = [key for key, value in results.items() if value is not condition]
            if failed_keys:
                print(f"{message} {', '.join(failed_keys)}")
    
        check_and_print(orthogonal_results, True, "Orthogonality test failed for")
        check_and_print(pre_orthogonal_results, False, "Pre-orthogonality test failed for")
    
    def test_energy_differences(aoX_energies_dict, pnaoX_energies_dict):
        non_zero_diffs = []
        for (key, a), (pkey, b) in zip(aoX_energies_dict.items(), pnaoX_energies_dict.items()):
            if not np.allclose(a, b, atol=1e-9):
                non_zero_diffs.append(f"{key} - {pkey}")
        
        if non_zero_diffs:
            print("Some energy differences are not zero:")
            print("\n".join(non_zero_diffs))
    
    def run_tests(aoX_trans, SAO, pnaoX_trans, pS, aoX_energies_dict, pnaoX_energies_dict):
        print("\nRunning orthogonality tests...")
        AO_ortho_results, AO_pre_ortho_results = test_orthogonality(aoX_trans, SAO)
        ortho_pass_fail(AO_ortho_results, AO_pre_ortho_results)
    
        PNAO_ortho_results, PNAO_pre_ortho_results = test_orthogonality(pnaoX_trans, pS)
        ortho_pass_fail(PNAO_ortho_results, PNAO_pre_ortho_results)
    
        print("\nChecking energy differences...")
        test_energy_differences(aoX_energies_dict, pnaoX_energies_dict)
    
    #Run All test
    run_tests(aoX_trans, SAO, pnaoX_trans, pS, aoX_energies_dict, pnaoX_energies_dict)
    
    """TESTS ENDS HERE HERE"""  
  

    return {
        'aoX_trans': aoX_trans,
        'opera_matrix': opera_matrix,
        'aoX_energies_dict': aoX_energies_dict,
        'pnaoX_trans': pnaoX_trans,
        'pnaoX_energies_dict': pnaoX_energies_dict,
        
    }

# Usage
result = process_all_shell(matrix_dict, is_open)
if result:
    # print("Processing completed successfully")
    aoX_trans = result['aoX_trans']
    pnaoX_trans = result['pnaoX_trans']
    opera_matrix = result['opera_matrix']    
    
    aoX_energies_dict = result['aoX_energies_dict']
    pnaoX_energies_dict = result['pnaoX_energies_dict']
    

else:
    print("Processing failed")


# Remove original PNAO overlap matrix
# The final PNAO OVERLAP in OPERA_MATRIX be the same as NBO generated SPNAO if user does not modify the AOPNAO
# Might slightly change if the AOPNAO is modified
# This is not an error. Make sure it is always removed
# opera_matrix.pop('PNAO_overlap_matrix', None)

# print(opera_matrix.keys())
# print()
# print(aoX_trans.keys())

total_processing_time = time.time() - start_time
# print(f"Total processing time: {total_processing_time:.4f} seconds")



def list_to_blo_bi(list_ele):
    unique_items = set(list_ele)
    result_list = [{'bflo': list_ele.index(item), 'bfhi': len(list_ele) - list_ele[::-1].index(item) - 1} for item in unique_items]
    return sorted(result_list, key=lambda x: x['bflo'])


def get_bf_low_high(fn, NBAS):
    
    
       
    with open(fn, 'r') as file:
        lines = file.readlines()

    NBAS = str(NBAS)
    pattern = r"[A-Za-z][\w/$@(]*\)"
    
    basis_dictionary = {}
    current_header = None

    for line in lines:
        line = line.strip()
        
        if any(line.startswith(header) for header in ['AO', 'NAO', 'NBO', 'NHO']):
            current_header = line.replace(NBAS, "").replace(" ", "").replace("ALPHA", " ALPHA").replace("BETA", " BETA")
            basis_dictionary[current_header] = []
        else:
            result = re.findall(pattern, line.replace(" ", ""))
            basis_dictionary[current_header].extend(result)

    basis_dictionary = {key: value for key, value in basis_dictionary.items() if "NBO" not in key}
    atom_basis_dict = {key: [value.split('(')[0] for value in value_list] for key, value_list in basis_dictionary.items()}
    
    is_alpha_beta = any("BETA" and "ALPHA" in key for key in atom_basis_dict.keys())
            
    ao_basis = atom_basis_dict['AO ALPHA'] if is_alpha_beta else atom_basis_dict['AO']

 
    ao_basis_list_dict = list_to_blo_bi(ao_basis)
    
    return ao_basis, ao_basis_list_dict


ao_basis, ao_basis_list_dict = get_bf_low_high(filename + '.46', NBAS)
# print(ao_basis)
# print(ao_basis_list_dict)




