

import numpy as np
import pandas as pd
import re
from collections import Counter
from sep_mat import  filename



spin_diff = {
    1: 1, 2: 0, 3: 1, 4: 0, 5: 1, 6: 2, 7: 3, 8: 2, 9: 1, 10: 0,
    11: 1, 12: 0, 13: 1, 14: 2, 15: 3, 16: 2, 17: 1, 18: 0, 19: 1, 20: 0,
    21: 1, 22: 2, 23: 3, 24: None, 25: 5, 26: 4, 27: 3, 28: 2, 29: None, 30: 0,
    31: 1, 32: 2, 33: 3, 34: 2, 35: 1, 36: 0, 37: 1, 38: 0, 39: 1, 40: 2,
    41: None, 42: None, 43: 5, 44: None, 45: None, 46: 0, 47: None, 48: 0, 49: 1, 50: 2,
    51: 3, 52: 2, 53: 1, 54: 0, 55: 1, 56: 0, 57: 1, 58: 2, 59: 3, 60: 4,
    61: 5, 62: 6, 63: 7, 64: 8, 65: 5, 66: 4, 67: 3, 68: 2, 69: 1, 70: 0,    
    71: 1, 72: 2, 73: 3, 74: 4, 75: 5, 76: 4, 77: 3, 78: None, 79: None, 80: 0,
    81: 1, 82: 2, 83: 3, 84: 2, 85: 1, 86: 0, 87: 1, 88: 0, 89: None, 90: None,
    91: None, 92: 4, 93: 5, 94: 6, 95: 7, 96: 8, 97: 5,  98: 4,  99: 3, 100: 2, 101: 1, 102: 0, 
    
    103: 1, 104: 2, 105: 3, 106: 4, 107: 5, 108: 4, 109: 3, 110: 2,  111: 1, 112: 0, 
    
    113: 1, 114: 2, 115: 3, 116: 2, 117: 1, 118: 0
}

atom_dict = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7,   "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16,   "Cl": 17, "Ar": 18, "K": 19, "Ca": 20,
    "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26,  "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40,
    "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56,  "La": 57, "Ce": 58, "Pr": 59, "Nd": 60,
    "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
    "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76,  "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
    "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96,  "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100,
    "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109,
    "Ds": 110, "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118
}

orbitals = ["1s", "2s", "2px", "2py", "2pz",
            
            "3s", "3px", "3py", "3pz", "4s", "3d1", "3d2", "3d3", "3d4", "3d5", "4px", "4py", "4pz", 
            
"5s", "4d1", "4d2", "4d3", "4d4", "4d5", "5px", "5py",  "5pz",  "6s", "4f1",  "4f2",  "4f3",  "4f4",  "4f5",  "4f6",  "4f7", 
            
"5d1",  "5d2",  "5d3",  "5d4",  "5d5", "6px", "6py", "6pz",  "7s",  "5f1", "5f2", "5f3", "5f4", "5f5", "5f6", "5f7", 

"6d1", "6d2", "6d3", "6d4", "6d5", "7px", "7py", "7pz"]



orbitalsValue = {"s": 1, "px": 1, "py": 1, "pz": 1,  "d1": 1, "d2": 1, "d3": 1, "d4": 1, "d5": 1, 
"f1": 1, "f2": 1, "f3": 1, "f4": 1, "f5": 1, "f6": 1, "f7": 1}


def electron_configuration(electrons):
    
    # if electron== 5:
    #     orbitals.replace("")
    
    
    result = []
    for orbital in orbitals:
        if electrons == 0:
            break
        value = min(1, electrons)
        result.append(f"{orbital}_{value}")
        electrons -= value
        if electrons == 0:
            break
    
    if len(result) == 0:
        return "1s_0"
    
    last_orbital = result[-1].split('_')[0]
    last_value = float(result[-1].split('_')[1])
    
    # Handle p orbitals
    if last_orbital.endswith('px'):
        average = last_value / 3
        result[-1] = f"{last_orbital}_{'%.5f' % average}"
        result.extend([f"{last_orbital[:-2]}py_{'%.5f' % average}", f"{last_orbital[:-2]}pz_{'%.5f' % average}"])
    elif last_orbital.endswith('py'):
        average = (1 + last_value) / 3
        result[-2] = f"{last_orbital[:-2]}px_{'%.5f' % average}"
        result[-1] = f"{last_orbital}_{'%.5f' % average}"
        result.append(f"{last_orbital[:-2]}pz_{'%.5f' % average}")
    
    # Handle d orbitals
    elif last_orbital.endswith('d1'):
        average = last_value / 5
        result[-1] = f"{last_orbital}_{'%.5f' % average}"
        result.extend([f"{last_orbital[:-2]}d2_{'%.5f' % average}", f"{last_orbital[:-2]}d3_{'%.5f' % average}", 
                       f"{last_orbital[:-2]}d4_{'%.5f' % average}", f"{last_orbital[:-2]}d5_{'%.5f' % average}"])
    elif last_orbital.endswith('d2'):
        average = (1 + last_value) / 5
        result[-2] = f"{last_orbital[:-2]}d1_{'%.5f' % average}"
        result[-1] = f"{last_orbital}_{'%.5f' % average}"
        result.extend([f"{last_orbital[:-2]}d3_{'%.5f' % average}", f"{last_orbital[:-2]}d4_{'%.5f' % average}", 
                       f"{last_orbital[:-2]}d5_{'%.5f' % average}"])
    elif last_orbital.endswith('d3'):
        average = (2 + last_value) / 5
        result[-3] = f"{last_orbital[:-2]}d1_{'%.5f' % average}"
        result[-2] = f"{last_orbital[:-2]}d2_{'%.5f' % average}"
        result[-1] = f"{last_orbital}_{'%.5f' % average}"
        result.extend([f"{last_orbital[:-2]}d4_{'%.5f' % average}", f"{last_orbital[:-2]}d5_{'%.5f' % average}"])
    elif last_orbital.endswith('d4'):
        average = (3 + last_value) / 5
        result[-4] = f"{last_orbital[:-2]}d1_{'%.5f' % average}"
        result[-3] = f"{last_orbital[:-2]}d2_{'%.5f' % average}"
        result[-2] = f"{last_orbital[:-2]}d3_{'%.5f' % average}"
        result[-1] = f"{last_orbital}_{'%.5f' % average}"
        result.append(f"{last_orbital[:-2]}d5_{'%.5f' % average}")
    
    # Handle f orbitals
    elif last_orbital.endswith('f1'):
        average = last_value / 7
        result[-1] = f"{last_orbital}_{'%.5f' % average}"
        result.extend([f"{last_orbital[:-2]}f2_{'%.5f' % average}", f"{last_orbital[:-2]}f3_{'%.5f' % average}", 
                       f"{last_orbital[:-2]}f4_{'%.5f' % average}", f"{last_orbital[:-2]}f5_{'%.5f' % average}", 
                       f"{last_orbital[:-2]}f6_{'%.5f' % average}", f"{last_orbital[:-2]}f7_{'%.5f' % average}"])
    elif last_orbital.endswith('f2'):
        average = (1 + last_value) / 7
        result[-2] = f"{last_orbital[:-2]}f1_{'%.5f' % average}"
        result[-1] = f"{last_orbital}_{'%.5f' % average}"
        result.extend([f"{last_orbital[:-2]}f3_{'%.5f' % average}", f"{last_orbital[:-2]}f4_{'%.5f' % average}", 
                       f"{last_orbital[:-2]}f5_{'%.5f' % average}", f"{last_orbital[:-2]}f6_{'%.5f' % average}", 
                       f"{last_orbital[:-2]}f7_{'%.5f' % average}"])
    elif last_orbital.endswith('f3'):
        average = (2 + last_value) / 7
        result[-3] = f"{last_orbital[:-2]}f1_{'%.5f' % average}"
        result[-2] = f"{last_orbital[:-2]}f2_{'%.5f' % average}"
        result[-1] = f"{last_orbital}_{'%.5f' % average}"
        result.extend([f"{last_orbital[:-2]}f4_{'%.5f' % average}", f"{last_orbital[:-2]}f5_{'%.5f' % average}", 
                       f"{last_orbital[:-2]}f6_{'%.5f' % average}", f"{last_orbital[:-2]}f7_{'%.5f' % average}"])
    elif last_orbital.endswith('f4'):
        average = (3 + last_value) / 7
        result[-4] = f"{last_orbital[:-2]}f1_{'%.5f' % average}"
        result[-3] = f"{last_orbital[:-2]}f2_{'%.5f' % average}"
        result[-2] = f"{last_orbital[:-2]}f3_{'%.5f' % average}"
        result[-1] = f"{last_orbital}_{'%.5f' % average}"
        result.extend([f"{last_orbital[:-2]}f5_{'%.5f' % average}", f"{last_orbital[:-2]}f6_{'%.5f' % average}", 
                       f"{last_orbital[:-2]}f7_{'%.5f' % average}"])
    elif last_orbital.endswith('f5'):
        average = (4 + last_value) / 7
        result[-5] = f"{last_orbital[:-2]}f1_{'%.5f' % average}"
        result[-4] = f"{last_orbital[:-2]}f2_{'%.5f' % average}"
        result[-3] = f"{last_orbital[:-2]}f3_{'%.5f' % average}"
        result[-2] = f"{last_orbital[:-2]}f4_{'%.5f' % average}"
        result[-1] = f"{last_orbital}_{'%.5f' % average}"
        result.extend([f"{last_orbital[:-2]}f6_{'%.5f' % average}", f"{last_orbital[:-2]}f7_{'%.5f' % average}"])
    elif last_orbital.endswith('f6'):
        average = (5 + last_value) / 7
        result[-6] = f"{last_orbital[:-2]}f1_{'%.5f' % average}"
        result[-5] = f"{last_orbital[:-2]}f2_{'%.5f' % average}"
        result[-4] = f"{last_orbital[:-2]}f3_{'%.5f' % average}"
        result[-3] = f"{last_orbital[:-2]}f4_{'%.5f' % average}"
        result[-2] = f"{last_orbital[:-2]}f5_{'%.5f' % average}"
        result[-1] = f"{last_orbital}_{'%.5f' % average}"
        result.append(f"{last_orbital[:-2]}f7_{'%.5f' % average}")
    
    
    return ' '.join(result)



def print_config(lst):
    for atomic_number in lst:
    
        # Calculate the number of alpha and beta electrons
        diff = spin_diff.get(atomic_number, 0)
        alpha_electrons = (atomic_number + diff) // 2
        beta_electrons = atomic_number - alpha_electrons

def read_file46(file_path):
    try:
        with open(file_path, 'r') as file:
            content = ""
            start_reading = False
            
            for line in file:
                if "NAO" in line:
                    start_reading = True
                    continue
                if start_reading and "NHO" not in line:
                    content += line
                if "NHO" in line:
                    break

            if content:
                # Extract strings between parentheses (orbitals)
                orbital_list = re.findall(r'\((.*?)\)', content)
                orbital_list = [orbital.strip() for orbital in orbital_list]
                
                # Extract atom labels with their numbers
                atom_labels = re.findall(r'([A-Za-z]+)\s*(\d+)\s*\(', content)
                atom_labels = [f"{atom}{num}" for atom, num in atom_labels]
                
                # Create a dictionary counting occurrences of each atom label
                atom_count = dict(Counter(atom_labels))
                
                # Create a dictionary with atom labels as keys and corresponding orbitals as values
                atom_orbital_dict = {}
                start_index = 0
                for atom, count in atom_count.items():
                    end_index = start_index + count
                    atom_orbital_dict[atom] = orbital_list[start_index:end_index]
                    start_index = end_index
                
                
                return orbital_list, atom_labels, atom_count, atom_orbital_dict
            else:
                print("No content found between 'NAO' and 'NHO'.")
                return [], [], {}, {}
                
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
    except IOError:
        print(f"Error: Unable to read file '{file_path}'.")
    
    return [], [], {}, {}



file_path = filename + '.46'
orb_lst, at_label, at_count, atom_orb_dict = read_file46(file_path)

at_lst = list(at_count.keys())


# Remove numbers from each atom label
at_symb = [re.sub(r'\d', '', atom) for atom in at_lst]

at_lst = at_symb

unique_at = list(dict.fromkeys(at_lst))


unique_at_orb_dict = {}
for key, value in atom_orb_dict.items():
    # Remove numbers from the key
    new_key = ''.join([i for i in key if not i.isdigit()])
    
    # If the key already exists, we don't need to add it again
    if new_key not in unique_at_orb_dict:
        unique_at_orb_dict[new_key] = value




at_prom_main = {}
atom_config_dict = {}
for atom, orbs in zip(unique_at, unique_at_orb_dict.values()):
    # Create DataFrames
    atom_df = pd.DataFrame(columns=orbs, index=[atom])
    atom_alpha = pd.DataFrame(columns=orbs, index=[atom])
    atom_beta = pd.DataFrame(columns=orbs, index=[atom])
    
    # Fill DataFrames with zeros
    atom_df.loc[atom] = 0
    atom_alpha.loc[atom] = 0
    atom_beta.loc[atom] = 0
    
        
    
    atomic_number = int(atom_dict.get(atom))
    
    if atomic_number == 24: #15, 9
        alpha_electrons = 15
        # alpha_config = """1s_1 2s_1 2px_1 2py_1 2pz_1 3s_1 3px_1 3py_1 3pz_1 4s_1 3d1_1.0000 3d2_1.0000 3d3_1.0000 3d4_1.0000 3d5_1.0000"""
        beta_electrons = 9
        beta_config =  """1s_1 2s_1 2px_1 2py_1 2pz_1 3s_1 3px_1 3py_1 3pz_1 """
        alpha_config = beta_config + """ 4s_1 3d1_1.0000 3d2_1.0000 3d3_1.0000 3d4_1.0000 3d5_1.0000"""

        
    elif atomic_number == 29: #15, 14 
        alpha_electrons = 15
        alpha_config = """1s_1 2s_1 2px_1 2py_1 2pz_1 3s_1 3px_1 3py_1 3pz_1 3d1_1.0000 3d2_1.0000 3d3_1.0000 3d4_1.0000 3d5_1.0000 4s_1 """
        beta_electrons = 14            
        # beta_config =  """1s_1 2s_1 2px_1 2py_1 2pz_1 3s_1 3px_1 3py_1 3pz_1 3d1_1.0000 3d2_1.0000 3d3_1.0000 3d4_1.0000 3d5_1.0000 4s_0"""
        beta_config = alpha_config.replace("4s_1", "4s_0")


    elif atomic_number == 41:  ###15, 14  = [Kr] 4d5 5s0 ,, Kr/2 = 18 
        alpha_electrons = 23
        # alpha_config = electron_configuration(18) + " 4d1_1 4d2_1 4d3_1 4d4_1 4d5_1 5s_0" #Verified with Orca
        beta_electrons = 18            
        beta_config =  electron_configuration(18)        
        alpha_config = beta_config + " 4d1_1 4d2_1 4d3_1 4d4_1 4d5_1 5s_0" #Verified with Orca

          
    elif atomic_number == 42:  ###15, 14  = [Kr] 4d5 5s1 ,, Kr/2 = 18 
        alpha_electrons = 24
        # alpha_config = electron_configuration(18) + " 4d1_1 4d2_1 4d3_1 4d4_1 4d5_1 5s_1" #Verified with Orca
        beta_electrons = 18            
        beta_config =  electron_configuration(18)
        alpha_config = beta_config + " 4d1_1 4d2_1 4d3_1 4d4_1 4d5_1 5s_1" #Verified with Orca

                         
    elif atomic_number == 44:  ###15, 14  = [Kr] 4d7 5s1         
        alpha_electrons = 24
        alpha_config = electron_configuration(18) + " 4d1_1 4d2_1 4d3_1 4d4_1 4d5_1 5s_1" #Verified with Orca
        beta_electrons = 20            
        beta_config =  electron_configuration(18) + " 4d1_0.40000 4d2_0.40000 4d3_0.40000 4d4_0.40000 4d5_0.40000 5s_0"
          
    elif atomic_number == 45:  ###15, 14  = [Kr] 4d9 5s0         
        alpha_electrons = 23
        alpha_config = electron_configuration(18) + " 4d1_1 4d2_1 4d3_1 4d4_1 4d5_1 5s_1" #Verified with Orca WITH DEF2QZVP
        beta_electrons = 22            
        beta_config =  electron_configuration(18) + " 4d1_0.80000 4d2_0.80000 4d3_0.80000 4d4_0.80000 4d5_0.80000 5s_0"
    
          
    elif atomic_number == 47:  ###15, 14  = [Kr] 4d10 5s1         
        alpha_electrons = 24
        alpha_config = electron_configuration(18) + " 4d1_1 4d2_1 4d3_1 4d4_1 4d5_1 5s_1" #Verified with Orca WITH DEF2QZVP
        beta_electrons = 23            
        beta_config = alpha_config.replace("5s_1", "5s_0")
    
                    
    elif atomic_number == 78:  ###15, 14  = [Xe] 4f14 5d9 6s1        
        alpha_electrons = 40
        alpha_config = electron_configuration(27) + " 6s_1 4f1_1 4f2_1 4f3_1 4f4_1 4f5_1 4f6_1 4f7_1 5d1_1 5d2_1 5d3_1 5d4_1 5d5_1"
        beta_electrons = 38            
        beta_config =  electron_configuration(27) + " 6s_0 4f1_1 4f2_1 4f3_1 4f4_1 4f5_1 4f6_1 4f7_1 5d1_0.8 5d2_0.8 5d3_0.8 5d4_0.8 5d5_0.8"

          
    elif atomic_number == 79:  ###15, 14  = [Xe] 4f14 5d10 6s1            
        alpha_electrons = 40
        alpha_config = electron_configuration(27) + " 6s_1 4f1_1 4f2_1 4f3_1 4f4_1 4f5_1 4f6_1 4f7_1 5d1_1 5d2_1 5d3_1 5d4_1 5d5_1"
        beta_electrons = 39            
        beta_config = alpha_config.replace("6s_1", "6s_0")
    

    else:
        # Calculate the number of alpha and beta electrons
        diff = spin_diff.get(atomic_number)
        alpha_electrons = (atomic_number + diff) // 2
        beta_electrons = atomic_number - alpha_electrons
        
        # Get configurations
        alpha_config = electron_configuration(alpha_electrons)
        beta_config = electron_configuration(beta_electrons)
        
        
    
   
    # Populate atom_alpha DataFrame
    for orb_config in alpha_config.split():
        orb, occupation = orb_config.split('_')
        if orb in atom_alpha.columns:
            atom_alpha.at[atom, orb] = float(occupation)
    
    # Populate atom_beta DataFrame
    for orb_config in beta_config.split():
        orb, occupation = orb_config.split('_')
        if orb in atom_beta.columns:
            atom_beta.at[atom, orb] = float(occupation)

    
    N = len(orbs)
    alpha_array = np.zeros((N, N))
    beta_array = np.zeros((N, N))
    
    for i, orb in enumerate(orbs):
        alpha_array[i, i] = atom_alpha.at[atom, orb]
        beta_array[i, i] = atom_beta.at[atom, orb]
    
    atom_config_dict[f'{atom}_alpha'] = alpha_array
    atom_config_dict[f'{atom}_beta'] = beta_array



def get_full_prom_den(atom_basis_ind_tup, at_prom_main):
    # Determine the total number of basis functions
    total_basis_functions = max(max(indices) for indices in atom_basis_ind_tup.values()) + 1

    alpha_matrix = np.zeros((total_basis_functions, total_basis_functions))
    beta_matrix = np.zeros((total_basis_functions, total_basis_functions))

    for atom_key, indices in atom_basis_ind_tup.items():
        atom_type = atom_key.split('_')[0]  # Extract atom type (C, O, H, etc.)
        
        # Get the corresponding promolecule density matrices
        alpha_block = at_prom_main[f"{atom_type}_alpha"]
        beta_block = at_prom_main[f"{atom_type}_beta"]

        # Ensure the block size matches the number of indices
        block_size = len(indices)
        if alpha_block.shape[0] > block_size:
            alpha_block = alpha_block[:block_size, :block_size]
            beta_block = beta_block[:block_size, :block_size]
        elif alpha_block.shape[0] < block_size:
            pad_size = block_size - alpha_block.shape[0]
            alpha_block = np.pad(alpha_block, ((0, pad_size), (0, pad_size)))
            beta_block = np.pad(beta_block, ((0, pad_size), (0, pad_size)))

        # Place the blocks in the appropriate positions
        for i, row_idx in enumerate(indices):
            for j, col_idx in enumerate(indices):
                alpha_matrix[row_idx, col_idx] = alpha_block[i, j]
                beta_matrix[row_idx, col_idx] = beta_block[i, j]

    return alpha_matrix, beta_matrix



def create_new_dict(atom_orb_dict):
    atom_basis_ind_tup = {}
    start = 0
    for atom, orbitals in atom_orb_dict.items():
        atom_sym = ''.join(char for char in atom if char.isalpha())
        atom_ind = ''.join(char for char in atom if char.isdigit())
        new_key = atom_sym + '_' + atom_ind
        end = start + len(orbitals)
        atom_basis_ind_tup[new_key] = tuple(range(start, end))
        start = end
    return atom_basis_ind_tup

atom_basis_ind_tup= create_new_dict(atom_orb_dict)

alp_promol_den_mat, bet_promol_den_mat = get_full_prom_den(atom_basis_ind_tup, atom_config_dict)


