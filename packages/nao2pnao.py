import re
from itertools import groupby
import numpy as np
import timeit
import os
from multiprocessing import Pool
from functools import partial
import numpy as np
import pandas as pd
from scipy.constants import physical_constants
import sys

bohr =  physical_constants['Bohr radius'][0]*1e10

def read_file47(filename):
    try:
        with open(filename, 'r') as file:
            content = file.read()
            return content
    except FileNotFoundError:
        print(f"File '{filename}' not found.")
    except IOError:
        print(f"Error reading file '{filename}'.")


# Check if the correct number of arguments are provided
if len(sys.argv) != 4:
    print("""Usage: python3 nao2pnao.py FILE.47 FILE.33 FILE.32a
FILE.32a is used so we don't overwrite the original NBO provided FILE.32
""")
    sys.exit(1)

# Assign the file names to variables
f47_file = sys.argv[1]
orbital_file = sys.argv[2]
output_file = sys.argv[3]
        
# f47_file = 'examples/close-shell/h2o-b3lyp-def2svp.47'# input("Enter FILE47 to get basis info: ")
# f47_file = "h2o-b3lyp-def2svp.47"# input("Enter FILE47 to get basis info: ")
filename = f47_file
file_content = read_file47(filename)
origin = [0, 0 , 0]


def system_info(content):
    pattern = r'\s+(\d+)\s+(\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)'
    matches = re.findall(pattern, content)
    
    atom_data = []
    for match in matches:
        atomic_number = int(match[0])
        charge = int(match[1])
        x = float(match[2])
        y = float(match[3])
        z = float(match[4])
        atom_data.append((atomic_number, charge, x, y, z))
    
    variable_names = ['EXP', 'CS', 'CP', 'CD', 'CF', 'CG', 'CH', 'CI', 'CJ']  # Add more variable names here if needed

    patterns = {}
    data = {}

    for variable in variable_names:
        pattern = rf'{variable}\s+=\s+((?:[-+]?\d+\.\d+E[+-]?\d+\s+)+)'
        patterns[variable] = pattern
        data[variable] = []


    for variable, pattern in patterns.items():
        matches = re.findall(pattern, content)
        for match in matches:
            values = match.split()
            data[variable].extend([float(value) for value in values])
            
    # print(data.get("CS"))
    
    for key, value in data.items():
        globals()[f'{key}'] = value
    

    variable_names = ['CENTER', 'LABEL', 'NSHELL', 'NEXP', 'NCOMP', 'NPRIM', 'NPTR']  # Add more variable names here if needed

    patterns = {}
    data = {}

    for variable in variable_names:
        pattern = rf'{variable}\s+=\s+([\d\s]+)'
        patterns[variable] = pattern
        data[variable] = []


    for variable, pattern in patterns.items():
        matches = re.findall(pattern, content)
        for match in matches:
            values = match.split()
            data[variable].extend([int(value) for value in values])
            
    # print(data.get("CS"))
    
    for key, value in data.items():
        globals()[f'{key}'] = value
        

    orb_val = []
    orb_type = []
    
    #Define the mapping of aotype to orb_type and orb_val
    orb_mapping = {
        1: ('s', 's'),
        
        101: ('px', 'px'),
        102: ('py', 'py'),
        103: ('pz', 'pz'),
        
        251: ('d_xy', 'ds2'),
        252: ('d_xz', 'ds1'),
        253: ('dd_yz', 'dc1'),
        254: ('d_x2-y2', 'dc2'),
        255: ('d_z2', 'd0'),
        
        351: ('fz(5z2-3r2)', 'f0'),
        352: ('fx(5z2-r2)', 'fc1'),
        353: ('fy(5z2-r2)', 'fs1'),
        354: ('fz(x2-y2)', 'fc2'),
        355: ('fxyz', 'fs2'),
        356: ('fx(x2-3y2)', 'fc3'), 
        357: ('f(3x2-y2)', 'fs3'),
        
        451: ('g0', 'g0'),
        452: ('gc1', 'gc1'),
        453: ('gs1', 'gs1'),
        454: ('gc2', 'gc2'),
        455: ('gs2', 'gs2'),
        456: ('gc3', 'gc3'),
        457: ('gs3', 'gs3'),
        458: ('gc4', 'gc4'),
        459: ('gs4', 'gs4'),
        
        551: ('h0', 'h0'),
        552: ('hc1', 'hc1'),
        553: ('hs1', 'hs1'),
        554: ('hc2', 'hc2'),
        555: ('hs2', 'hs2'),
        556: ('hc3', 'hc3'),
        557: ('hs3', 'hs3'),
        558: ('hc4', 'hc4'),
        559: ('hs4', 'hs4'),
        560: ('hc5', 'hc5'),
        561: ('hs5', 'hs5'),
        
        651: ('i0', 'i0'),
        652: ('ic1', 'ic1'),
        653: ('is1', 'is1'),
        654: ('ic2', 'ic2'),
        655: ('is2', 'is2'),
        656: ('ic3', 'ic3'),
        657: ('is3', 'is3'),
        658: ('ic4', 'ic4'),
        659: ('is4', 'is4'),
        660: ('ic5', 'ic5'),
        661: ('is5', 'is5'),
        662: ('ic6', 'ic6'),
        663: ('is6', 'is6'),
        
        751: ('j0', 'j0'),
        752: ('jc1', 'jc1'),
        753: ('js1', 'js1'),
        754: ('jc2', 'jc2'),
        755: ('js2', 'js2'),
        756: ('jc3', 'jc3'),
        757: ('js3', 'js3'),
        758: ('jc4', 'jc4'),
        759: ('js4', 'js4'),
        760: ('jc5', 'jc5'),
        761: ('js5', 'js5'),
        762: ('jc6', 'jc6'),
        763: ('js6', 'js6'),
        764: ('jc7', 'jc7'),
        765: ('js7', 'js7')
    }
    
    
    for aotype in LABEL:
        # Get the mapping from the dictionary, default to ('other', 'other') if not found
        type_val_pair = orb_mapping.get(aotype, ('other', 'other'))
        orb_type.append(type_val_pair[0])
        orb_val.append(type_val_pair[1])
       
     
    shell_num = []
    count = 0
    orb_count = 0  # Variable to keep track of the number of orb_type elements in the current shell

    for i in range(len(orb_type)):
        if i > 0 and orb_type[i][0] in ['p', 'd', 'f', 'g', 'h', 'i', 'j'] and orb_type[i-1][0] == orb_type[i][0]:
            if orb_count < 4 and orb_type[i][0] == 'p':
                shell_num.append(shell_num[-1])
            elif orb_count < 6 and orb_type[i][0] == 'd':
                shell_num.append(shell_num[-1])
            elif orb_count < 8 and orb_type[i][0] == 'f':
                shell_num.append(shell_num[-1])
            elif orb_count < 10 and orb_type[i][0] == 'g':
                shell_num.append(shell_num[-1])
            elif orb_count < 12 and orb_type[i][0] == 'h':
                shell_num.append(shell_num[-1])
            elif orb_count < 14 and orb_type[i][0] == 'i':
                shell_num.append(shell_num[-1])
            elif orb_count < 16 and orb_type[i][0] == 'j':
                shell_num.append(shell_num[-1])
            
            else:
                count += 1
                shell_num.append(count)
                orb_count = 1
            
        else:
            count += 1
            shell_num.append(count)
            orb_count = 1

        orb_count += 1
              
    shell_num_list = [list(v) for k, v in groupby(shell_num)]
        
    NPTR_new = []
    NPRIM_new = []
    
    for shell, ptr, prim in zip(shell_num_list, NPTR, NPRIM):
        NPTR_new.append([ptr] * len(shell))
        NPRIM_new.append([prim] * len(shell))
        pass
    
    PTR = [item for sublist in NPTR_new for item in sublist]
    PRIM = [item for sublist in NPRIM_new for item in sublist]
    
    
    # bas_info_dict = [{}]
    bas_info_dict = []
    N = [i for i in range(1, len(CENTER) + 1)]
   
    for i, PTR, PRIM  in zip (range(len(N)), PTR, PRIM):
        # print(i, PTR, PRIM)
        info = {
            "N" :N[i],
            "CENTER": CENTER[i],
            "LABEL":LABEL[i],
            "shell_num": shell_num[i],
            "type":orb_type[i],
            "orb_val":orb_val[i],
            "exps":EXP[PTR-1:PTR-1+PRIM]
            }
        
        coeffs = []
        coord = []
        m = []
        l = []
        n = []
        # r_square = [] # The value of r2 = x2 + y2 + z2
        
        if len(CS) > 0:  # and all(abs(elem) > 0 for elem in CS):
            coeffs.extend([elem for elem in CS[PTR - 1:PTR - 1 + PRIM] if elem != 0.0])
        if len(CP) > 0:  # and all(abs(elem) > 0 for elem in CP):
            coeffs.extend([elem for elem in CP[PTR - 1:PTR - 1 + PRIM] if elem != 0.0])
        if len(CD) > 0:  # and all(elem == 0 for elem in CD):
           coeffs.extend([elem for elem in CD[PTR - 1:PTR - 1 + PRIM] if elem != 0.0])
        if len(CF) > 0:  # and all(elem == 0 for elem in CF):
            coeffs.extend([elem for elem in CF[PTR - 1:PTR - 1 + PRIM] if elem != 0.0])
        if len(CG) > 0:  # and all(elem == 0 for elem in CG):
            coeffs.extend([elem for elem in CG[PTR - 1:PTR - 1 + PRIM] if elem != 0.0])
        if len(CH) > 0:  # and all(elem == 0 for elem in CH):
            coeffs.extend([elem for elem in CH[PTR - 1:PTR - 1 + PRIM] if elem != 0.0])
        if len(CI) > 0:  # and all(elem == 0 for elem in CI):
            coeffs.extend([elem for elem in CI[PTR - 1:PTR - 1 + PRIM] if elem != 0.0])
        if len(CJ) > 0:  # and all(elem == 0 for elem in CJ):
            coeffs.extend([elem for elem in CJ[PTR - 1:PTR - 1 + PRIM] if elem != 0.0])
        
        info["coeffs"] = coeffs
        
        coord.extend(atom_data[CENTER[i] -1][2:5])
       
        bas_info_dict.append(info)
            
    return bas_info_dict, atom_data


basis_info_dict , atom_data = system_info(file_content)


data =  basis_info_dict
atom_info = atom_data
atom_data = {
    'coordinates': [tuple(coord[2:]) for coord in atom_data],
    'atom_charge': [tuple(coord[0:2]) for coord in atom_data]
}

# for info in basis_info_dict:
#         for key, value in info.items():
#             print(f"{key}: {value}")
#         print("---------------------------")

coordinates = (atom_data.get("coordinates"))
atom_charge = np.array(atom_data.get("atom_charge"))

NBAS = len(data)


# Assuming basis_info_dict is your list of dictionaries
N = len(basis_info_dict)

# Create an N*N numpy zero array
zero_array = np.zeros((N, N))

# Modify the zero_array
for col, info in enumerate(basis_info_dict):
    current_orb_val = info['orb_val'][0]
    for row, check_info in enumerate(basis_info_dict):
        if check_info['orb_val'][0] == current_orb_val:
            zero_array[row, col] = 1

# Print the modified array
# print(zero_array)

# Print the modified array
# print(zero_array)
# print(pd.DataFrame(zero_array))


###############################################################################
def read_and_extract_mat(filename, keyword, length):
    with open(filename, 'r') as file:
        words = file.read().split()
    try:
        index = words.index(keyword)
        numbers = [float(word) for word in words[index+1:index+length+1]]
    except ValueError:
        print(f"The word '{keyword}' was not found in the file.")
        numbers = []
    return numbers

def create_symm_mat(numbers, N):
    matrix = np.zeros((N, N))
    k = 0
    for i in range(N):
        for j in range(i+1):
            matrix[i][j] = float(numbers[k])
            k += 1
    
    mat_res = matrix + matrix.T - np.diag(np.diag(matrix))
    # print(len(mat_res), "pppppppppppppppppppppppppppp")
    return mat_res

length = int(NBAS*(NBAS+1)/2)

over_num = read_and_extract_mat(filename, "$OVERLAP" , length)
overlap_mat = create_symm_mat(over_num, NBAS)

###############################################################################


# orbital_file = 'examples/close-shell/h2o-b3lyp-def2svp.33'
# orbital_file = "h2o-b3lyp-def2svp.33"#input("Enter NBO key file to get orbitals: ")

def get_cmos():
    try:
        with open(orbital_file, "r") as file:
            lines = file.readlines()
    
            # Check if the file has at least two lines
            if len(lines) >= 2:
                second_line = lines[1].strip()
                orbital_type = second_line.split()[0]
                # print(orbital_type, "in AO basis")
    
            if "ALPHA" in lines[3].strip():
                print("This is an open-shell system")
                alpha_or_beta = input("Enter A or B to select Alpha or Beta spin: ")
    
                if alpha_or_beta.lower() == 'a':
                    start_line = 4
                elif alpha_or_beta.lower() == 'b':
                    # Initialize a variable to store the line number
                    for i, line in enumerate(lines):
                        if "BETA" in line:
                            start_line = i + 1  # Adding 1 because line numbers start from 1
                            break
            else:
                print("This is a closed shell system...")
                start_line = 3
    
            words = []
            for line in lines[start_line:]:
                line_elems = line.split()
    
                for elem in line_elems:
                    try:
                        if len(words) < NBAS * NBAS:
                            float_elem = float(elem)
                            words.append(float_elem)
    
                    except ValueError:
                        pass
    
            num_cmos = int(len(words) / NBAS)
            # print(len(words), num_cmos)
            orbital_arr = np.array(words).reshape(NBAS, num_cmos)
            
            # print(pd.DataFrame(orbital_arr))
            
            
            # x = 23
            # s = overlap_matdd
            # c = cmos
            # # print(pd.DataFrame(cmos.T @ overlap_matdd @ cmos))
            # print(c[x].T @ s @ c[x])
            # print(pd.DataFrame(s))
            
    except FileNotFoundError:
        print("The file does not exist.")
    
    except IOError: 
        print("Error reading the file.")

    return orbital_arr

orbital_arr = get_cmos()

# Setting the display option for float format
# pd.options.display.float_format = '{:.4f}'.format
# Setting display options to show all rows
pd.set_option('display.max_rows', None)
def pre_orb_from_orb(orb_mat, over_mat):    
    # print("\n\n")
    # print(pd.DataFrame(orb_mat[73:74,73]))
    # print("\n\n")
    # 
    # print("---------------------------------------")
    # res = orb_mat.T @ over_mat @ orb_mat
    # print(pd.DataFrame(np.diag(res)))
    
    # print(pd.DataFrame(orb_mat))
    
    res = np.dot(orb_mat.T, over_mat)
    res = np.dot(res, orb_mat)
    # print(pd.DataFrame(np.diag(res)))
    # print("---------------------------------------")
    
    
    # Initialize an empty list to store the start and end indices
    blow_bhigh = []
    
    # Initialize variables to track the current center and start index
    current_center = basis_info_dict[0]['CENTER']
    start_index = 1
    
    # Iterate over the info dictionaries in basis_info_dict
    for i, info in enumerate(basis_info_dict, start=1):
        # If the center has changed, append the previous range to blow_bhigh
        # and update current_center and start_index
        if info['CENTER'] != current_center:
            blow_bhigh.append({'start': start_index, 'end': i-1})
            current_center = info['CENTER']
            start_index = i
    
    # Append the last range to blow_bhigh
    blow_bhigh.append({'start': start_index, 'end': i})
    
    # Convert start_end to 0-based indexing
    start_ends = [{'start': se['start']-1, 'end': se['end']-1} for se in blow_bhigh]
    # print(start_ends)
    
    # Create a mask of all zeros with the same shape as nums
    mask = np.zeros_like(orb_mat)

    # Set the elements at the indices specified in start_end to 1 in the mask
    for se in start_ends:
        mask[se['start']:se['end']+1, se['start']:se['end']+1] = 1

    # # Multiply nums by the mask to get the modified array
    # print(pd.DataFrame(mask), 'mask ===================\n')
    # print(pd.DataFrame(zero_array), "zero_array++++++++++++++++++++\n")
    
    # new_mask =  zero_array * mask
    # print(pd.DataFrame(new_mask), "new_mask\n")
    
    # print(pd.DataFrame(orb_mat), "orb_mat ####################\n")
    
    # orb_typ_mul_new_mask = orb_mat * new_mask
    
    # print(pd.DataFrame(orb_typ_mul_new_mask), "%%%%%%%%%%%%%%%%%%%\n")
    
    
    
    
    repaired_orb_type = orb_mat * mask
    # print(pd.DataFrame(repaired_orb_type), "repaired_orb_type")
    
    # print(mask)
    
    # print(pd.DataFrame(repaired_orb_type))
    # print("\n\n")
    
    # repaired_orb_type = orb_mat
    # repaired_orb_type[117:124] = repaired_orb_type[117:124] * 0
    # print(pd.DataFrame(np.array(repaired_orb_type[82][83:124])))

    
    res = repaired_orb_type.T @ over_mat @ repaired_orb_type
    

    
    # Get the diagonal of res
    diag_res = np.diag(res)

    # Compute the inverse square root of the diagonal elements
    inv_sqrt_diag_res = 1 / np.sqrt(diag_res)
    # print(inv_sqrt_diag_res)
    
    # print(pd.DataFrame(inv_sqrt_diag_res), "inv_sqrt_diag_res \n")
    
    # print(pd.DataFrame(repaired_orb_type), "repaired_orb_type\n")
    
    renorm_repaired_orb = repaired_orb_type * inv_sqrt_diag_res
    
    # print(pd.DataFrame(renorm_repaired_orb), "renorm_repaired_orb\n")
    
    # print(pd.DataFrame(renorm_repaired_orb), "renorm_repaired_orb !!!!!!!!!!!!!!!!!!\n")
    
    
    # c = renorm_repaired_orb
    # c = repaired_orb_type
    # s = over_mat
    
    # print(pd.DataFrame(np.diag(c.T @ s @ c)))
    # print(pd.DataFrame(renorm_repaired_orb))
    
    
    return renorm_repaired_orb.T





pnao = pre_orb_from_orb(orbital_arr.T, overlap_mat)
# print(pd.DataFrame(pnao))
# print(renorm_repaired_orb)


# print(pnao[0], len(pnao), type(pnao), pnao.shape)

pnao = pd.DataFrame(pnao).T
# print(pnao)

# fn = "h2o-b3lyp-def2svp.32"
fn = output_file
with open(fn, 'w') as file:
    file.write(f"""NAO2PNAO file: {fn}
Orbitals in AO Basis
-------------------------------------------------------------------------------
""")   
   
   
    for col in pnao.columns:
        column_data = pnao[col].values.tolist()
        for i, value in enumerate(column_data):
            file.write(str("{:.13e}".format(value)))
            if (i + 1) % 5 == 0:
                file.write('\n')
            else:
                file.write(' ')
        file.write('\n')
print( fn, "created and saved...\n")


# 
