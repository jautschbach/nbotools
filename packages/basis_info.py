
import numpy as np
import re
from itertools import groupby
from scipy.constants import physical_constants
import os


def parse_file47(filename):
    print(f"Parsing {filename} as a .47 file")

    def read_file47(filename):
        try:
            with open(filename, 'r') as file:
                return file.read()
        except FileNotFoundError:
            raise FileNotFoundError(f"File '{filename}' not found.")
        except IOError:
            raise IOError(f"Error reading file '{filename}'.")

    file_content = read_file47(filename)

    def parse_array_from_block(varname, content, dtype=float):
        pattern = re.compile(rf'{varname}\s+=\s+((?:[-+]?\d+\.\d+(?:E[+-]?\d+)?\s+)+)')
        matches = pattern.findall(content)
        values = []
        for match in matches:
            parts = match.split()
            converted = [dtype(v) for v in parts]
            values.extend(converted)
        return values

    def parse_int_array(varname, content):
        pattern = re.compile(rf'{varname}\s+=\s+([\d\s]+)')
        matches = pattern.findall(content)
        values = []
        for match in matches:
            parts = match.split()
            converted = [int(v) for v in parts]
            values.extend(converted)
        return values

    def is_ungrouped(ncomp):
        return all(x == 1 for x in ncomp)

    def process_orbital_labels(label, ncomp, orb_mapping):
        if not sum(ncomp) == len(label):
            raise ValueError(f"NCOMP sum ({sum(ncomp)}) does not match LABEL length ({len(label)})")
        
        orb_type = [orb_mapping.get(l, ('unknown', 'unknown'))[0] for l in label]
        orb_val = [orb_mapping.get(l, ('unknown', 'unknown'))[1] for l in label]
        shell_num = []

        # Primary shell assignment based on NCOMP
        idx = 0
        for shell_idx, nc in enumerate(ncomp, 1):
            group = orb_type[idx:idx + nc]
            if nc > 1:
                # Validate that orbitals in the group share the same base type
                base_type = group[0][0] if group else None
                if not all(t[0] == base_type for t in group) or len(group) != nc:
                    raise ValueError(f"Invalid NCOMP grouping at shell {shell_idx}: {group} does not match NCOMP={nc}")
            for _ in range(nc):
                shell_num.append(shell_idx)
            idx += nc

        if is_ungrouped(ncomp):
            shell_num = list(range(1, len(label) + 1))
        else:
            # Secondary validation with type-based grouping
            type_limits = {'p': 3, 'd': 5, 'f': 7, 'g': 9, 'h': 11, 'i': 13, 'j': 15}
            temp_shell_num = []
            count = 1
            orb_count = 0
            for i, t in enumerate(orb_type):
                base_type = t[0]
                if i > 0 and base_type in type_limits and orb_type[i-1][0] == base_type:
                    if orb_count < type_limits[base_type]:
                        temp_shell_num.append(temp_shell_num[-1])
                        orb_count += 1
                    else:
                        count += 1
                        temp_shell_num.append(count)
                        orb_count = 1
                else:
                    count += 1
                    temp_shell_num.append(count)
                    orb_count = 1

            # Cross-validate NCOMP-based and type-based shell assignments
            ncomp_shells = [len(list(g)) for k, g in groupby(shell_num)]
            type_shells = [len(list(g)) for k, g in groupby(temp_shell_num)]
            if ncomp_shells != type_shells:
                print(f"Warning: NCOMP-based shells {ncomp_shells} differ from type-based shells {type_shells}. Using NCOMP-based.")

        return orb_type, orb_val, shell_num

    def system_info(content):
        bohr_to_ang = physical_constants['Bohr radius'][0] * 1e10  # Bohr to Ångström
        use_bohr = "BOHR" in content.upper()
        to_bohr = 1 if use_bohr else bohr_to_ang

        # Parse atom coordinates: (atomic_number, charge, x, y, z)
        coord_pattern = r'\s+(\d+)\s+(\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)'
        atom_matches = re.findall(coord_pattern, content)
        atom_data = [
            (int(z), int(chg), float(x), float(y), float(z_))
            for z, chg, x, y, z_ in atom_matches
        ]

        # Parse basis data
        variable_names_float = ['EXP', 'CS', 'CP', 'CD', 'CF', 'CG', 'CH', 'CI', 'CJ']
        variable_names_int = ['CENTER', 'LABEL', 'NSHELL', 'NEXP', 'NCOMP', 'NPRIM', 'NPTR']
        
        parsed_float = {var: parse_array_from_block(var, content, float) for var in variable_names_float}
        parsed_int = {var: parse_int_array(var, content) for var in variable_names_int}

        EXP = parsed_float['EXP']
        CS, CP, CD, CF, CG, CH, CI, CJ = [parsed_float.get(v, []) for v in variable_names_float[1:]]
        CENTER = parsed_int['CENTER']
        LABEL = parsed_int['LABEL']
        NCOMP = parsed_int['NCOMP']
        NPRIM = parsed_int['NPRIM']
        NPTR = parsed_int['NPTR']

        # Orbital mapping
        orb_mapping = {
            1: ('s', 's'), 51: ('s', 's'), 101: ('px', 'px'), 102: ('py', 'py'), 103: ('pz', 'pz'),
            151: ('px', 'px'), 152: ('py', 'py'), 153: ('pz', 'pz'),
            251: ('d_xy', 'ds2'), 252: ('d_xz', 'ds1'), 253: ('d_yz', 'dc1'),
            254: ('d_x2-y2', 'dc2'), 255: ('d_z2', 'd0'),
            351: ('fz(5z2-3r2)', 'f0'), 352: ('fx(5z2-r2)', 'fc1'), 353: ('fy(5z2-r2)', 'fs1'),
            354: ('fz(x2-y2)', 'fc2'), 355: ('fxyz', 'fs2'), 356: ('fx(x2-3y2)', 'fc3'),
            357: ('f(3x2-y2)', 'fs3'),
            451: ('g0', 'g0'), 452: ('gc1', 'gc1'), 453: ('gs1', 'gs1'), 454: ('gc2', 'gc2'),
            455: ('gs2', 'gs2'), 456: ('gc3', 'gc3'), 457: ('gs3', 'gs3'), 458: ('gc4', 'gc4'),
            459: ('gs4', 'gs4'),
            551: ('h0', 'h0'), 552: ('hc1', 'hc1'), 553: ('hs1', 'hs1'), 554: ('hc2', 'hc2'),
            555: ('hs2', 'hs2'), 556: ('hc3', 'hc3'), 557: ('hs3', 'hs3'), 558: ('hc4', 'hc4'),
            559: ('hs4', 'hs4'), 560: ('hc5', 'hc5'), 561: ('hs5', 'hs5'),
            651: ('i0', 'i0'), 652: ('ic1', 'ic1'), 653: ('is1', 'is1'), 654: ('ic2', 'ic2'),
            655: ('is2', 'is2'), 656: ('ic3', 'ic3'), 657: ('is3', 'is3'), 658: ('ic4', 'ic4'),
            659: ('is4', 'is4'), 660: ('ic5', 'ic5'), 661: ('is5', 'is5'), 662: ('ic6', 'ic6'),
            663: ('is6', 'is6'),
            751: ('j0', 'j0'), 752: ('jc1', 'jc1'), 753: ('js1', 'js1'), 754: ('jc2', 'jc2'),
            755: ('js2', 'js2'), 756: ('jc3', 'jc3'), 757: ('js3', 'js3'), 758: ('jc4', 'jc4'),
            759: ('js4', 'js4'), 760: ('jc5', 'jc5'), 761: ('js5', 'js5'), 762: ('jc6', 'jc6'),
            763: ('js6', 'js6'), 764: ('jc7', 'jc7'), 765: ('js7', 'js7')
        }

        orb_type, orb_val, shell_num = process_orbital_labels(LABEL, NCOMP, orb_mapping)

        # Map NPRIM and NPTR to basis functions based on shells
        if not is_ungrouped(NCOMP):
            if len(NPRIM) != len(NCOMP) or len(NPTR) != len(NCOMP):
                raise ValueError(f"NPRIM ({len(NPRIM)}) or NPTR ({len(NPTR)}) length does not match NCOMP ({len(NCOMP)})")
            # Replicate NPRIM and NPTR for each basis function in a shell
            NPRIM_expanded = []
            NPTR_expanded = []
            for shell_idx, nc in enumerate(NCOMP):
                NPRIM_expanded.extend([NPRIM[shell_idx]] * nc)
                NPTR_expanded.extend([NPTR[shell_idx]] * nc)
        else:
            NPRIM_expanded = NPRIM
            NPTR_expanded = NPTR

        # Validate lengths
        if len(NPRIM_expanded) != len(LABEL) or len(NPTR_expanded) != len(LABEL):
            raise ValueError(f"Expanded NPRIM ({len(NPRIM_expanded)}) or NPTR ({len(NPTR_expanded)}) does not match LABEL ({len(LABEL)})")

        # Build basis info dictionary
        bas_info_dict = []
        for i in range(len(LABEL)):
            atom_idx = CENTER[i] - 1
            prim = NPRIM_expanded[i]
            ptr = NPTR_expanded[i]

            info = {
                "N": i + 1,
                "CENTER": CENTER[i],
                "LABEL": LABEL[i],
                "shell_num": shell_num[i],
                "type": orb_type[i],
                "orb_val": orb_val[i],
                "exps": EXP[ptr - 1: ptr - 1 + prim]
            }

            coeffs = []
            for coeff_array in [CS, CP, CD, CF, CG, CH, CI, CJ]:
                slice_ = coeff_array[ptr - 1: ptr - 1 + prim]
                coeffs.extend([c for c in slice_ if c != 0.0])
            info["coeffs"] = coeffs

            atom_coord = atom_data[atom_idx][2:5]  # x, y, z
            info["xcenter"] = atom_coord[0] / to_bohr
            info["ycenter"] = atom_coord[1] / to_bohr
            info["zcenter"] = atom_coord[2] / to_bohr

            bas_info_dict.append(info)

        # Convert atom coordinates to Ångström for output
        atom_data_ang = [
            (z, charge, x * bohr_to_ang, y * bohr_to_ang, z_ * bohr_to_ang)
            if use_bohr else (z, charge, x, y, z_)
            for (z, charge, x, y, z_) in atom_data
        ]

        return bas_info_dict, atom_data_ang, to_bohr

    basis_info_dict, atom_data, to_bohr = system_info(file_content)

    # Prepare output
    coordinates = [atom[2:] for atom in atom_data]
    atom_info = [(atom[0],) + tuple(atom[2:]) for atom in atom_data]

    return basis_info_dict, coordinates, atom_data, to_bohr


def get_basis_info(filename):
    try:
        with open(filename, 'r') as file:
            content = file.read()
            return content
    except FileNotFoundError:
        print(f"File '{filename}' not found.")
    except IOError:
        print(f"Error reading file '{filename}'.")
        

filename = [f for f in os.listdir('.') if f.endswith('.47')][0]
basis_info_dict, coordinates, atom_data, to_bohr = parse_file47(filename)
data =  basis_info_dict



def pre_orb_from_orb(orb_mat, over_mat):    
    res = orb_mat.T @ over_mat @ orb_mat
    res = np.dot(orb_mat.T, over_mat)
    res = np.dot(res, orb_mat)
    
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
    
    # Create a mask of all zeros with the same shape as nums
    mask = np.zeros_like(orb_mat)

    # Set the elements at the indices specified in start_end to 1 in the mask
    for se in start_ends:
        mask[se['start']:se['end']+1, se['start']:se['end']+1] = 1

    # Multiply nums by the mask to get the modified array
    repaired_orb_type = orb_mat * mask
    
    res = repaired_orb_type.T @ over_mat @ repaired_orb_type
    
    # Get the diagonal of res
    diag_res = np.diag(res)

    # Compute the inverse square root of the diagonal elements
    inv_sqrt_diag_res = 1 / np.sqrt(diag_res)
    
    renorm_repaired_orb = repaired_orb_type * inv_sqrt_diag_res
    
    # c = repaired_orb_type
    # s = over_mat
    
    return renorm_repaired_orb