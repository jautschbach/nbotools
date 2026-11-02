
import pandas as pd
import numpy as np
from basis_info import atom_data, basis_info_dict
from sep_mat import is_open, aoX_trans,  opera_matrix, filename, ao_basis_list_dict, ao_basis
import re
import builtins
import sys
<<<<<<< HEAD
from icecream import ic
=======
#from icecream import ic
>>>>>>> 03b3a81 (Citation and menu list updated)



def process_nbo46(file_path, start_marker1="NBO", start_marker2=None, atom_index=1):
    copy_lines = []
    copying = False
    
    with open(file_path, 'r') as file:
        for line in file:
            if not copying:
                if start_marker2:
                    if start_marker1 in line and start_marker2 in line:
                        copying = True
                    continue
                else:
                    if start_marker1 in line:
                        copying = True
                    continue
            else:
                words = line.split()
                if len(words) < 4:
                    if ')' not in line:
                        break
                copy_lines.append(line.rstrip('\n'))

    def combine_labels(words):
        combined = []
        i = 0
        while i < len(words):
            word = words[i]
            if '-' in word:
                label = word
                i += 1
                while i < len(words):
                    label += words[i]
                    if re.search(r'[\d*]$', words[i]):
                        if i+1 < len(words) and re.match(r'^[A-Za-z]', words[i+1]):
                            break
                    i += 1
                combined.append(label)
                i += 1
            else:
                label = word
                while not (label.endswith(')') or label.endswith('*')):
                    i += 1
                    if i >= len(words):
                        break
                    label += words[i]
                combined.append(label)
                i += 1
        return combined

    def split_labels(concatenated):
        pattern = r'(?<=[\d*)])(?=[A-Za-z])'
        return re.split(pattern, concatenated)

    all_labels = []
    for line in copy_lines:
        words = line.split()
        combined_labels = combine_labels(words)
        for label in combined_labels:
            all_labels.extend(split_labels(label))

    def find_labels_with_indices(all_labels, atom_index):
        result = []
        pattern = re.compile(r'^[A-Za-z]+(\d+)\(')
        for idx, label in enumerate(all_labels, start=1):
            match = pattern.match(label)
            if match and match.group(1) == str(atom_index):
                result.append((idx, label))
        return result

    cr_lp_pattern = re.compile(r'\((cr|lp)\)')
    ry_pattern = re.compile(r'\(ry\)')

    matches = find_labels_with_indices(all_labels, atom_index)


    atom_cr_lp_idx = [idx for idx, label in matches if cr_lp_pattern.search(label)]
    atom_ry_idx    = [idx for idx, label in matches if ry_pattern.search(label)]

    all_ry_idx     = {idx for idx, label in enumerate(all_labels, start=1) if ry_pattern.search(label)}
    all_cr_lp_idx  = set(atom_cr_lp_idx)

    other_idx      = [idx for idx, label in enumerate(all_labels, start=1)
                      if idx not in all_ry_idx and idx not in all_cr_lp_idx]

    non_atom_ry_idx = [idx for idx, label in enumerate(all_labels, start=1)
                       if ry_pattern.search(label) and idx not in set(atom_ry_idx)]
    
    return atom_cr_lp_idx, atom_ry_idx, other_idx, non_atom_ry_idx

unique_atoms = []
for label in ao_basis:
    atom = label  
    if atom not in unique_atoms:
        unique_atoms.append(atom)

atom_labels = {i+1: f"{atom.rstrip('0123456789')} {i+1}" for i, atom in enumerate(unique_atoms)}

def get_angmom_ranges(basis_info_dict):
    atom_angmom_ranges = {}
    
    for index, info in enumerate(basis_info_dict, start=0):
        center = info["CENTER"]
        orb_type = info["type"].upper()  # Normalize type names
        
        # Map orbital types to labels
        if orb_type.startswith("S"):
            angmom = "S"
        elif orb_type.startswith("P"):
            angmom = "P"
        elif orb_type.startswith("D"):
            angmom = "D"
        elif orb_type.startswith("F"):
            angmom = "F"
        elif orb_type.startswith("G"):
            angmom = "G"
        elif orb_type.startswith("H"):
            angmom = "H"
        elif orb_type.startswith("I"):
            angmom = "I"
        elif orb_type.startswith("J"):
            angmom = "J"
        else:
            continue  
        
        if center not in atom_angmom_ranges:
            atom_angmom_ranges[center] = {}
        
        if angmom not in atom_angmom_ranges[center]:
            atom_angmom_ranges[center][angmom] = (index, index)
        else:
            atom_angmom_ranges[center][angmom] = (atom_angmom_ranges[center][angmom][0], index)
    
    return atom_angmom_ranges
atom_angmom_ranges = get_angmom_ranges(basis_info_dict)


def nlmo_atom_angmom_partition_normalized(orb_mat, s_mat):
    sc = s_mat @ orb_mat
    nloc = orb_mat.shape[1]
    results = []
    for iorb in range(nloc):
        atom_contribs = {}
        for atom, angmom_dict in atom_angmom_ranges.items():
            angmom_contribs = {}
            atom_total = 0.0
            for angmom, (lo, hi) in angmom_dict.items():
                val = sum(
                    orb_mat[u, iorb] * sc[u, iorb]
                    for u in range(lo, hi + 1)
                )
                angmom_contribs[angmom] = abs(val)
                atom_total += abs(val)
                
            for angmom in angmom_dict:
                if atom_total > 0:
                    angmom_contribs[angmom] = 100 * angmom_contribs[angmom] / atom_total
                else:
                    angmom_contribs[angmom] = 0.0
            angmom_contribs['total'] = atom_total * 100
            atom_contribs[atom] = angmom_contribs
        results.append(atom_contribs)
    return results

def print(*args, **kwargs):
    if 'file' not in kwargs:
        kwargs['file'] = sys.stderr
    builtins.print(*args, **kwargs)


file_path = filename + ".46"

def print_atom_populations_table(atom_populations_per_orbital, cutoff=1e-4, top_n=None):
    max_pop_warnings = 0
    total_warnings = 0

    for i, popdict in enumerate(atom_populations_per_orbital, start=1):
        items = [(atom, pop) for atom, pop in popdict.items() if abs(pop) >= cutoff]
        items.sort(key=lambda x: x[1], reverse=True)

        if top_n is not None:
            items = items[:top_n]

        total = sum(popdict.values())
        max_pop = max(popdict.values()) if popdict else 0.0

        print(f"\nOrbital {i:>3}   total = {total:.6f}")
        if max_pop > 1.00:
            print(f"WARNING: contribution exceeds 100% (max = {max_pop:.6f})")
            max_pop_warnings += 1
        if round(total, 2) > 1.00:
            print(f"WARNING: orbital total exceeds 100% (total = {total:.6f})")
            total_warnings += 1

        print("-" * 42)
        print(f"{'Atom':>6}  {'Population':>14}")
        print("-" * 42)

        for atom, pop in items:
            print(f"{atom:>6}  {pop:>14.8f}")

        if not items:
            print("(no contributions above cutoff)")

    print("\nSummary")
    print("-" * 42)
    print(f"Orbitals with max contribution > 100%: {max_pop_warnings}")
    print(f"Orbitals with total > 100%:          {total_warnings}")



def analys_nlmo(c, s, density, fock, occ, m_cr_lp, m_ry, l_don_nlmo, l_ry, tmi):
    
    nuclear_charge = atom_data[tmi - 1][1]
   
    basis = [{}] + ao_basis_list_dict  # Index basis from 1
    natom = len(ao_basis_list_dict)

    n_form_occ = int(round(np.trace(occ) if is_open else np.trace(occ) / 2))
    nloc =  int(s.shape[0])
    iloc = np.arange(nloc + 1, dtype=int)

    

    orb_occ = np.diag(occ)
    sc = s @ c

    # Precompute atom populations for all orbitals
    atom_populations_per_orbital = []
    for ss in range(1, nloc + 1):
        s_index = iloc[ss]
        atom_populations = {}
        for atom in range(1, natom + 1):
            qas = sum(
                c[u, s_index - 1] * sc[u, s_index - 1]
                for u in range(basis[atom]['bflo'], basis[atom]['bfhi'] + 1)
            )
            atom_populations[atom] = abs(qas)
        atom_populations_per_orbital.append(atom_populations)
        

    atom_populations_per_orb = atom_populations_per_orbital[:]
                 
    tm_atom_idx = tmi
    tm_total_pop = sum(
        atom_populations_per_orbital[idx][tm_atom_idx] * orb_occ[idx]
        for idx in range(len(atom_populations_per_orbital))
    )
    
    print(f"Total electron occupation on transition metal (all NLMOs): {tm_total_pop:.6f}")


  
    print("""The NLMOs are partitioned into Four Groups: m_cr_lp, l_don_nlmo, m_ry and l_ry.""")
   
    
    m_cr_lp_occ = sum(
        orb_occ[idx - 1] for idx in m_cr_lp
    )
    
    m_ry_occ = sum(
        orb_occ[idx - 1] for idx in m_ry
    )
    
    
    l_don_nlmo_occ = sum(
        orb_occ[idx - 1] for idx in l_don_nlmo
    )
    
    l_ry_occ = sum(
        orb_occ[idx - 1] for idx in l_ry
    )
    
    tot_e = (m_ry_occ + m_cr_lp_occ + l_ry_occ  + l_don_nlmo_occ)
    
    print(f"NLMO occupations: Metal core and LP ={m_cr_lp_occ:.6f}, metal rydberg={m_ry_occ:.6f}, ligand donor={l_don_nlmo_occ:.6f}, ligand rydberg={l_ry_occ:.6f}")





    print(f"Total electron count obtained from summing the four NLMO groups: {tot_e:.6f}")
    # print(m_ry_occ + m_cr_lp_occ + l_ry_occ  + l_don_nlmo_occ)

    
    sigma_donated_elect = sum(
        atom_populations_per_orbital[idx - 1][tm_atom_idx] * orb_occ[idx - 1]
        for idx in l_don_nlmo
    )
    print(f"Sigma Donated electrons = {sigma_donated_elect:.6f}")
    
    pi_donated_elect = sum(
        sum(
            atom_populations_per_orbital[idx - 1][atom_idx]
            for atom_idx in range(1, natom + 1) if atom_idx != tm_atom_idx
        ) * orb_occ[idx - 1]
        for idx in m_cr_lp
        )
    
    print(f"Pi Donated electrons = {pi_donated_elect:.6f}")
    
    
    print('\nPost Validation')
    if is_open:
        nlmo_charge = nuclear_charge/2 - tm_total_pop
        print(f"NLMO charge (nuclear_charge/2 - tm_total_pop) = {nuclear_charge/2 - tm_total_pop:.6f}")
    else:
        nlmo_charge = nuclear_charge - tm_total_pop
        print(f"NLMO charge (nuclear_charge - tm_total_pop) = {nuclear_charge - tm_total_pop:.6f}")
        
    print("Note: NLMO charge is the same as the NPA charge for HF and DFT calculations.")


    
    l_cr_lp = [x for x in l_don_nlmo if x <= n_form_occ]
    
    
    print("""
Calculate the electron loss or gained due to correlation on to the Rydberg  Orbitals.""")
    if is_open:
        M_corr =  m_cr_lp_occ + m_ry_occ - len(m_cr_lp)
        L_corr =  l_don_nlmo_occ + l_ry_occ - len(l_cr_lp)
        
        print(f"M_corr (m_cr_lp_occ + m_ry_occ - len(m_cr_lp)) = {m_cr_lp_occ + m_ry_occ - len(m_cr_lp):.6f}")
        print(f"L_corr (l_don_nlmo_occ + l_ry_occ - len(l_cr_lp)) = {l_don_nlmo_occ + l_ry_occ - len(l_cr_lp):.6f}")

    else:
        M_corr =  m_cr_lp_occ + m_ry_occ - 2 * len(m_cr_lp)
        L_corr =  l_don_nlmo_occ + l_ry_occ - 2* len(l_cr_lp)
        
        print(f"M_corr (m_cr_lp_occ + m_ry_occ - 2 * len(m_cr_lp)) = {m_cr_lp_occ + m_ry_occ - 2 * len(m_cr_lp):.6f}")
        print(f"L_corr (l_don_nlmo_occ + l_ry_occ - 2 * len(l_cr_lp)) = {l_don_nlmo_occ + l_ry_occ - 2 * len(l_cr_lp):.6f}")

    
    print("\nGet oxidation number")
    
    oxid_num = sigma_donated_elect + nlmo_charge - pi_donated_elect + M_corr
    print(f"oxid_num (sigma_donated_elect + nlmo_charge - pi_donated_elect + M_corr) = {sigma_donated_elect + nlmo_charge - pi_donated_elect + M_corr:.6f}")
    
    print("\nMetal electron count")
    print(f"metal_e_count (m_cr_lp_occ + m_ry_occ + sigma_donated_elect - pi_donated_elect) = {m_cr_lp_occ + m_ry_occ + sigma_donated_elect - pi_donated_elect:.6f}")
    
    
    
    
    ############################
    
    M_in_l_ry = sum(
        atom_populations_per_orbital[idx - 1][tm_atom_idx] * orb_occ[idx - 1]
        for idx in l_ry
    )
    print(f"Metal Population in Ligand RY = {M_in_l_ry:.6f}")
    
    
    L_in_M_ry = sum(
        sum(
            atom_populations_per_orbital[idx - 1][atom_idx]
            for atom_idx in range(1, natom + 1) if atom_idx != tm_atom_idx
        ) * orb_occ[idx - 1]
        for idx in m_ry
        )
    
    print(f"Ligand population in metal RY = {L_in_M_ry:.6f}")
    

    return pi_donated_elect, sigma_donated_elect, tm_total_pop, atom_populations_per_orb, oxid_num
      

def extract_matrices(opera_matrix, aoX_trans, spin=None):
    if spin == 'alpha':
        fao = pd.DataFrame(opera_matrix["FOCK_ALPHA"]).values
        dm  = pd.DataFrame(opera_matrix["DENSITY_ALPHA"]).values
        nlmo = pd.DataFrame(aoX_trans["AONLMOs_ALPHA"]).values
        aonao = pd.DataFrame(aoX_trans["AONAOs_ALPHA"]).values
    elif spin == 'beta':
        fao = pd.DataFrame(opera_matrix["FOCK_BETA"]).values
        dm  = pd.DataFrame(opera_matrix["DENSITY_BETA"]).values
        nlmo = pd.DataFrame(aoX_trans["AONLMOs_BETA"]).values
        aonao = pd.DataFrame(aoX_trans["AONAOs_BETA"]).values
    else:
        fao = pd.DataFrame(opera_matrix["FOCK"]).values
        dm  = pd.DataFrame(opera_matrix["DENSITY"]).values
        nlmo = pd.DataFrame(aoX_trans["AONLMOs"]).values
        aonao = pd.DataFrame(aoX_trans["AONAOs"]).values
    sao = pd.DataFrame(opera_matrix["OVERLAP"]).values
    return sao, fao, dm, nlmo, aonao

def get_occupation(nlmo, dm):
    return np.linalg.inv(nlmo) @ dm @ np.linalg.inv(nlmo.T)

def get_naonlmo(aonao, nlmo):
    return np.linalg.inv(aonao) @ nlmo

def transform_to_nao(aonao, matrix):
    return aonao.T @ matrix @ aonao

def compute_charges(atom_data, ao_basis_list_dict, nao_alpha_occ=None, nao_beta_occ=None, closed_shell=False):
    charges = []
    # for i, (atom, ao_range) in enumerate(zip(atom_data, ao_basis_list_dict)):
    for atom, ao_range in zip(atom_data, ao_basis_list_dict):
        atomic_num = atom[0]
        nuclear_charge = atom[1]
        idx_lo = ao_range['bflo']
        idx_hi = ao_range['bfhi'] + 1
        if closed_shell:
            occ = np.sum(nao_alpha_occ[idx_lo:idx_hi])
            charge = nuclear_charge - occ
            charges.append({
                'atomic_num': atomic_num, 'nuclear_charge': nuclear_charge,
                'charge': charge, 'net_charge': charge,
                'occ': occ, 'total_nao_occ': occ
            })
        else:
            occ_a = np.sum(nao_alpha_occ[idx_lo:idx_hi]) if nao_alpha_occ is not None else 0
            occ_b = np.sum(nao_beta_occ[idx_lo:idx_hi])  if nao_beta_occ is not None else 0
            total_occ = occ_a + occ_b
            alpha_charge = (nuclear_charge / 2) - occ_a
            beta_charge = (nuclear_charge / 2) - occ_b
            net_charge = nuclear_charge - total_occ
            charges.append({
                'atomic_num': atomic_num, 'nuclear_charge': nuclear_charge,
                'alpha_charge': alpha_charge, 'beta_charge': beta_charge,
                'net_charge': net_charge, 'occ_alpha': occ_a, 'occ_beta': occ_b,
                'total_nao_occ': total_occ
            })
    return charges

def debug_block(debug_data):
    print("\n--- Debug Information ---")

    label_width = max(len(k) for k in debug_data)
    value_width = 10

    for key, value in debug_data.items():
        if isinstance(value, (int, float)):
            print(f"{key:<{label_width}} : {value:>{value_width}.5f}")
        else:
            print(f"{key:<{label_width}} : {value}")

    print("-" * (label_width + value_width + 3) + "\n")

def part_run():
    
    transition_metal_index = int(input("Enter index of transition metal: "))
    print(f"Index of first transition metal atom: {transition_metal_index}")
    
    
    try:
        if is_open:
            print("\nOpen shell system")
            print("Processing Alpha Spin")
            # Alpha extraction & calculation
            sao_a, fao_a, dm_a, nlmo_a, aonao_a = extract_matrices(opera_matrix, aoX_trans, 'alpha')
            occupation_a = get_occupation(nlmo_a, dm_a)
            naonlmo_a = get_naonlmo(aonao_a, nlmo_a)
            snao_a = transform_to_nao(aonao_a, sao_a)
            dmnao_a = transform_to_nao(aonao_a, dm_a)
            fnao_a = transform_to_nao(aonao_a, fao_a)
    
            m_cr_lp_a, m_ry_a, l_cr_lp_a, l_ry_a = process_nbo46(file_path, "NBO", "ALPHA", transition_metal_index)      
            a_pi_don, a_sigma_don, a_tot_nlmo, a_atom_pop_orb, a_oxid_num = analys_nlmo(naonlmo_a, snao_a, dmnao_a, fnao_a, occupation_a, m_cr_lp_a, m_ry_a, l_cr_lp_a, l_ry_a, transition_metal_index)

            
            
            print("\nProcessing Beta Spin")
            sao_b, fao_b, dm_b, nlmo_b, aonao_b = extract_matrices(opera_matrix, aoX_trans, 'beta')
            occupation_b = get_occupation(nlmo_b, dm_b)
            naonlmo_b = get_naonlmo(aonao_b, nlmo_b)
            snao_b = transform_to_nao(aonao_b, sao_b)
            dmnao_b = transform_to_nao(aonao_b, dm_b)
            fnao_b = transform_to_nao(aonao_b, fao_b)
            m_cr_lp_b, m_ry_b, l_cr_lp_b, l_ry_b = process_nbo46(file_path, "NBO", "BETA", transition_metal_index)      
            b_pi_don, b_sigma_don, b_tot_nlmo, b_atom_pop_orb, b_oxid_num = analys_nlmo(naonlmo_b, snao_b, dmnao_b, fnao_b, occupation_b, m_cr_lp_b, m_ry_b, l_cr_lp_b, l_ry_b, transition_metal_index)
            
            print("\nProcessing Data")
            nao_alpha_occ = np.diag(np.linalg.inv(aonao_a) @ dm_a @ np.linalg.inv(aonao_a.T))
            nao_beta_occ  = np.diag(np.linalg.inv(aonao_b) @ dm_b @ np.linalg.inv(aonao_b.T))
            charges = compute_charges(atom_data, ao_basis_list_dict, nao_alpha_occ=nao_alpha_occ, nao_beta_occ=nao_beta_occ, closed_shell=False)
    
            tm_info = next(
                c for c in charges
                if c['atomic_num'] == atom_data[transition_metal_index - 1][0]
            )
            
          
            # Compile calculation results
            total_sigma_donation = a_sigma_don + b_sigma_don
            total_pi_donation    = a_pi_don + b_pi_don
            total_nlmo_occ       = a_tot_nlmo + b_tot_nlmo
            total_nao_occ        = tm_info['total_nao_occ'] #sum(tm['total_nao_occ'] for tm in transition_metal_charges)
            deviat               = total_nao_occ - total_nlmo_occ
            net_ON = a_oxid_num + b_oxid_num
            
            debug_block({
                'Alpha L-M donation': a_sigma_don,
                'Alpha M-L donation': a_pi_don,
                'Beta L-M donation': b_sigma_don,
                'Beta M-L donation': b_pi_don,
                'Total L-M donation': total_sigma_donation,
                'Total M-L donation': total_pi_donation,
                'Alpha NPA charge': tm_info['alpha_charge'],
                'Beta NPA charge': tm_info['beta_charge'],
                'Net NPA charge': tm_info['net_charge'],
                'Alpha O.N.':  a_oxid_num,
                'Beta O.N.': b_oxid_num,
                'Net O.N.': net_ON,                
                'total NAO occ. - total NLMO occ.': deviat,
            })
            
           
        
        else:
            
            print("\nClosed shell system")
            sao, fao, dm, nlmo, aonao = extract_matrices(opera_matrix, aoX_trans)
            occupation = get_occupation(nlmo, dm)
            naonlmo = get_naonlmo(aonao, nlmo)
            snao = transform_to_nao(aonao, sao)
            dmnao = transform_to_nao(aonao, dm)
            fnao = transform_to_nao(aonao, fao)
            
            m_cr_lp, m_ry, l_cr_lp, l_ry = process_nbo46(file_path, "NBO", None, transition_metal_index)
            pi_don, sigma_don, tot_nlmo_occ, atom_pop_orb, oxid_num  = analys_nlmo(naonlmo, snao, dmnao, fnao, occupation, m_cr_lp, m_ry, l_cr_lp, l_ry, transition_metal_index)

            
            print("\nProcessing Data")
            nao_occ_arr = np.diag(np.linalg.inv(aonao) @ dm @ np.linalg.inv(aonao.T))
            charges = compute_charges(atom_data, ao_basis_list_dict, nao_alpha_occ=nao_occ_arr, closed_shell=True)

            tm_info = next(
            c for c in charges
            if c['atomic_num'] == atom_data[transition_metal_index - 1][0]
            )

            total_sigma_donation = sigma_don
            total_pi_donation    = pi_don
            total_nlmo_occ       = tot_nlmo_occ
            total_nao_occ = tm_info['occ']
            
            total_nao_occ        = tm_info['total_nao_occ'] 
            deviat               = total_nao_occ - total_nlmo_occ

            debug_block({
                'L-M donation': sigma_don,
                'M-L donation': pi_don,
                'NPA charge': tm_info['net_charge'],
                'ON': oxid_num,
                'NAO occ. - NLMO occ.': deviat,
            })
            
    
    except Exception as e:
        print(f"Error: {e}")
        
        print("0,0,0,0,0,0,0,0,0,0,0" if is_open else "0,0,0,0,0,", file=sys.stdout)


