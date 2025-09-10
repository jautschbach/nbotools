
     
import numpy as np
import pandas as pd
from sep_mat import is_open
from sep_mat import aoX_trans

if is_open:
    AONAO = pd.DataFrame(aoX_trans["AONAOs_ALPHA"]).values
else:
    AONAO = np.array(aoX_trans.get('AONAOs'))

def pop_(loc_orb, smat, ao_basis_list_dict):
    num_col = loc_orb.shape[1]
    print("Number of Columns: ", num_col)
    sorted_c = loc_orb
    for orb in range(num_col):
        print(f"\n\nANALYSIS FOR ORBITAL {orb + 1}")
    
        c1 = sorted_c[:, orb]
        sc = smat @ c1
        
        # Compute atom-wise contributions and absolute values
        contributions = []
        abs_contributions = []
        contribution = 0
        
        for atom in ao_basis_list_dict:
            bflo = atom['bflo']
            bfhi = atom['bfhi']
            start_idx = bflo
            end_idx = bfhi + 1
            
            sc_m = sc[start_idx:end_idx]  
            c1_m = c1[start_idx:end_idx]   
            
            contribution = c1_m.T @ sc_m            
            contributions.append(contribution.item() if contribution.item() > 0 else 0)
            abs_contributions.append(abs(contribution.item()))  # Take absolute value
        
        # Calculate total absolute contribution
        total_abs_contribution = np.sum(abs_contributions)
        
        # Compute percentages
        percentages = [(abs_cont / total_abs_contribution) * 100 for abs_cont in abs_contributions]
        
        # Print results
        for i, atom in enumerate(ao_basis_list_dict):
            print(f"Atom {i+1} contribution: {contributions[i]:.6f}")
            print(f"Atom {i+1} percentage: {percentages[i]:.2f}%")
            print("-----")
            
def population_analy(c, s, basis_info, fock, NAO_mat=False):
    if NAO_mat == False:
        print("\n Orbitals in AO basis. Transforming to NAO basis...\n")
        
        NAOLMO = np.linalg.inv(AONAO) @ c
        SNAO = AONAO.T @ s @ AONAO #Should be a unit matrix
        FNAO = AONAO.T @ fock @ AONAO #Should be a unit matrix
        
        fock = FNAO
        c = NAOLMO
        s = SNAO
    else:
        #This is for when systems is done in reduced set and called
        print("\n Orbitals already in NAO basis. No transformations needed...\n")
        
        pass
                
    basis = [{}] + basis_info #We want to make basis to be indexed from 1
    natom = int(len(basis_info))       
    nloc = c.shape[1] 
    iloc = np.array([i for i in range(nloc + 1)], dtype=int)        
    list_size = 100000  # 
    list_array = [0] * list_size  # Create an empty list with size 10
    
    pop_size = 1000000  # Example size of the pop array
    pop_array = [0] * pop_size  # Create an empty list with size 10
    
    orb_energ = np.diag(c.T @ fock @ c)

    sc = s @ c
    print("Orbital  | Energy/Ha   | Atom: Percent contribution  ")
    print("---------+-------------+-----------------------------")   
    
    for ss in range(1, nloc + 1):
        s = iloc[ss]
        nlist = 0
        
        for a in range(1, natom + 1):
            qas = 0        
            bflo = basis[a]['bflo']
            bfhi = basis[a]['bfhi']
            
            for u in range(bflo, bfhi + 1):
                pass
                qas = qas + c[u , s - 1] * sc[u, s - 1]
            if abs(qas) > 0.001:
                nlist = nlist + 1
                list_array[nlist] = a
                pop_array[nlist] = qas
            
        for u in range(1, nlist+1):
            for t in range(1, u + 1):
                if abs(pop_array[t]) < abs(pop_array[u]):
                    tmp = pop_array[u]
                    pop_array[u] = pop_array[t]
                    pop_array[t] = tmp
                    tt = list_array[u]
                    list_array[u] = list_array[t]
                    list_array[t] = tt
       
        eval_array = orb_energ
        eval_energy = eval_array[s - 1]
        
        # Create a list of tuples and sort by the first element (the identifier)
        pop_dist = [(list_array[a], pop_array[a]) for a in range(1, nlist + 1)]
        pop_dist.sort(key=lambda x: x[0])  # Sort by the identifier (a)
        
        # Calculate electron distribution
        electron_dist = [(list_array[a], pop_array[a]) for a in range(1, nlist + 1)]
        electron_dist.sort(key=lambda x: x[0])  # Sort by atom number
        
        # Set negative values to zero 
        electron_dist_non_negative = [
            (t[0], max(0, t[1]))
            for t in electron_dist
        ]
        
        # Calculate the sum of the second elements (now all non-negative)
        total = sum(t[1] for t in electron_dist_non_negative)

        # Create new list with percentages
        new_electron_dist = [
            (*t, (t[1] / total) * 100 if total > 0 else 0)
            for t in electron_dist_non_negative
        ]
        
        electron_str = ", ".join(f"{a}: {perc_dist:.3f}%" for a, count, perc_dist in new_electron_dist if perc_dist > 0.01)
        print(f"{s:^8} | {eval_energy:^11.7f} | {electron_str}")