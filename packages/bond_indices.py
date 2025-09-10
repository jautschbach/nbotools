# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 10:42:02 2024

@author: aoben
"""

import numpy as np
from typing import Dict, Tuple, List
from dataclasses import dataclass

from basis_info import atom_data
from sep_mat import is_open as isOpenShell, matrix_dict, NBAS, ao_basis_list_dict
from generate_promolecule import alp_promol_den_mat, bet_promol_den_mat



@dataclass
class BondIndices:
    mayer: Dict[Tuple[str, str], float]
    wiberg: Dict[Tuple[str, str], float]
    nalewajski_mrozek: Dict[Tuple[str, str], float]

class BondOrderAnalyzer:
    def __init__(self):
        self.atom_list = self._generate_atom_list()
        self.atom_basis_ind_tup = self._generate_atom_basis_ind_tup()

    def _generate_atom_list(self) -> List[str]:
       
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

        atomic_number_to_symbol = {v: k for k, v in atom_dict.items()}
        return [atomic_number_to_symbol[atom[0]] for atom in atom_data]

    def _generate_atom_basis_ind_tup(self) -> Dict[str, Tuple[int, ...]]:
        return {
            f"{atom}_{i}": tuple(range(basis['bflo'], basis['bfhi'] + 1))
            for i, (atom, basis) in enumerate(zip(self.atom_list, ao_basis_list_dict), start=1)
        }

    def get_matrix(self, key: str) -> np.ndarray:
        return matrix_dict.get(key, np.zeros((NBAS, NBAS)))

    def is_zero_matrix(self, matrix: np.ndarray, precision: float = 1e-9) -> bool:
        return matrix.size == 0 or np.all(np.abs(matrix) < precision)

    def wiberg_bond(self, density: np.ndarray) -> Dict[Tuple[str, str], float]:
        wiberg_bond_ind = {}
        atoms = list(self.atom_basis_ind_tup.keys())
        for i, atom1 in enumerate(atoms):
            for atom2 in atoms[i+1:]:
                indices1, indices2 = self.atom_basis_ind_tup[atom1], self.atom_basis_ind_tup[atom2]
                double_sum = sum(density[a,b]**2 for a in indices1 for b in indices2)
                wiberg_bond_ind[(atom1, atom2)] = 2 * double_sum if isOpenShell else double_sum
        return wiberg_bond_ind

    def mayer_bond(self, density: np.ndarray, overlap: np.ndarray) -> Dict[Tuple[str, str], float]:
        PS = density @ overlap
        mayer_bond_ind = {}
        atoms = list(self.atom_basis_ind_tup.keys())
        for i, atom1 in enumerate(atoms):
            for atom2 in atoms[i+1:]:
                indices1, indices2 = self.atom_basis_ind_tup[atom1], self.atom_basis_ind_tup[atom2]
                double_sum = sum(PS[a,b] * PS[b,a] for a in indices1 for b in indices2)
                mayer_bond_ind[(atom1, atom2)] = double_sum
        return mayer_bond_ind

    def run_analysis(self) -> BondIndices:
        if isOpenShell:
            return self._analyze_open_shell()
        else:
            return self._analyze_closed_shell()

    def _analyze_open_shell(self) -> BondIndices:
        overlap_mat = self.get_matrix("OVERLAP")
        BODM_ALPHA, BODM_BETA = self.get_matrix("DENSITY_ALPHA"), self.get_matrix("DENSITY_BETA")
        SAO = np.array(overlap_mat)

        mayer_bo = self._calculate_mayer_bo_open_shell(BODM_ALPHA, BODM_BETA, SAO)
        wiberg_bo = self._calculate_wiberg_bo_open_shell(SAO, BODM_ALPHA, BODM_BETA)
        nm_bo = self._calculate_nm_bo_open_shell(SAO, BODM_ALPHA, BODM_BETA)

        return BondIndices(mayer_bo, wiberg_bo, nm_bo)

    def _analyze_closed_shell(self) -> BondIndices:
        overlap_mat = self.get_matrix("OVERLAP")
        SAO = np.array(overlap_mat)
        BODM = self.get_matrix("DENSITY")

        mayer_bo = self.mayer_bond(BODM, SAO)
        wiberg_bo = self._calculate_wiberg_bo_closed_shell(SAO, BODM)
        nm_bo = self._calculate_nm_bo_closed_shell(SAO, BODM)

        return BondIndices(mayer_bo, wiberg_bo, nm_bo)

    def _calculate_mayer_bo_open_shell(self, BODM_ALPHA: np.ndarray, BODM_BETA: np.ndarray, SAO: np.ndarray) -> Dict[Tuple[str, str], float]:
        mayer_bo_alpha = self.mayer_bond(BODM_ALPHA, SAO)
        mayer_bo_beta = self.mayer_bond(BODM_BETA, SAO)
        return {key: 2 * (value1 + value2) for (key, value1), (_, value2) in zip(mayer_bo_alpha.items(), mayer_bo_beta.items())}





    def _calculate_wiberg_bo_open_shell(self, SAO: np.ndarray, BODM_ALPHA: np.ndarray, BODM_BETA: np.ndarray) -> Dict[Tuple[str, str], float]:
        NAO = self.get_matrix("AONAOs")
        AODM_ALPHA, AODM_BETA = SAO @ BODM_ALPHA @ SAO, SAO @ BODM_BETA @ SAO
        DMNAO_alpha, DMNAO_beta = NAO.T @ AODM_ALPHA @ NAO, NAO.T @ AODM_BETA @ NAO
        wiberg_bo_alpha = self.wiberg_bond(DMNAO_alpha)
        wiberg_bo_beta = self.wiberg_bond(DMNAO_beta)
        return {key: value1 + value2 for (key, value1), (_, value2) in zip(wiberg_bo_alpha.items(), wiberg_bo_beta.items())}

    def _calculate_wiberg_bo_closed_shell(self, SAO: np.ndarray, BODM: np.ndarray) -> Dict[Tuple[str, str], float]:
        AODM = self.get_matrix('AODM')
        NAO = self.get_matrix('AONAOs')
        DMNAO = NAO.T @ AODM @ NAO
        return self.wiberg_bond(DMNAO)
    
    
    def _calculate_nm_bo_open_shell(self, SAO: np.ndarray, BODM_ALPHA: np.ndarray, BODM_BETA: np.ndarray) -> Dict[Tuple[str, str], float]:
        NAO = self.get_matrix("AONAOs")
        AODM_ALPHA, AODM_BETA = SAO @ BODM_ALPHA @ SAO, SAO @ BODM_BETA @ SAO
        DMNAO_alpha, DMNAO_beta = NAO.T @ AODM_ALPHA @ NAO, NAO.T @ AODM_BETA @ NAO
        
        alpha_den_one_cent_ion, beta_den_one_cent_ion = DMNAO_alpha, DMNAO_beta
        
        mole_alpha_ele = np.trace(alpha_den_one_cent_ion)
        mole_beta_ele = np.trace(beta_den_one_cent_ion)        
        mole_elec_cnt = np.round((mole_alpha_ele + mole_beta_ele),4)
        
        
        promole_alpha_ele = np.trace(alp_promol_den_mat)
        promole_beta_ele = np.trace(bet_promol_den_mat)        
        promole_elec_cnt = np.round((promole_alpha_ele + promole_beta_ele),4)
        
        
        # if promole_elec_cnt != mole_elec_cnt:
        #     delta_P_alpha = alpha_den_one_cent_ion - alp_promol_den_mat * mole_alpha_ele / promole_alpha_ele
        #     delta_P_beta = beta_den_one_cent_ion - bet_promol_den_mat * mole_beta_ele / promole_beta_ele
        #     # alp_promol_den_mat = * mole_alpha_ele / promole_alpha_ele
        #     # bet_promol_den_mat = bet_promol_den_mat * mole_beta_ele / promole_beta_ele
        
        # else:
        #     delta_P_alpha = alpha_den_one_cent_ion - alp_promol_den_mat
        #     delta_P_beta = beta_den_one_cent_ion - bet_promol_den_mat
            
        delta_P_alpha = alpha_den_one_cent_ion - alp_promol_den_mat
        delta_P_beta = beta_den_one_cent_ion - bet_promol_den_mat
        
        print(f"delta_P_alpha: {np.trace(delta_P_alpha):.4f}, delta_P_beta: {np.trace(delta_P_beta):.4f}, Total delta_P: {np.trace(delta_P_alpha + delta_P_beta):.4f}")            
        # print("Molecule: ", mole_elec_cnt, 'Promolecule: ', promole_elec_cnt)

            
        ion_one_cent_cont_s3 = self._nm_BI_sets_one_center_ion(alpha_den_one_cent_ion, beta_den_one_cent_ion, delta_P_alpha, delta_P_beta)
        cov_one_cent_cont_s3 = self._nm_BI_sets_one_center_cov(alpha_den_one_cent_ion, beta_den_one_cent_ion)
        two_cent_cont_s3, sum_two_cent_cont_s3 = self._nm_BI_sets_two_center(DMNAO_alpha, DMNAO_beta)
        
        weighing_factors_set_3 = self._calculate_weighing_factors(two_cent_cont_s3, sum_two_cent_cont_s3)
        total_one_center_s3 = {atom: ion_one_cent_cont_s3[atom] + cov_one_cent_cont_s3[atom] for atom in cov_one_cent_cont_s3}
        bond_indices_s3 = self._calculate_bond_indices(two_cent_cont_s3, total_one_center_s3, weighing_factors_set_3)
        
        return bond_indices_s3
    
    def _calculate_nm_bo_closed_shell(self, SAO: np.ndarray, BODM: np.ndarray) -> Dict[Tuple[str, str], float]:
        NAO = self.get_matrix("AONAOs")
        AODM = SAO @ BODM @ SAO
        DMNAO = NAO.T @ AODM @ NAO
        
        DMNAO_alpha_close = DMNAO
        DMNAO_beta_close = DMNAO_alpha_close * 0
        
        alpha_den_one_cent_ion = beta_den_one_cent_ion = DMNAO * 0.5
        delta_P_alpha = alpha_den_one_cent_ion - alp_promol_den_mat
        delta_P_beta = beta_den_one_cent_ion - bet_promol_den_mat
        
        ion_one_cent_cont_s3 = self._nm_BI_sets_one_center_ion(alpha_den_one_cent_ion, beta_den_one_cent_ion, delta_P_alpha, delta_P_beta)
        cov_one_cent_cont_s3 = self._nm_BI_sets_one_center_cov(DMNAO_alpha_close,  DMNAO_beta_close)
        two_cent_cont_s3, sum_two_cent_cont_s3 = self._nm_BI_sets_two_center(DMNAO_alpha_close,  DMNAO_beta_close)
        
        weighing_factors_set_3 = self._calculate_weighing_factors(two_cent_cont_s3, sum_two_cent_cont_s3)
        total_one_center_s3 = {atom: ion_one_cent_cont_s3[atom] + cov_one_cent_cont_s3[atom] for atom in cov_one_cent_cont_s3}
        bond_indices_s3 = self._calculate_bond_indices(two_cent_cont_s3, total_one_center_s3, weighing_factors_set_3)
        
        return bond_indices_s3
    
    def _nm_BI_sets_one_center_ion(self, alpha_den_one_cent_ion: np.ndarray, beta_den_one_cent_ion: np.ndarray, 
                                   delta_P_alpha: np.ndarray, delta_P_beta: np.ndarray) -> Dict[str, float]:
        one_cent_cont_s3 = {}
        for atom, indices in self.atom_basis_ind_tup.items():
            ion_cont_set3 = 0           
            for a in indices:
                ion_cont_set3 += (alpha_den_one_cent_ion[a, a] * delta_P_alpha[a, a] + 
                                  beta_den_one_cent_ion[a, a] * delta_P_beta[a, a])
            one_cent_cont_s3[atom] = ion_cont_set3
            
        
        print("\nOne Center Ionic Contribution for each atom")
        # print(one_cent_cont_s3)
        print({k: f'{v:.5f}' for k, v in one_cent_cont_s3.items()})

        return one_cent_cont_s3
    
    def _nm_BI_sets_one_center_cov(self, alpha_den_one_cent_ion: np.ndarray, beta_den_one_cent_ion: np.ndarray) -> Dict[str, float]:
        one_cent_cont_s3 = {}

        for atom, indices in self.atom_basis_ind_tup.items():
            cov_cont_set3 = 0
            for a in indices:
                for a_prime in indices:
                    if a < a_prime:  # Condition to exclude diagonal elements
                        if isOpenShell:
                            cov_cont_set3 += (alpha_den_one_cent_ion[a, a_prime]**2 + 
                                              beta_den_one_cent_ion[a, a_prime]**2)
                        else:
                            cov_cont_set3 += (alpha_den_one_cent_ion[a, a_prime]**2 + 
                                              beta_den_one_cent_ion[a, a_prime]**2)
            if isOpenShell:
                cov_cont_set3 *= 2
            one_cent_cont_s3[atom] = cov_cont_set3
            
        print("\nOne Center Covalent Contribution for each atom")
        # print(one_cent_cont_s3)
        print({k: f'{v:.5f}' for k, v in one_cent_cont_s3.items()})
        return one_cent_cont_s3            
            

    
    def _nm_BI_sets_two_center(self, DMNAO_alpha: np.ndarray, DMNAO_beta: np.ndarray) -> Tuple[Dict[Tuple[str, str], float], Dict[str, float]]:
        two_cent_cont_s3 = {}

        sum_two_cent_cont_s3 = {atom: 0 for atom in self.atom_basis_ind_tup}
        atoms = list(self.atom_basis_ind_tup.keys())
    
        for i, atom1 in enumerate(atoms):
            indices1 = self.atom_basis_ind_tup[atom1]
            for atom2 in atoms[i+1:]:
                indices2 = self.atom_basis_ind_tup[atom2]
                                             
                cov_cont_set3 = 0
                for a in indices1:
                    for b in indices2:
                        if isOpenShell:
                            cov_cont_set3 += (DMNAO_alpha[a, b]**2 + DMNAO_beta[a, b]**2)
                        else:
                            cov_cont_set3 += ( DMNAO_alpha[a, b]**2 + DMNAO_beta[a, b]**2)
                
                if isOpenShell:
                    cov_cont_set3 *= 2
                
                two_cent_cont_s3[(atom1, atom2)] = cov_cont_set3
                sum_two_cent_cont_s3[atom1] += cov_cont_set3
                sum_two_cent_cont_s3[atom2] += cov_cont_set3
    
        print_significant_bonds(two_cent_cont_s3, "G-J bond Indices", 0.2)
        
                
                
    
        return two_cent_cont_s3, sum_two_cent_cont_s3
    
    def _calculate_weighing_factors(self, two_center_contributions: Dict[Tuple[str, str], float], 
                                    total_contributions: Dict[str, float]) -> Dict[Tuple[str, str], Dict[str, float]]:
        weighing_fact = {
            (atom1, atom2): {
                atom1: contribution / total_contributions[atom1],
                atom2: contribution / total_contributions[atom2]
            }
            for (atom1, atom2), contribution in two_center_contributions.items()
        }
        
        print_significant_weighting_factors(weighing_fact,  0.2)
        return weighing_fact
    
    def _calculate_bond_indices(self, two_center_contributions: Dict[Tuple[str, str], float], 
                                total_one_center_contributions: Dict[str, float], 
                                weighing_factors: Dict[Tuple[str, str], Dict[str, float]]) -> Dict[Tuple[str, str], float]:
        return {
            (atom1, atom2): V_AB_2 + weighing_factors[(atom1, atom2)][atom1] * total_one_center_contributions[atom1] + 
                            weighing_factors[(atom1, atom2)][atom2] * total_one_center_contributions[atom2]
            for (atom1, atom2), V_AB_2 in two_center_contributions.items()
        }    

    

def print_significant_bonds(bond_pairs: Dict[Tuple[str, str], float], title: str, threshold: float = 0.1):
    print(f"\nSignificant {title} (absolute value > {threshold}):")
    print("-" * 50)
    print(f"{'Atom Pair':<20}{'Bond Index':>15}")
    print("-" * 50)
    for (atom1, atom2), value in bond_pairs.items():
        if abs(value) > threshold:
            atom_pair = f"{atom1:<4} - {atom2:<4}"
            print(f"{atom_pair:<20}{value:15.4f}")
    print("-" * 50)

def print_significant_weighting_factors(weighing_factors: Dict[Tuple[str, str], Dict[str, float]], threshold: float = 0.00002):
    print(f"\nSignificant weighting factors (at least one value > {threshold}):")
    print("-" * 65)
    print(f"{'Atom1':<10}{'Weight1':<15}{'Atom2':<10}{'Weight2':<15}")
    print("-" * 65)
    for (atom1, atom2), weights in weighing_factors.items():
        if any(abs(weight) > threshold for weight in weights.values()):
            weight1, weight2 = weights[atom1], weights[atom2]
            print(f"{atom1:<10}{weight1:<15.4f}{atom2:<10}{weight2:<15.4f}")
    print("-" * 65)

def bnd_run():
    analyzer = BondOrderAnalyzer()
    results = analyzer.run_analysis()
    
    print_significant_bonds(results.nalewajski_mrozek, "Nalewajski-Mrozek Bond Indices", threshold=0.2)
    print_significant_bonds(results.mayer, "Mayer Bond Indices", threshold=0.2)
    print_significant_bonds(results.wiberg, "Wiberg Bond Indices", threshold=0.2)
    
    

# if __name__ == "__main__":
#     bnd_run()