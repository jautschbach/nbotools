#!/usr/bin/env python3


from enum import Enum, auto
from typing import Callable, Dict

# import nlmo_analysis
# import val_core_ortho

class MenuOption(Enum):
    """Enumeration of menu options for better type safety and readability."""
    EXIT = 0
    TRANSFORMATION_ENERGIES = auto()
    ORBITAL_TRANSFORMATION = auto()
    OPERATOR_MATRICES = auto()
    PNAO_ORBITALS = auto()
    PNAO_ENERGIES = auto()
    # GEE_SOLVER = auto()
    ORBITAL_LOCALIZATION = auto()
    BOND_ORDER_ANALYSIS = auto()
    RHE_SOLVER = auto()

class NboToolsMenu:
    """A comprehensive menu-driven quantum chemistry analysis tool."""
    
    def __init__(self, actions: Dict[MenuOption, Callable]):
        """
        Initialize the menu with a dictionary of actions.
        
        Args:
            actions (Dict[MenuOption, Callable]): Mapping of menu options to their corresponding functions.
        """
        self._actions = actions
    
    def _display_menu(self) -> None:
        """Display the menu options with clear, descriptive text."""
        print("\n===== NBOTOOLS Menu =====")
        menu_descriptions = {
            MenuOption.TRANSFORMATION_ENERGIES: "Computes energies such as AOMO energies.",
            MenuOption.ORBITAL_TRANSFORMATION: "Performs transformations like AOMO",
            MenuOption.OPERATOR_MATRICES: "Retrieves matrices, such as the Fock matrix in the AO basis (if available)",
            MenuOption.PNAO_ORBITALS: "Computes orbitals in the PNAO basis (e.g., PNAOMO).",
            MenuOption.PNAO_ENERGIES: "Calculates orbital energies in the PNAO basis (e.g., PNAOMO energy)",
            # MenuOption.GEE_SOLVER: "Generalized Eigenvalue Equation Solver use option 9",
            MenuOption.ORBITAL_LOCALIZATION: "Performs Pipek-Mezey orbital localization",
            MenuOption.BOND_ORDER_ANALYSIS: "Bond order analysis based on NBO data",
            MenuOption.RHE_SOLVER: "Solves generalized eigenvalue equation",
            MenuOption.EXIT: "Terminates the program"
        }
        
        for option, description in menu_descriptions.items():
            print(f"{option.value}. {description}")
    
    def run(self) -> None:
        """
        Run the interactive menu-driven interface.
        Handles user interactions and executes selected actions.
        """
        while True:
            try:
                self._display_menu()
                choice = MenuOption(int(input("\nEnter your choice: ")))
                
                if choice == MenuOption.EXIT:
                    print("Exiting the program. Goodbye!")
                    break
                
                action = self._actions.get(choice)
                if action:
                    action()
                else:
                    print("Invalid option selected.")
            
            except ValueError:
                print("Please enter a valid numeric option.")
            except Exception as e:
                print(f"An error occurred: {e}")

def main():
    """Main function to set up and run Nbotools menu."""
    from my_func import vi_save_mat_from_dict
    from sep_mat import (
        aoX_trans, aoX_energies_dict, opera_matrix, 
        pnaoX_trans, pnaoX_energies_dict
    )
    # from program_solver import solve_gee
    from pm_localization import pm_run
    from bond_indices import bnd_run
    from val_core_ortho import solve_rhe

    # Define actions mapping
    actions = {
        MenuOption.TRANSFORMATION_ENERGIES: lambda: vi_save_mat_from_dict(aoX_energies_dict),
        MenuOption.ORBITAL_TRANSFORMATION: lambda: vi_save_mat_from_dict(aoX_trans),
        MenuOption.OPERATOR_MATRICES: lambda: vi_save_mat_from_dict(opera_matrix),
        MenuOption.PNAO_ORBITALS: lambda: vi_save_mat_from_dict(pnaoX_trans),
        MenuOption.PNAO_ENERGIES: lambda: vi_save_mat_from_dict(pnaoX_energies_dict),
        # MenuOption.GEE_SOLVER: solve_gee,
        MenuOption.ORBITAL_LOCALIZATION: pm_run,
        MenuOption.BOND_ORDER_ANALYSIS: bnd_run,
        MenuOption.RHE_SOLVER: solve_rhe
    }

    # Create and run the menu
    menu = NboToolsMenu(actions)
    menu.run()

if __name__ == "__main__":
    main()

# Record the end time
# end_time = time.time()
# 

# print(start_time, "--------", end_time)
# Calculate the elapsed time
# elapsed_time = end_time - start_time
# 
# Print the elapsed time
# print("Elapsed time: {:.2f} seconds".format(elapsed_time))

# Main program
