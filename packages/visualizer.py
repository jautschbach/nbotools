

import numpy as np
import pandas as pd
from sep_mat import  opera_matrix


def arrange_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Move all-zero rows to the end of the DataFrame.
    """
    return (
        df.assign(all_zero=(df == 0).all(axis=1))
          .sort_values("all_zero")
          .drop(columns="all_zero")
          .reset_index(drop=True)
    )


def diag_one_zero(df: pd.DataFrame, precision: float = 1e-6) -> bool:
    """
    Check if a DataFrame is:
    - square
    - diagonal
    - with only 0 or 1 allowed on the diagonal (within tolerance).
    """
    if df.shape[0] != df.shape:
        return False

    # Off-diagonal elements must be ~0
    off_diag_ok = np.allclose(df.values - np.diag(np.diag(df.values)), 0, atol=precision)

    # Diagonal elements must be ~0 or ~1
    diag = np.diag(df.values)
    diag_ok = np.all(
        [np.isclose(val, 0, atol=precision) or np.isclose(val, 1, atol=precision) for val in diag]
    )

    return off_diag_ok and diag_ok



def _write_matrix(file, matrix: pd.DataFrame, values_per_line: int = 5, pad_value: float = 1e-10):
    """
    Utility function:
    Write matrix values to open file stream with controlled formatting.
    If the matrix is not square, pad extra columns with `pad_value` (default 1e-20).
    """
    if matrix.shape[0] != matrix.shape:
        num_rows, num_cols = matrix.shape
        if num_cols < num_rows:
            # Add padding columns to make the matrix square
            new_cols = num_rows - num_cols
            pad_cols = pd.DataFrame(
                np.full((num_rows, new_cols), pad_value),  # use 1e-20 instead of zeros
                index=matrix.index
            )
            # Generate new column names
            max_col = max(matrix.columns, default=0)
            new_col_names = [max_col + i + 1 for i in range(new_cols)]
            pad_cols.columns = new_col_names
            matrix = pd.concat([matrix, pad_cols], axis=1)

    # Debug print (can remove later)
    print(matrix)

    # Write values
    for col in matrix.columns:
        column_data = matrix[col].to_numpy()
        for i, value in enumerate(column_data):
            file.write(f"{value:.13E}")
            if (i + 1) % values_per_line == 0:
                file.write("\n")
            else:
                file.write(" ")
        file.write("\n")


def save_loc_data(matrix: pd.DataFrame, file_name: str, text: str, matrix2: pd.DataFrame = None):
    """
    Save matrix/matrices to a formatted text file for NBO Tools.
    
    Parameters
    ----------
    matrix : pd.DataFrame
        Primary matrix (always saved).
    file_name : str
        Name of the output file.
    text : str
        Header text description to include in the output file.
    matrix2 : pd.DataFrame, optional
        If provided, writes as ALPHA (matrix) and BETA (matrix2).
        If not provided, writes just one block with no spin labels.
    """
    with open(file_name, "w") as file:
        file.write(
f"""NBO Tools: {file_name.replace('.40', '')}
{text}
-------------------------------------------------------------------------------
"""
        )

        if matrix2 is not None:  
            # Case: Spin-polarized (Alpha + Beta)
            file.write("ALPHA SPIN\n")
            _write_matrix(file, matrix)

            file.write("BETA  SPIN\n")
            _write_matrix(file, matrix2)

        else:
            # Case: No spin separation, just write one matrix
            _write_matrix(file, matrix)

    print(f"{file_name} saved...")


def sub_to_full(matrix_input: pd.DataFrame) -> pd.DataFrame:
    """
    Expand a submatrix into the full AO overlap matrix space.  
    After expansion, drop all-zero columns (orbitals not present).

    Parameters
    ----------
    matrix_input : pd.DataFrame
        Submatrix to be expanded.

    Returns
    -------
    pd.DataFrame
        Full matrix with submatrix embedded and zero-columns removeGGGd.
    """
    # Full overlap matrix (from opera_matrix dictionary)
    full_overlap = pd.DataFrame(opera_matrix["OVERLAP"])
    n = full_overlap.shape[0]

    # Prepare empty full-dimension DataFrame
    full_data = pd.DataFrame(np.zeros((n, n)), index=range(1, n + 1), columns=range(1, n + 1))

    # Ensure integer indices for the submatrix
    matrix = matrix_input.copy()
    matrix.index = matrix.columns = [int(c) for c in matrix.columns]

    # Insert submatrix into the larger blank matrix
    for row_index, row in matrix.iterrows():
        full_data.loc[row_index, row.index] = row.values

    # Drop all-zero columns (orbitals not present in subspace)
    full_data = full_data.loc[:, (full_data != 0).any(axis=0)]

    return full_data
