import numpy as np
import pandas as pd


def put_elements_greater_than_val_in_qoutes(df, val):
    # Set pandas display options to print the full matrix
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.expand_frame_repr', False)

    df = df.map(lambda x: "'" + str(x) + "'" if abs(x) > val else x)
    print(df)
    return df

def display_available_keys(matrix_dict):
    keys = list(matrix_dict.keys())
    key_dict = {i + 1: key for i, key in enumerate(keys)}

    print("Available keys:")
    for num, key in key_dict.items():
        print(f"{num}. {key}")
    print("-1. Back")


def get_user_input(prompt):
    user_input = input(prompt)
    return user_input


def validate_key_selection(user_input, key_dict):
    if user_input.isdigit():
        num = int(user_input)
        if num in key_dict:
            return key_dict[num]
        else:
            return None
    else:
        return None


def get_range_input(prompt, min_value, max_value):
    user_input = int(get_user_input(prompt))

    while user_input < min_value or user_input > max_value:
        print(f"Invalid input. Please enter a number between {min_value} and {max_value}.")
        user_input = int(get_user_input(prompt))

    return user_input


def print_full_matrix(matrix):
    full = pd.DataFrame(matrix)
    full.index = full.index + 1
    full.columns = full.columns + 1
    # put_elements_greater_than_val_in_qoutes(full, 0.01)
    print(full)
    # transform_dataframe(full, 0.01)
    # display_dataframe(full)
    return full


def print_selected_range(matrix, start_row, end_row, start_col, end_col):
    reduced_matrix = matrix.loc[start_row:end_row, start_col:end_col]
    reduced_matrix.index = reduced_matrix.index + 1
    reduced_matrix.columns = reduced_matrix.columns + 1
    print(reduced_matrix)
    return reduced_matrix


def save_matrix(matrix, selected_key, start_row=None, end_row=None):
    filename = f"{selected_key}"
    if start_row is not None and end_row is not None:
        filename += f"_{start_row}-{end_row}"
    filename += "_matrix.out"
    matrix.to_csv(filename, index=True, encoding='utf-8')
    print("Result saved successfully.")

def parse_index(index):
    if "-" in index:
        start, end = map(int, index.split("-"))
        return list(range(start - 1, end ))
    else:
        return [int(index) - 1]

def submatrix(matrix, *args):
    columns = []
    rows = []
    for arg in args:
        if isinstance(arg, int):
            columns.extend(parse_index(str(arg).strip()))
            rows.extend(parse_index(str(arg).strip()))
        elif isinstance(arg, str):
            if "-" in arg:
                columns.extend(parse_index(arg.strip()))
                rows.extend(parse_index(arg.strip()))
            else:
                columns.extend(parse_index(str(arg).strip()))
                rows.extend(parse_index(str(arg).strip()))
                
    is_square = matrix.shape[0] == matrix.shape[1]
    if is_square:
        submatrix = matrix.iloc[rows, columns]

    else:
        submatrix = matrix.iloc[rows, :]

    return submatrix

def select_cols_only(matrix, *args):
    columns = []
    rows = []
    for arg in args:
        if isinstance(arg, int):
            columns.extend(parse_index(str(arg).strip()))
            rows.extend(parse_index(str(arg).strip()))
        elif isinstance(arg, str):
            if "-" in arg:
                columns.extend(parse_index(arg.strip()))
                rows.extend(parse_index(arg.strip()))
            else:
                columns.extend(parse_index(str(arg).strip()))
                rows.extend(parse_index(str(arg).strip()))
                
    is_square = matrix.shape[0] == matrix.shape[1]
    if is_square:
        submatrix = matrix.iloc[:, columns]

    else:
        # print(rows, columns)
        submatrix = matrix.iloc[:, columns]

    return submatrix



def view_matrix(matrix):
    user_input = get_user_input("Print full or select range? (f/s): ")

    if user_input.lower() == 'f':
        full_matrix = print_full_matrix(matrix)
        # print(full_matrix)

        save_input = get_user_input("Do you want to save the result? (y/n): ")

        if save_input.lower() == 'y':
            file_name = matrix
            
            print(file_name)
            # save_matrix(full_matrix, file_name)
        else:
            print("Result not saved.")
        
        return full_matrix
        
    elif user_input.lower() == 's':
        # Accept user input for the index argument
        index_input = input("Enter the index argument: ")

        # Split the input by commas to get individual arguments
        index_arguments = index_input.split(",")
        full_matrix = pd.DataFrame(matrix)
        full_matrix.index = full_matrix.index + 1
        full_matrix.columns = full_matrix.columns + 1
        # print(full_matrix)
        reduced_matrix = submatrix(full_matrix, *index_arguments)
        print(reduced_matrix)
        # print("_jjjj________________________________")
        # print(type(reduced_matrix))
        
        return pd.DataFrame(reduced_matrix)
        
        save_input = get_user_input("Do you want to save the result? (y/n): ")

        if save_input.lower() == 'y':
            file_name = input("Enter file name without extension: ")
            file_name = file_name +".csv"
            # file_name = f"{matrix}_{index_arguments}"
            # file_name = file_name.replace("'","").replace("[", "").replace("]","")
            # file_name = file_name.replace(" ", "").replace(",", "_")
            
            save_matrix(reduced_matrix, file_name)
            print(file_name, " saved!!!")
        else:
            print("File not saved.")
       
    

def vi_save_mat_from_dict(matrix_dict):
    keys = list(matrix_dict.keys())
    key_dict = {i + 1: key for i, key in enumerate(keys)}

    while True:
        display_available_keys(matrix_dict)
        user_input = get_user_input("Enter the number corresponding to a key: ")

        if user_input == "-1":
            break

        selected_key = validate_key_selection(user_input, key_dict)

        if selected_key is None:
            print("Invalid input. Number does not correspond to a key.")
            continue

        print(f"Matrix for key '{selected_key}':")
        
        matrix_data = (matrix_dict[selected_key])
        # view_matrix(matrix_data)

        user_input = get_user_input("Print full or select range? (f/s): ")

        if user_input.lower() == 'f':
            full_matrix = print_full_matrix(matrix_dict[selected_key])
            # print(full_matrix)
           

            save_input = get_user_input("Do you want to save the result? (y/n): ")

            if save_input.lower() == 'y':
                # save_matrix(full_matrix, selected_key)
                file_name = f"{selected_key}.csv"
                full_matrix.to_csv(file_name, index=True)
                print(file_name, " saved!!!")
            else:
                print("Result not saved.")

        elif user_input.lower() == 's':
            # Accept user input for the index argument
            index_input = input("Enter the index argument: ")

            # Split the input by commas to get individual arguments
            index_arguments = index_input.split(",")
            full_matrix = pd.DataFrame(matrix_dict[selected_key])
            full_matrix.index = full_matrix.index + 1
            full_matrix.columns = full_matrix.columns + 1
            # print(full_matrix)
            reduced_matrix = submatrix(full_matrix, *index_arguments)
            print(reduced_matrix)
            # print(type(reduced_matrix))
            
           
            
            save_input = get_user_input("Do you want to save the result? (y/n): ")

            if save_input.lower() == 'y':
                file_name = input("Enter file name without extension: ")
                file_name = file_name +".csv"
                # file_name = f"{selected_key}_{index_arguments}.csv"
                # file_name = file_name.replace("'","").replace("[", "").replace("]","")
                # file_name = file_name.replace(" ", "").replace(",", "_")
                
                # save_matrix(reduced_matrix, file_name)
                reduced_matrix.to_csv(file_name, index=True)
                print(file_name, " saved!!!")
            else:
                print("File not saved.")
           



            
# PNAOPNHOs_1_3_5-6.csv
def get_matrix(matrix):
    user_input = get_user_input("Print full or select range? (f/s): ")

    if user_input.lower() == 'f':
        full_matrix = print_full_matrix(matrix)
        # transform_dataframe(full_matrix, 0.001)
        # print(full_matrix)
        
        return full_matrix

        
        
    elif user_input.lower() == 's':
        # Accept user input for the index argument
        index_input = input("Enter the index argument: ")

        # Split the input by commas to get individual arguments
        index_arguments = index_input.split(",")
        full_matrix = pd.DataFrame(matrix)
        full_matrix.index = full_matrix.index + 1
        full_matrix.columns = full_matrix.columns + 1
        # print(full_matrix)
        reduced_matrix = submatrix(full_matrix, *index_arguments)
        # transform_dataframe(reduced_matrix, 0.001)
        # print(reduced_matrix)
        # print(type(reduced_matrix))
        
        return pd.DataFrame(reduced_matrix)
    
def check_all_zero(df, tol=1e-7):
    result = np.isclose(df, 0, atol=tol).all().all()
    if result:
        print("All elements are zero.")
    else:
        print("Not all elements are zero.")
    return result


def is_unit_matrix(df, tol=1e-7):
    if not isinstance(df, pd.DataFrame):
        raise ValueError("Input must be a DataFrame.")

    n_rows, n_cols = df.shape

    if n_rows != n_cols:
        return False

    diagonal = df.values.diagonal()

    if not np.allclose(diagonal, 1, atol=tol):
        return False

    off_diagonal = df.values[~np.eye(n_rows, dtype=bool)]

    if not np.allclose(off_diagonal, 0, atol=tol):
        return False

    return True


