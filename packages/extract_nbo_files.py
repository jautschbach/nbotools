# import os
# import numpy as np
# import pandas as pd
# import time
# from multiprocessing import Pool
# import re
# import warnings

# def extract_floats_numpy(content, start_index, count):
#     return np.fromstring(content[start_index:], sep=' ', count=count)

# def create_symmetric_matrix_vectorized(lower_triangular, n):
#     matrix = np.zeros((n, n))
#     matrix[np.tril_indices(n)] = lower_triangular
#     return matrix + matrix.T - np.diag(matrix.diagonal())






# def process_47_file(file_path: str, nbas: int) -> tuple[bool, dict[str, np.ndarray]]:
    
   
#     # Read file content and detect calculation type
#     with open(file_path, 'r') as file:
#         content = file.read()
#     lines = content.split('\n')
#     is_open_shell = 'OPEN' in lines[0].upper()
    
#     # Calculate triangular matrix element count
#     float_count = int(nbas * (nbas + 1) / 2)
#     target_keywords = ['$OVERLAP', '$DENSITY', '$FOCK']
#     keyword_dict = {}

#     for keyword in target_keywords:
#         try:
#             # Locate keyword section
#             start_idx = content.index(keyword) + len(keyword)
            
#             if keyword in ['$DENSITY', '$FOCK'] and is_open_shell:
#                 # Process spin-polarized matrices
#                 spin_data = extract_floats_numpy(content, start_idx, 2 * float_count)
                
#                 alpha = create_symmetric_matrix_vectorized(spin_data[:float_count], nbas)
#                 beta = create_symmetric_matrix_vectorized(spin_data[float_count:], nbas)
                
#                 keyword_dict[f"{keyword[1:]}_ALPHA"] = alpha
#                 keyword_dict[f"{keyword[1:]}_BETA"] = beta
#             else:
#                 # Process closed-shell/single matrices
#                 data = extract_floats_numpy(content, start_idx, float_count)
#                 matrix = create_symmetric_matrix_vectorized(data, nbas)
#                 dict_key = keyword[1:] if keyword.startswith('$') else keyword
#                 keyword_dict[dict_key] = matrix
                
#         except ValueError:
#             # Handle missing keywords by creating zero matrices and printing warnings
#             zero_data = np.zeros(float_count)
#             if keyword in ['$DENSITY', '$FOCK'] and is_open_shell:
#                 keyword_dict[f"{keyword[1:]}_ALPHA"] = create_symmetric_matrix_vectorized(zero_data, nbas)
#                 keyword_dict[f"{keyword[1:]}_BETA"] = create_symmetric_matrix_vectorized(zero_data, nbas)
#             else:
#                 dict_key = keyword[1:] if keyword.startswith('$') else keyword
#                 keyword_dict[dict_key] = create_symmetric_matrix_vectorized(zero_data, nbas)

#             # Print warning message for missing keywords
#             print(f"\n{keyword.replace('$', '')} matrix is not available. Some analysis might not work.")

#     return is_open_shell, keyword_dict





# def extract_floats_from_line(line):
#     return [float(x) for x in line.split() if is_float(x)]

# def is_float(string):
#     try:
#         float(string)
#         return True
#     except ValueError:
#         return False



# def process_other_file(file_path, nbas, is_open_shell):
    
#     def extract_closed_shell():
#         for i, line in enumerate(lines):
#             if "in the" in line or "matrix" in line.lower():
#                 current_key = line.strip()
#                 is_matrix = "matrix" in line.lower()
#                 total_floats = int(nbas * (nbas + 1) / 2) if is_matrix else nbas * nbas
                
#                 start_line = i + 2  # Two lines after the keyword
#                 matrix_data = []
                
#                 float_count = 0
#                 for j in range(start_line, len(lines)):
#                     floats = [float(x) for x in lines[j].split() if x]
#                     matrix_data.extend(floats)
#                     float_count += len(floats)
#                     if float_count >= total_floats:
#                         break
                
#                 matrix_data = matrix_data[:total_floats]
                
#                 current_key = current_key.replace("in","").replace("the", "").replace("basis:", "")                        
#                 c1 = current_key.strip().split()[0]
#                 c2 = current_key.strip().split()[1]                                                
#                 current_key = f"{c2}{c1}"
#                 current_key = current_key.replace("densityAO", "AODM")
#                 current_key = current_key.replace("overlapPNAO", "PNAO_overlap_matrix")


                
#                 if is_matrix:
#                     matrix = create_symmetric_matrix_vectorized(matrix_data, nbas)
#                 else:
#                     matrix = np.array(matrix_data).reshape(nbas, nbas)
#                     matrix = matrix.T
                
#                 keyword_dict[current_key] = matrix
    
#     with open(file_path, 'r') as file:
#         lines = file.readlines()

#     keyword_dict = {}
    
#     if is_open_shell:
#         i = 0
#         while i < len(lines):
#             line = lines[i]
#             if "in the" in line or "matrix" in line.lower():
#                 current_key = line.strip()
#                 is_matrix = "matrix" in line.lower()
                
#                 i += 1
#                 while i < len(lines):
#                     next_line = lines[i].strip()
                    
#                     if "ALPHA" in next_line.upper() or "BETA" in next_line.upper():
#                         spin = "ALPHA" if "ALPHA" in next_line.upper() else "BETA"
                        
#                         current_key = current_key.replace("in","").replace("the", "").replace("basis:", "")                        
#                         c1 = current_key.strip().split()[0]
#                         c2 = current_key.strip().split()[1]                                                
#                         # full_key = f"{c2}{c1}_{spin}.upper()"
#                         full_key = f"{c2}{c1}_{spin.upper()}"
#                         full_key = full_key.replace("densityAO", "AODM")

                        
                        
                        
#                         total_floats = int(nbas * (nbas + 1) / 2) if is_matrix else nbas * nbas
#                         extracted_numbers = []
                        
#                         while len(extracted_numbers) < total_floats and i < len(lines) - 1:
#                             i += 1
#                             extracted_numbers.extend(extract_floats_from_line(lines[i]))
                        
#                         if len(extracted_numbers) != total_floats:
#                             print(f"Warning: Mismatch in number of floats for {full_key}. Expected {total_floats}, got {len(extracted_numbers)}")
                        
#                         if is_matrix:
#                             matrix = create_symmetric_matrix_vectorized(extracted_numbers, nbas)
#                         else:
#                             matrix = np.array(extracted_numbers).reshape(nbas, nbas)
#                             matrix = matrix.T
                            
#                         keyword_dict[full_key] = matrix
                    
#                     elif "in the" in next_line or "matrix" in next_line.lower():
#                         i -= 1  # Go back one line to process this new key
#                         break
                    
#                     i += 1
            
#             i += 1

#     if not is_open_shell or (is_open_shell and file_path.endswith("32") or file_path.endswith("33")):
#         extract_closed_shell()

#     return keyword_dict



# def process_file_wrapper(args):
#     file_path, nbas, is_open_shell = args
#     try:
#         if file_path.endswith('.47'):
#             return file_path, process_47_file(file_path, nbas)
#         else:
#             return file_path, process_other_file(file_path, nbas, is_open_shell)
#     except Exception as e:
#         print(f"Error processing {file_path}: {str(e)}")
#         return file_path, None

# def find_and_process_files():
#     current_directory = os.getcwd()
#     print(f"Current working directory: {current_directory}")

#     total_start_time = time.time()
    
#     # Find the .47 file in the current directory
#     file_47 = next((f for f in os.listdir(current_directory) if f.endswith('.47')), None)
#     if not file_47:
#         print("Error: No .47 file found in the current directory.")
#         return

#     basename = os.path.splitext(file_47)[0]
#     # print(f"Processing files with basename: {basename}")

#     with open(file_47, 'r') as file:
#         content = file.read()
#     nbas = int(next(line.split('NBAS=')[1].split()[0] for line in content.split('\n') if 'NBAS=' in line))
#     is_open_shell = 'OPEN' in content.split('\n')[0].upper()

#     file_extensions = ['47', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '48', '49']

#     file_paths = [
#         os.path.join(current_directory, f"{basename}.{ext}")
#         for ext in file_extensions
#         if os.path.exists(os.path.join(current_directory, f"{basename}.{ext}"))
#     ]

#     with Pool() as pool:
#         results = pool.map(process_file_wrapper, [(fp, nbas, is_open_shell) for fp in file_paths])

#     all_data = {}
#     for file_path, result in results:
#         if result:
#             if file_path.endswith('.47'):
#                 _, data = result
#             else:
#                 data = result
#             all_data[os.path.basename(file_path)] = data

#     # Print results
#     print(f"\nShell type: {'Open shell' if is_open_shell else 'Closed shell'}")
#     variables_mat = {}
#     for filename, data in all_data.items():
#         # print(f"\nProcessed file: {filename}")
#         for key, matrix in data.items():
#             # print(f"\n{key}:")
#             # print(pd.DataFrame(matrix))
#             variables_mat[key] = matrix

#     total_processing_time = time.time() - total_start_time
#     print(f"\nTotal files processed: {len(file_paths)}")
#     # print(f"Total processing time: {total_processing_time:.4f} seconds")
    
#     return variables_mat, basename, nbas



# def main():
#     return find_and_process_files()

# if __name__ == '__main__':
#     matrix_dict, filename, NBAS = main()
#     variables = matrix_dict
# else:
#     matrix_dict, filename, NBAS, variables = None, None, None, None

# # Allow other scripts to run the processing
# def run_processing():
#     return main()



import os
import numpy as np
from multiprocessing import Pool


def extract_floats_numpy(content, start_index, count):
    return np.fromstring(content[start_index:], sep=' ', count=count)

def create_symmetric_matrix_vectorized(lower_triangular, n):
    matrix = np.zeros((n, n))
    matrix[np.tril_indices(n)] = lower_triangular
    return matrix + matrix.T - np.diag(matrix.diagonal())

def process_47_file(file_path: str, nbas: int) -> tuple[bool, dict[str, np.ndarray]]:
    # Read file content and detect calculation type
    with open(file_path, 'r') as file:
        content = file.read()
    lines = content.split('\n')
    is_open_shell = 'OPEN' in lines[0].upper()
    
    # Calculate triangular matrix element count
    float_count = int(nbas * (nbas + 1) / 2)
    target_keywords = ['$OVERLAP', '$DENSITY', '$FOCK']
    keyword_dict = {}

    for keyword in target_keywords:
        try:
            # Locate keyword section
            start_idx = content.index(keyword) + len(keyword)
            
            if keyword in ['$DENSITY', '$FOCK'] and is_open_shell:
                # Process spin-polarized matrices
                spin_data = extract_floats_numpy(content, start_idx, 2 * float_count)
                
                alpha = create_symmetric_matrix_vectorized(spin_data[:float_count], nbas)
                beta = create_symmetric_matrix_vectorized(spin_data[float_count:], nbas)
                
                keyword_dict[f"{keyword[1:]}_ALPHA"] = alpha
                keyword_dict[f"{keyword[1:]}_BETA"] = beta
            else:
                # Process closed-shell/single matrices
                data = extract_floats_numpy(content, start_idx, float_count)
                matrix = create_symmetric_matrix_vectorized(data, nbas)
                dict_key = keyword[1:] if keyword.startswith('$') else keyword
                keyword_dict[dict_key] = matrix
                
        except ValueError:
            # Handle missing keywords by creating zero matrices and printing warnings
            zero_data = np.zeros(float_count)
            if keyword in ['$DENSITY', '$FOCK'] and is_open_shell:
                keyword_dict[f"{keyword[1:]}_ALPHA"] = create_symmetric_matrix_vectorized(zero_data, nbas)
                keyword_dict[f"{keyword[1:]}_BETA"] = create_symmetric_matrix_vectorized(zero_data, nbas)
            else:
                dict_key = keyword[1:] if keyword.startswith('$') else keyword
                keyword_dict[dict_key] = create_symmetric_matrix_vectorized(zero_data, nbas)

            print(f"\n{keyword.replace('$', '')} matrix is not available. Some analysis might not work.")

    return is_open_shell, keyword_dict

def extract_floats_from_line(line):
    return [float(x) for x in line.split() if is_float(x)]

def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def process_other_file(file_path, nbas, is_open_shell):
    def extract_closed_shell():
        for i, line in enumerate(lines):
            if "in the" in line or "matrix" in line.lower():
                current_key = line.strip()
                is_matrix = "matrix" in line.lower()
                total_floats = int(nbas * (nbas + 1) / 2) if is_matrix else nbas * nbas
                
                start_line = i + 2  # Two lines after the keyword
                matrix_data = []
                
                float_count = 0
                for j in range(start_line, len(lines)):
                    floats = [float(x) for x in lines[j].split() if x]
                    matrix_data.extend(floats)
                    float_count += len(floats)
                    if float_count >= total_floats:
                        break
                
                matrix_data = matrix_data[:total_floats]
                
                current_key = current_key.replace("in","").replace("the", "").replace("basis:", "")                        
                c1 = current_key.strip().split()[0]
                c2 = current_key.strip().split()[1]                                                
                current_key = f"{c2}{c1}"
                current_key = current_key.replace("densityAO", "AODM")
                current_key = current_key.replace("overlapPNAO", "PNAO_overlap_matrix")

                if is_matrix:
                    matrix = create_symmetric_matrix_vectorized(matrix_data, nbas)
                else:
                    matrix = np.array(matrix_data).reshape(nbas, nbas)
                    matrix = matrix.T
                
                keyword_dict[current_key] = matrix
    
    with open(file_path, 'r') as file:
        lines = file.readlines()

    keyword_dict = {}
    
    if is_open_shell:
        i = 0
        while i < len(lines):
            line = lines[i]
            if "in the" in line or "matrix" in line.lower():
                current_key = line.strip()
                is_matrix = "matrix" in line.lower()
                
                i += 1
                while i < len(lines):
                    next_line = lines[i].strip()
                    
                    if "ALPHA" in next_line.upper() or "BETA" in next_line.upper():
                        spin = "ALPHA" if "ALPHA" in next_line.upper() else "BETA"
                        
                        current_key = current_key.replace("in","").replace("the", "").replace("basis:", "")                        
                        c1 = current_key.strip().split()[0]
                        c2 = current_key.strip().split()[1]                                                
                        full_key = f"{c2}{c1}_{spin.upper()}"
                        full_key = full_key.replace("densityAO", "AODM")

                        total_floats = int(nbas * (nbas + 1) / 2) if is_matrix else nbas * nbas
                        extracted_numbers = []
                        
                        while len(extracted_numbers) < total_floats and i < len(lines) - 1:
                            i += 1
                            extracted_numbers.extend(extract_floats_from_line(lines[i]))
                        
                        if len(extracted_numbers) != total_floats:
                            print(f"Warning: Mismatch in number of floats for {full_key}. Expected {total_floats}, got {len(extracted_numbers)}")
                        
                        if is_matrix:
                            matrix = create_symmetric_matrix_vectorized(extracted_numbers, nbas)
                        else:
                            matrix = np.array(extracted_numbers).reshape(nbas, nbas)
                            matrix = matrix.T
                            
                        keyword_dict[full_key] = matrix
                    
                    elif "in the" in next_line or "matrix" in next_line.lower():
                        i -= 1  # Go back one line to process this new key
                        break
                    
                    i += 1
            
            i += 1

    if not is_open_shell or (is_open_shell and file_path.endswith("32") or file_path.endswith("33")):
        extract_closed_shell()

    return keyword_dict

def process_file_wrapper(args):
    file_path, nbas, is_open_shell = args
    try:
        if file_path.endswith('.47'):
            return file_path, process_47_file(file_path, nbas)
        else:
            return file_path, process_other_file(file_path, nbas, is_open_shell)
    except Exception as e:
        print(f"Error processing {file_path}: {str(e)}")
        return file_path, None

def find_and_process_files():
    current_directory = os.getcwd()
    print(f"Current working directory: {current_directory}")

    # total_start_time = time.time()
    
    # Find the .47 file in the current directory
    file_47 = next((f for f in os.listdir(current_directory) if f.endswith('.47')), None)
    if not file_47:
        print("Error: No .47 file found in the current directory.")
        return

    basename = os.path.splitext(file_47)[0]
    with open(file_47, 'r') as file:
        content = file.read()
    nbas = int(next(line.split('NBAS=')[1].split()[0] for line in content.split('\n') if 'NBAS=' in line))
    is_open_shell = 'OPEN' in content.split('\n')[0].upper()

    file_extensions = ['47', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '48', '49']

    file_paths = [
        os.path.join(current_directory, f"{basename}.{ext}")
        for ext in file_extensions
        if os.path.exists(os.path.join(current_directory, f"{basename}.{ext}"))
    ]

    with Pool() as pool:
        results = pool.map(process_file_wrapper, [(fp, nbas, is_open_shell) for fp in file_paths])

    all_data = {}
    for file_path, result in results:
        if result:
            if file_path.endswith('.47'):
                _, data = result
            else:
                data = result
            all_data[os.path.basename(file_path)] = data

    print(f"\nShell type: {'Open shell' if is_open_shell else 'Closed shell'}")
    variables_mat = {}
    for filename, data in all_data.items():
        for key, matrix in data.items():
            variables_mat[key] = matrix

    # total_processing_time = time.time() - total_start_time
    print(f"\nTotal files processed: {len(file_paths)}")
    
    return variables_mat, basename, nbas

def main():
    return find_and_process_files()

if __name__ == '__main__':
    matrix_dict, filename, NBAS = main()
else:
    matrix_dict, filename, NBAS = None, None, None

# Allow other scripts to run the processing
def run_processing():
    return main()