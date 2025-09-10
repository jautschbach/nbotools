import re
import os
import sys

def process_file(input_file_name, output_file_path):
    # function to check if a string can be converted to float
    def is_float(s):
        try:
            float(s)
            return True
        except ValueError:
            return False
    
    # function to format a line by adding commas after float numbers
    def format_line(line, next_line):
        float_pattern = r'(?<!\w)[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?(?!\w)'
        numbers = list(re.finditer(float_pattern, line))
        numbers.sort(key=lambda x: x.start(), reverse=True)
        for match in numbers:
            start, end = match.span()
            line = line[:end] + ',' + line[end:]
        # Remove trailing comma if it's the last number in the section
        if line.strip().endswith(',') and (not next_line or not is_float(next_line.strip().split()[0])):
            line = line.rstrip().rstrip(',') + '\n'
        return line

    try:
        with open(input_file_name, 'r') as input_file:
            lines = input_file.readlines()

        coord_sect = None
        coord_sect_end = None
        
        gennbo_sect = None
        # gennbo_sect_end = None
        
        basis_sect = None
        
        z_nucl_charge = []
        modified_lines = []
        
        job_name = None
        # overlap_sect = None

        for line_number, line in enumerate(lines):
            words = line.strip().split()
            
            # Skip lines containing "ORCA Job:" or Job name
            if "ORCA Job:" in line:
                # print(f"Line {line_number + 1}: Skipped (contains 'ORCA')")
                # print(line)
                # job_name = line
                
                pass
            
            if words:
                if words[0].startswith('$') and words[0] != "$END":
                    line = line.replace("$END", "/")
                    line = line.replace(words[0][0], "&")
                    # print(f"{words[0]} in Line {line_number + 1}: Starts with $")
                    
                    if words[0] == "$COORD":
                        line = line.replace("&COORD", "&molecule znuc =")
                        coord_sect = line_number
                        #Job name comes after $COORD"
                        job_name = lines[coord_sect + 1].strip()
                        
                    
                    elif words[0] == "$GENNBO":
                        gennbo_sect = line_number 
                        #Set the end of gennbo line to the start 
                        #Since they are always on the same line
                        # gennbo_sect_end = gennbo_sect 
                        
                    
                    elif words[0] == "$BASIS":
                        basis_sect = line_number 
                        
                    # elif words[0] == '$OVERLAP':
                    #     overlap_sect = line_number
                                                               
                    
                    modified_lines.append(line)                   
                    
                    
                
                elif words[-1] == "$END":
                    line = line.replace("$END", "/")
                    # print(f"Line {line_number + 1}: $END replaced with /")
                    modified_lines.append(line)
                    if coord_sect is not None and coord_sect_end is None:
                        coord_sect_end = line_number
                    # print("\n------------------")
                
                else:
                    next_line = lines[line_number + 1] if line_number + 1 < len(lines) else ""
                    formatted_line = format_line(line, next_line)
                    modified_lines.append(formatted_line)
                    
                    # Extract nuclear charges in coordinate section
                    if coord_sect is not None and coord_sect_end is None:
                        # if line_number == coord_sect + 1:
                        #     #Skip the 
                        #     continue
                        # else:                        
                        # z_nuc_val = re.findall(r'-?\d+(?:\.\d+)?', line)
                        #Do not use the above. Might break if job name is not included in file
                        
                        coord_eles = re.split(r'[,\s]+', line.strip())
                        
        
                        # Filter out non-numeric coord_eles
                        z_nuc_val = [coord_ele for coord_ele in coord_eles if re.match(r'^-?\d+(?:\.\d+)?$', coord_ele)]
                        
                        #The sec
                        if len(z_nuc_val) > 4:
                            z_nucl_charge.append(str(float(z_nuc_val[1])))
                        
            else:
                modified_lines.append(line)
                
  
        #Do some specific replacement of files
        modified_lines = [
            re.sub(r'\bEXP\b', 'ZETA',
            re.sub(r'\bCS\b', 'CCOE(:,0)',
            re.sub(r'\bCP\b', 'CCOE(:,1)',
            re.sub(r'\bCD\b', 'CCOE(:,2)',                   
            re.sub(r'\bCF\b', 'CCOE(:,3)',                   
            re.sub(r'\bCG\b', 'CCOE(:,4)',
            re.sub(r'\bCH\b', 'CCOE(:,5)',
            re.sub(r'\bCI\b', 'CCOE(:,6)', line))))))))
            for line in modified_lines ]   

        with open(output_file_path, 'w') as output_file:
            i = 0
            skip_next = False
            # Initialize a flag to indicate if we've passed the &OVERLAP line
            passed_overlap = False
            passed_contract = False
            while i < len(modified_lines):
                line = modified_lines[i]
        
                if skip_next:
                    skip_next = False
                    if line.strip() == '' or '/' in line:
                        i += 1
                        continue
        
                if i == gennbo_sect:
                    keywords = [("UPPER", "UPPER=T"), ("BOHR", "BOHR=T"), ("OPEN", "UNREST=T"), ("BODM", "BODM=T")]
                    format_keywords = ["FORMAT=PRECISE", "FORMAT"]
                    # Convert the entire line to uppercase before processing
                    line = line.upper()

                    for keyword, replacement in keywords:
                        if keyword in line:
                            line = line.replace(keyword, replacement)
                        else:
                            words = line.split()
                            words.insert(-1, replacement.split('=')[0] + "=F")
                            line = ' '.join(words) + '\n'
                    
                    # Handle FORMAT and FORMAT=PRECISE
                    if any(format_keyword in line for format_keyword in format_keywords):
                        for format_keyword in format_keywords:
                            line = line.replace(format_keyword, "FORM=T")
                    else:
                        words = line.split()
                        words.insert(-1, "FORM=F")
                        line = ' '.join(words) + '\n'
                    
                    output_file.write(line)                
                
                
                #Skip the &NBO section continue to the molecule section.
                elif '&NBO' in line:
                    skip_next = True
                                
                elif i == coord_sect:
                    output_file.write(f"&molecule\n jobtitle=\"{job_name}\"\n znuc ={','.join(z_nucl_charge)}\n")


                elif coord_sect < i < coord_sect_end:                  
                    if i == coord_sect + 1:
                        #Replace the Job title lin with coord=
                        mod_coord_sect = " coord="
                    else:
                        mod_coord_sect = ','.join(line.split(',')[2:]).strip() 
                        
                    if mod_coord_sect:
                        output_file.write(f"{' '*8 if i > coord_sect + 1 else ''}{mod_coord_sect}\n")
              
                

                # elif coord_sect < i < coord_sect_end:
                #     indent = ''
                #     mod_coord_sect = ''
                
                #     if i == coord_sect + 1:
                #         line_parts = line.split()
                #         #Avoid the title line by all means
                #         if len(line_parts) == 5 and all(part.replace('.', '').isdigit() for part in line_parts[2:]):
                #             coord_data = ' '.join(line_parts[2:])
                #             mod_coord_sect = f" coord={coord_data}"
                #         else:
                #             mod_coord_sect = " coord="
                #     else:
                #         mod_coord_sect = ','.join(line.split(',')[2:]).strip() + '\n'
                        
                #     if i > coord_sect + 1:
                #         indent = ' ' * 8
                #         #Handle the special case the second line only
                #         if i == coord_sect + 2:
                #             indent = ' ' * 1
                    
                #     if mod_coord_sect:
                #         output_file.write(f"{indent}{mod_coord_sect}")                 
                                        
                                    
                elif '&BASIS' in line:
                    output_file.write("&BASIS\n")
    
                elif '&CONTRACT' in line:
                    # Read the next two lines
                    nshell_line = modified_lines[i+1] if i+1 < len(modified_lines) else ""
                    nexp_line = modified_lines[i+2] if i+2 < len(modified_lines) else ""
    
                    # Write NSHELL and NEXP before &CONTRACT
                    if "NSHELL" in nshell_line:
                        output_file.write(nshell_line)
                    if "NEXP" in nexp_line:
                        output_file.write(nexp_line)
                    
                    output_file.write("/\n&CONTRACT\n")
    
                    # Skip the next two lines as we've already processed them
                    i += 2
                    
                    passed_contract = True
                    
            
                                                      
                # Don't write "/"  after this line automatically. They will be placed manually where necessary
                elif i > basis_sect and line.strip()[0] == "/":
                    pass  
                
    
                elif any(keyword in line for keyword in ['&OVERLAP', '&DENSITY', '&FOCK', '&LCAOMO', '&KINETIC', '&NUCLEAR']):
                    #Add an two lines before Overlap= since $overlap  is always in all file47 files
                    if '&OVERLAP' in line:
                        #This add's 
                        
                        
                        output_file.write("/\n&MATRICES\n")
                        passed_overlap = True  # Set the flag to True
                    keyword = line.split('&')[1].split()[0]
                    output_file.write(line.replace(f"&{keyword}", f"{keyword}="))
                    
                    
                # New elif block to handle other '&' keywords after &OVERLAP
                elif passed_overlap and line.strip().startswith('&'):
                    keyword = line.strip()[1:].split()[0]  # Remove '&' and get the keyword
                    output_file.write(line.replace(f"&{keyword}", f"{keyword}="))
                    
                
                else:
                    output_file.write(line)
    
                i += 1
            
            #Add the last / to close the file
            output_file.write('/\n')
        

       
        print(f"Processing complete. Output written to {output_file_path}")
                                                                           
    except FileNotFoundError:
        print(f"Error: The file '{input_file_name}' was not found.")
    except IOError:
        print(f"Error: An I/O error occurred while processing the files.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")



if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 script.py input_filename")
        sys.exit(1)
    
    input_file_name = sys.argv[1]
    output_file_name = os.path.splitext(input_file_name)[0] + '.namelists'
    process_file(input_file_name, output_file_name)



