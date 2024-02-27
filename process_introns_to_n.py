import sys

def process_fasta_with_utrs(input_file_path, output_file_path, utr_boolean):
    """
    Processes a FASTA file to replace introns (lowercase characters) with 'N's, 
    while retaining the UTR sequences (initial and final sets of lowercase characters) 
    and writes the output to a new FASTA file.

    :param input_file_path: Path to the input FASTA file.
    :param output_file_path: Path to the output FASTA file.
    """
    with open(input_file_path, 'r') as input_file:
        lines = input_file.readlines()
    
    processed_lines = []

    for line in lines:
        if line.startswith('>'):
            processed_lines.append(line)
        else:
            start_utr = ''
            end_utr = ''
            middle_sequence = line.strip()
            
            for char in middle_sequence:
                if char.islower():
                    start_utr += char
                else:
                    break
            middle_sequence = middle_sequence[len(start_utr):]

            for char in reversed(middle_sequence):
                if char.islower():
                    end_utr = char + end_utr
                else:
                    break
            middle_sequence = middle_sequence[:-len(end_utr)] if end_utr else middle_sequence

            processed_middle_sequence = ''.join(['N' if char.islower() else char for char in middle_sequence])
            if utr_boolean:
            	processed_sequence = start_utr + processed_middle_sequence + end_utr + '\n'
            else:
            	processed_sequence = processed_middle_sequence + '\n'
            processed_lines.append(processed_sequence)
    
    with open(output_file_path, 'w') as output_file:
        output_file.writelines(processed_lines)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script_name.py <input_fasta_file> <output_fasta_file> <retain UTRs True/False>")
        sys.exit(1)
    print(sys.argv)
    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]
    utr_boolean = sys.argv[3]
	
    process_fasta_with_utrs(input_file_path, output_file_path, utr_boolean)
