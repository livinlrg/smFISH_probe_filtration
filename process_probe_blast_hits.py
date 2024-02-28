import sys
import operator

def process_probe_blast_hits(input_file_path, output_file_path, gene_name):
    """
    Processes a blast hit file to fine the number of off-target hits

    :param input_file_path: Path to the input FASTA file.
    :param output_file_path: Path to the output FASTA file.
    :param gene_name
    
    Will separate the alignment file into a description section (gene hit stats) and a
    probe_section (alignments), and build a dictionary to be able to interact with these
    files.
    
    """
    '''
    RID: XA6CJ91K016	
	Job Title:48 sequences (Probe_1)				
	Program: BLASTN 
	Database: refseq_rna NCBI Transcript Reference Sequences
	Query #1: Probe_1 Query ID: lcl|Query_1571020 Length: 20
	
	Sequences producing significant alignments:
	                                                                  Scientific      Common                     Max    Total Query   E   Per.   Acc.                        
	Description                                                       Name            Name            Taxid      Score  Score cover Value Ident  Len        Accession        
	
	then
	Alignments:
	
	Query #2: Probe_2 Query ID: lcl|Query_1571021 Length: 20
	
	Sequences producing significant alignments:
	                                                                  Scientific      Common                     Max    Total Query   E   Per.   Acc.                        
	Description                                                       Name            Name            Taxid      Score  Score cover Value Ident  Len        Accession        
	Caenorhabditis briggsae Protein CBR-CEP-1 (Cbr-cep-1), partial... Caenorhabdit... NA              6238       40.1   40.1  100%  1e-04 100.00 1977       XM_002639437.2     
    '''
    with open(input_file_path, 'r') as input_file:
        lines = input_file.readlines()
    
    out_line = ''
    
    desc_header = ["Description", "Name", "Name", "Taxid", "Score", "Score", "cover", "Value", "Ident", "Len", "Accession"]
    desc_line_space = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    probe_off_target = []
    probe_number = []
    probe_length = []
    probe_hits_line = []
    hit_dict_list = {}
    gene_hit_dict = {}
    i = -1
    probe_section = False
    previous_query_over = True
    
    for line in lines:
    	print(str(len(line)) + " " + line.strip())
    	if line.startswith('Query #'):
    		# Determines the length and name of the probe
    		# Should always be at the start of a probe section
    		i += 1
    		split_line = line.split(':')
    		#print(line)
    		probe_number.append(split_line[1].split(' ')[1])
    		probe_length.append(int(split_line[3][1:]))
    		#print(probe_number[i] + ' ' + str(probe_length[i]))
    		probe_hits_line.append([])
    		hit_dict_list[probe_number[i]] = []
    	elif line.startswith('Description'):
    		# Calculate line spacing
    		probe_section = True
    		
    		for j in range(1, len(desc_header)):
    			desc_line_space[j] = line.index(desc_header[j], desc_line_space[j - 1] + 1)
    			
    		continue
    	elif line.startswith('Alignments'):
    		# denote the start of the alignment section
    		probe_section = False
    		hit_dict = {}
    		hit_number = -1
    	    
    	if probe_section and len(line) > 1:
    		# denotes the probe_section
    		print(probe_section)
    		line_split = []
    		for j in range(10):
    			line_split.append(line[desc_line_space[j]:desc_line_space[j + 1]].strip())
    		print(line_split)
    		print(i)
    		probe_hits_line[i].append(line_split)
    	if not probe_section and line.startswith('>') and previous_query_over:
    		# process the full gene name of the alignment
    		previous_query_over = False
    		if ", partial mRNA" in line:
    			hit_name = line[24:line.index(', partial mRNA')]
    		elif ", mRNA" in line:
    			hit_name = line[24:line.index(', mRNA')]
    		hit_number += 1
    		
    		probe_hits_line[i][hit_number][0] = hit_name # update hit hit name to match
    		hit_dict_list[probe_number[i]].append(0)
    	if not probe_section and line.startswith('>') and not previous_query_over:
    		# Deal with multiple hits to paralogs
    		continue
    	if not probe_section and line.startswith('Range'):
    		# Extracts the length of the alignment
    		previous_query_over = True
    		hit_length = int(line[line.index(':'):].split('to')[1].strip()) - int(line[line.index(':') + 2:].split('to')[0].strip()) + 1
    		hit_dict_list[probe_number[i]][hit_number] = max(hit_dict_list[probe_number[i]][hit_number], hit_length)
    
    for i in range(len(probe_hits_line)):
    	for hit in probe_hits_line[i]:
    		if hit[0] not in gene_hit_dict.keys():
    			gene_hit_dict[hit[0]] = [0] * len(probe_hits_line)
	
    probe_out_list = []
    for i in range(len(probe_hits_line)):
    	
    	out_list = [probe_number[i], probe_length[i], 0, 0, 0, 0, 0]
    	j = 0
    	for hit in probe_hits_line[i]:
    		max_hit_length = hit_dict_list[probe_number[i]][j]
    		if gene_name in hit[0]:
    			out_list[2] += 1
    			out_list[3] = max(out_list[3], round(float(hit[6][:-1])/100 * int(probe_length[i]) * (float(hit[8])/100)))
    		elif gene_name not in hit[0]:
    			out_list[4] += 1
    			out_list[5] = max(out_list[5], round(float(hit[6][:-1])/100 * int(probe_length[i]) * (float(hit[8])/100)))
    			out_list[6] = max(out_list[6], round(max_hit_length * float(hit[8][:-1])/100))

    			gene_hit_dict[hit[0]][i] = max(gene_hit_dict[hit[0]][i], out_list[6])
    		j += 1
    	probe_out_list.append(out_list)
    #print(dict(sorted(gene_hit_dict.items(), key=operator.itemgetter(1))))
    
    for key in gene_hit_dict.keys():
    	gene_hit_dict[key]
    	if sum(x > 15 for x in gene_hit_dict[key]) > 2:
    		print(key + ": " + str([idx + 1 for idx, val in enumerate(gene_hit_dict[key]) if val > 15]))
    
    with open(output_file_path, 'a') as output_file:
    	output_file.write("probe_number" + "\t" + "probe_length" + "\t" + "num_on_target_hits" + "\t" + "max_on_target_length" + "\t" + "num_off_target_hits" + "\t" + "max_off_hit_length_query" + "\t" + "max_off_hit_length_actual" + "\t" +"\n")
    	for i in range(len(probe_out_list)):
    		out_text = ""
    		for j in range(len(probe_out_list[i])):
    			out_text = out_text + str(probe_out_list[i][j]) + "\t"
    		out_text += "\n"
    		output_file.write(out_text)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script_name.py <input_fasta_file> <output_fasta_file>")
        sys.exit(1)

input_file_path = sys.argv[1]
output_file_path = sys.argv[2]
gene_name = sys.argv[3]

process_probe_blast_hits(input_file_path, output_file_path, gene_name)
