# smFISH probe design using Stellaris

### For creating smFISH probes that are filtered for off-target hits to a species genome.

1. Download unspliced mRNA transcript sequence from wormbase.
Replace the introns with N’s and possibly remove the UTR’s using the python script: process_introns_to_n.py
```
python process_introns_to_n.py <input_file_path (.fasta)> <output_file_path (.fasta)> <retain UTR (True/False)>
```

Example:
```
python process_introns_to_n.py sequences/Cbr-gld-1.fasta sequences/Cbr-gld-1_n_introns.fasta False
```

3. Input N replaced sequence into the Stellaris probe designer (https://www.biosearchtech.com/stellaris-designer).
Copy the probe set into a google sheet file for later
Using either google sheets or a text editor like BBEdit, construct a fasta file from the probe number and the probe sequences. An example of this for BBEdit with the grep box checked would be:

Find and replace ^ with >Probe_
Find and replace \t with \n

The final output would look something like:
```
>Probe_1
cttccgtgttcattttgatg
>Probe_2
acgctaagttgagaatccgt
>Probe_3
tcgagttgtatttccatcga
```

3. Copy the newly created fasta sequences and paste into NCBI blast (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome). Select the organism and change the database to refseq_rna.

4. Download the BLAST hits as a Text file.

5. Score the probe set by the number of off-target hits using the python script process_probe_blast_hits.py. The program tries to parse the BLAST output and score the probes by how many on-target hits it has (to the gene_name you provide), how many off-target hits it has, and the length of the off-target hits weighted by the percent match to the probe. In the terminal window, it’ll also output the names of genes that were hit by at least 16 bp of a probe, and hit at least 3 probes. The number of these probes will be output (see below). I’m sure there will be weird edge cases that will break the script, so just let me know if it fails.
```
python process_probe_blast_hits.py <input_blast_hits> <output_probe_scores> gene_name
```
Example:
```
python process_probe_blast_hits.py sequences/Cel-cep-1_Alignment.txt sequences/Cel-cep-1_probe_scores.txt cep-1
```
Example of terminal output:
```
 Protein CBG18920 (CBG18920): [1, 3, 30, 32]
 Protein CBR-EGL-19 (Cbr-egl-19): [1, 6, 16]
 Protein CBR-BET-2 (Cbr-bet-2): [1, 23, 34, 48]
```
The output file has a number of columns:
probe_number		probe_length	num_on_target_hits	max_on_target_length	num_off_target_hits	max_off_hit_length_query	max_off_hit_length_actual	

Most of these are self explanatory. The max_off_hit_length_query is the probe_length times the query cover and the percent ID. It can be off because BLAST gives multiple hits to the same gene the maximum query cover. Instead, the max_off_hit_length_actual takes the maximum single off-target hit, then multiples that by the percent ID. For filtering probes, I imagine looking at the on-target hit rate (sanity check), the terminal output, and the max_off_hit_length_actual would be the best filters.
