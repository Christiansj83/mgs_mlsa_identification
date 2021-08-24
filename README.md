# mgs_mlsa_identification
Repository containing a script taking in WGS (Whole genome sequencing) sequences from Mitis Group Streptococci and returns the concatenated MultiLocus Sequence Analysis (MLSA) sequence*. 

in your root dir make folders called 'blast_db', 'output' and 'blast_result' together with your "fasta-query-sequence"-file (NOT using the ".fasta"-extension, use instead eg. the ".fa"-extension. Put your fasta sequences (You need to use the ".fasta"-extension, NOT ".fna", "faa" or anything else) in the root folder and run:

mlsa_script.py -q [fasta query sequence]

In the "output-folder it returns the individual genes for each strain in a fasta-format, together with the concatenated genes appearing in the order found at the "fasta query sequence".

Dependencies: Blast+, Python3 (including pandas)

*see PMID 19171050

