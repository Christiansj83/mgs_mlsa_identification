# mgs_mlsa_identification
Repository containing a script taking in WGS (Whole genome sequencing) sequences from Mitis Group Streptococci and returns the concatenated MultiLocus Sequence Analysis (MLSA) sequence*. 

in your root dir make folders called 'blast_db', 'output' and 'blast_result'. Put your fasta sequences (You need to use the ".fasta"-extension, NOT ".fna", "faa" or anything else) in the root folder and run:

mlsa_script.py -q [fasta query sequence]



Dependencies: Blast+, Python3 (including pandas)

*see PMID 19171050

