#!/usr/bin/env python

import argparse
import pandas as pd
import glob
import os 

#Denne her er med de nye gener brugt til den nye mlsa

def ReverseComplement1(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])

def import_genome_names():
	files = glob.glob('*.fasta')
	return files

def gene_length_query(query):
    '''takes a fasta file and returns a dictionary with genes as keys and gene length as values'''
    fasta = open('%s' %(query), 'r')
    fasta_read = fasta.readlines()
    gene_length = {}
    for line in fasta_read:
        if line[0]=='>':
            length = ''
            gene = line[1:]
            gene = gene.rstrip()
            gene_length[gene] = None 
        else:
            value = line
            value = value.rstrip()
            length +=value
            gene_length[gene] = len(length)
    return gene_length

def number_of_genes(query):
    '''takes a fastefile and returns the number of genes in the fastafile'''
    fasta = open('%s' %(query), 'r')
    fasta_read = fasta.readlines()
    no_of_genes = 0
    for line in fasta_read:
        if line[0]=='>':
            no_of_genes = no_of_genes + 1
    return no_of_genes

def blasting(genome, query):
    '''takes a genome and make a blast search against the query'''
    #First making a blast database
    db_genome = '%s' %(genome)
    out_folder = 'blast_db/%s_db' %(genome)
    command_db = 'makeblastdb -in %s -out %s -dbtype nucl' %(db_genome, out_folder)
    os.system(command_db)
    #then i make the blast command and saves it in the folder blat_results
    db = 'blast_db/%s_db' %(genome)
    command = 'blastn -db %s -query %s -out blast_result/%s.txt -outfmt 7 -task blastn-short' %(db, query, genome)
    os.system(command)

def import_blast_result(genome, no_of_genes):
    #print(genome)
    header=['query_id', 'subject_id', 'perc_identity', 'alignment length', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end','evalue', 'bit_score', 'sequence']
    blast_df = pd.read_csv('blast_result/%s.txt' %(genome), sep='\t', names=header, comment='#')
    blast_df = blast_df.drop_duplicates(subset='query_id', keep="first")
    #Denne her har jeg lige ændret så de fjerner alle < 200 bit_Score
    #blast_df = blast_df.drop(blast_df.loc[blast_df['bit_score']<200].index)
    if len(blast_df)!=no_of_genes:
        print('%s has a number of gens different from expected' %(genome))
    return blast_df
    
def fasta_to_dict(genome):
    '''functions opens fasta/genome and returns a dictionery with headers as keys and sequences as value'''
    fasta = open('%s' %(genome), 'r')
    fasta_read = fasta.readlines()
    dict_fasta = {}
    
    for line in fasta_read:
        if line[0]=='>':
            line = line[1:]
            line = line.split(' ')
            key = line[0]
            key = key.rstrip()
            #key = line[1:]
            #key = key.split(' ')
            #print(key)
            #key= key[1]
            dict_fasta[key]= ''
        else:
            value = line
            value = value.rstrip()
            dict_fasta[key] +=value
            
    return dict_fasta

def extract_gene_sequences_with_concatenated(dict_fasta, blast_result, gene_length):
    '''Same as "extract_gene_sequences" (see above) but it also makes a concatenated fasta with all sequences from each
    strain concatenated'''
    rows = range(0,len(blast_result))
    dict_genes = {}
    dict_genes['concatenated'] = ''
    for n in rows:
        contig = blast_result.iloc[n]['subject_id']
        start = blast_result.iloc[n]['s_start']
        stop = blast_result.iloc[n]['s_end']
        start_query = blast_result.iloc[n]['q_start']
        stop_query = blast_result.iloc[n]['q_end']
        gene = blast_result.iloc[n]['query_id']
        #-1 fordi der mangler en base ved alle gener... uvist hvorfor.
        #Først ordner jeg de hits der er i "rigtig" rækkefølge, stop>start 
        start = int(start)
        #start_minus = start -1
        start_extended_right = start - start_query
        stop_extended_right = stop + (gene_length[gene]-stop_query)
        #stop_gene_length = start_minus + gene_length[gene]
        #Så ordner jeg dem der er i forkert rækkefølge, dvs. start>stop
        stop = blast_result.iloc[n]['s_end']
        stop = int(stop)
        stop_extended_wrong = stop + (blast_result.iloc[n]['q_end']-gene_length[gene]-1)
        start_extended_wrong = start + (blast_result.iloc[n]['q_start']-1)
        #stop_minus = stop - blast_result.iloc[n]['q_start']
        #reverse_start = stop_minus + gene_length[gene]
        #stop_minus = stop -1
        #Jeg går lige en base længere tilbage da der ved alignment mangler en base i forhold til reference
        #Jeg skal lige tjekke om jeg skal vende de gener om som den ikke finder
        if stop > start:
            sequence = dict_fasta[contig]
            #print(start_extended_right)
            #print(stop_extended_right)
            if start_extended_right < 0:
                start_extended_right = 0
            sequence = sequence[start_extended_right:stop_extended_right]
            dict_genes[gene] = sequence
            dict_genes['concatenated'] +=sequence
        elif stop < start:
            sequence = dict_fasta[contig]
            #print(stop_extended_wrong)
            #print(start_extended_wrong)
            if stop_extended_wrong < 0:
                stop_extended_wrong = 0
            sequence = sequence[stop_extended_wrong:start_extended_wrong]
            #sequence = sequence[stop_minus:start]
            sequence = ReverseComplement1(sequence)
            dict_genes[gene] = sequence
            dict_genes['concatenated'] +=sequence
    return dict_genes

def append_genes_to_fasta(gene_sequences, genome):
    for gene in list(gene_sequences.keys()):
        sequence = gene_sequences[gene]
        file = open('output/%s.txt' %(gene), 'a')
        file.write(">%s" '\n' '%s' '\n' %(genome, sequence))   
    return None

def main(args):
    list_of_genomes = import_genome_names()
    for genome in list_of_genomes:
        gene_length = gene_length_query(args.q)
        no_of_genes = number_of_genes(args.q)
        blasting(genome, args.q)
        blast_result = import_blast_result(genome, no_of_genes)
        dict_fasta_file = fasta_to_dict(genome)
        gene_sequences = extract_gene_sequences_with_concatenated(dict_fasta_file, blast_result, gene_length)
        append_genes_to_fasta(gene_sequences, genome)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts MLSA gene sequences')
    parser.add_argument('-q', help='query file containing genes in MLSA scheme', required=True)

    args = parser.parse_args()

    main(args)

    

