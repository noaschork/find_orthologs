import subprocess
import sys
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO
import json

def combine_fasta(filename1, filename2): #DM
    '''
    This function combines two FASTA files into once

    It takes in two different FASTA files, one with the human sequence and the
    other with the nonhuman sequences and combines the two so that they can be
    aligned using MUSCLE

    Args:
    filename1 - the non-human FASTA file
    filename2 - the human FASTA file

    Returns:
    combined.fasta - a new file with the combined sequences
    '''
    with open("combined.fasta", "w") as output:
        with open(filename1, "r") as f1:
            line = f1.readline()
            while line:
                output.write(line)
                line = f1.readline()
        with open(filename2, "r") as f2:
            output.write('\n')
            line = f2.readline()
            while line:
                output.write(line)
                line = f2.readline()

def get_hum_access(filename2):
    '''
    This function parses out the accession number of the human gene from
    filename2. This will later be used to search the ClinVar file for only
    relevant mutations.

    Args:
    filename2 - the human FASTA file

    Returns:
    hum_access - the accession of the human gene
    '''
    hum_access = ''
    with open(filename2, 'r') as f2:
        for line in f2:
            if line[0] == '>':
                line = line.strip('>')
                line = line.split(' ')
                print(line[0])
                hum_access = hum_access + line[0]
    return hum_access

def import_variant_summary():
    '''
    Reads in the clin_var_data_dict from the file parse_clin_var.json. This file
    is the parsed version of the variant summary text.

    Returns:
    clin_var_data_dict - dictionary with accession numbers and data
    '''
    with open("parse_clin_var.json") as var_sum:
        clin_var_data_dict = json.load(var_sum)
    return clin_var_data_dict

def gene_of_interest(clin_var_data_dict,hum_access):
    var_interest = []
    for gene in clin_var_data_dict.keys():
        if gene[:9] == hum_access:
            var_interest.append(clin_var_data_dict[gene])
    return var_interest

def read_in_file():
    with open('aligned_combined.aln', 'r') as aligned_seq:
        file_list = []
        for line in aligned_seq:
            line = line.strip('\n')
            file_list.append(line)
    return file_list

def convert_to_string(file_list,hum_access):
    '''
    MUSCLE outputs the alignment in type Multipe Sequence Alignment (MSA), which is
    not mutable and so not very helpful for our needs. This function converts
    the alignment into a string.

    Args:
    align - the MSA output from MUSCLE

    Returns:
    symbols - a list of the alignment symbols (*,:, or .)
    '''
    header = ''
    hum_sequence = ''
    symbols = ''
    for pos,line in enumerate(file_list):
        if line.startswith('NM'):
            if line.startswith(hum_access):
                hum_sequence += line[20:]
        else:
            if file_list[pos-1].startswith('NM'):
                symbols += line[20:]

    return symbols, hum_sequence


def full_conserve_CVs(var_interest, symbols, hum_access): #KG&AP&AC
    '''
    Identify Conserved Clinical Variants Between Protein Sequences

    This method identifies the clinical variants conserved between the non-human
    protein sequence and the human protein sequence. Creates a list of the cv
    index positions. Determines if matches are conserved clinical variants and
    formats a new alignment string. Identify the conserved clinical variants and
    their parameters to generate data for a conserved clinical variant table.

    Args:
    clin_var_data_dict - dictionary of accession numbers and associated clinical
    variants
    symbols - list of alignment symbols
    hum_access - the accession number of the human gene of interest

    Returns:
    conserved_var_output - nested list of CONSERVED clinical variants and parameters


    '''
    index = 0
    full_index = []
    symbol_list = list(symbols[0])
    for character in symbol_list:
        if character == "*":
            position = symbol_list.index(character) + 1
            full_index.append(position)
            symbol_list[symbol_list.index(character)] = '/'

    hum_positions = []
    hum_list = symbols[1]
    counter = 0
    for char in hum_list:
        if char == "-":
            hum_positions.append(counter)
            counter += 0
        if char!= "-":
            hum_positions.append(counter)
            counter += 1

    full_hum_index = []
    for pos in full_index:
        for loc in hum_positions:
            if pos == hum_positions.index(loc):
                full_hum_index.append(loc)

    full_var_output = []
    for index in full_hum_index:
        for line in var_interest:
            for cv in line:
                if index == int(cv[1]):
                    full_var_output.append(cv)
    return full_var_output

def strong_conserve_CVs(var_interest, symbols, hum_access): #KG&AP&AC
    '''
    Identify Conserved Clinical Variants Between Protein Sequences

    This method identifies the clinical variants conserved between the non-human
    protein sequence and the human protein sequence. Creates a list of the cv
    index positions. Determines if matches are conserved clinical variants and
    formats a new alignment string. Identify the conserved clinical variants and
    their parameters to generate data for a conserved clinical variant table.

    Args:
    clin_var_data_dict - dictionary of accession numbers and associated clinical
    variants
    symbols - list of alignment symbols
    hum_access - the accession number of the human gene of interest

    Returns:
    conserved_var_output - nested list of CONSERVED clinical variants and parameters


    '''
    index = 0
    strong_index = []
    symbol_list = list(symbols[0])
    for character in symbol_list:
        if character == ":":
            position = symbol_list.index(character) + 1
            strong_index.append(position)
            symbol_list[symbol_list.index(character)] = '/'

    hum_positions = []

    hum_list = symbols[1]
    counter = 0
    for char in hum_list:
        if char == "-":
            hum_positions.append(counter)
            counter += 0
        if char!= "-":
            hum_positions.append(counter)
            counter += 1

    strong_hum_index = []
    for pos in strong_index:
        for loc in hum_positions:
            if pos == hum_positions.index(loc):
                strong_hum_index.append(loc)

    strong_var_output = []
    for index in strong_hum_index:
        for line in var_interest:
            for cv in line:
                if index == int(cv[1]):
                    strong_var_output.append(cv)
    return strong_var_output

def weak_conserve_CVs(var_interest, symbols, hum_access): #KG&AP&AC
    '''
    Identify Conserved Clinical Variants Between Protein Sequences

    This method identifies the clinical variants conserved between the non-human
    protein sequence and the human protein sequence. Creates a list of the cv
    index positions. Determines if matches are conserved clinical variants and
    formats a new alignment string. Identify the conserved clinical variants and
    their parameters to generate data for a conserved clinical variant table.

    Args:
    clin_var_data_dict - dictionary of accession numbers and associated clinical
    variants
    symbols - list of alignment symbols
    hum_access - the accession number of the human gene of interest

    Returns:
    conserved_var_output - nested list of CONSERVED clinical variants and parameters


    '''
    index = 0
    weak_index = []
    symbol_list = list(symbols[0])
    for character in symbol_list:
        if character == ".":
            position = symbol_list.index(character) + 1
            weak_index.append(position)
            symbol_list[symbol_list.index(character)] = '/'

    hum_positions = []
    hum_list = symbols[1]
    counter = 0
    for char in hum_list:
        if char == "-":
            hum_positions.append(counter)
            counter += 0
        if char!= "-":
            hum_positions.append(counter)
            counter += 1

    weak_hum_index = []
    for pos in weak_index:
        for loc in hum_positions:
            if pos == hum_positions.index(loc):
                weak_hum_index.append(loc)

    weak_var_output = []
    for index in weak_hum_index:
        for line in var_interest:
            for cv in line:
                if index == int(cv[1]):
                    weak_var_output.append(cv)
    return weak_var_output

def export_file(full_var_output, strong_var_output, weak_var_output):
    with open("clinvar.txt", "w") as f:
        f.write("Clinical Variants in Fully Conserved Regions\n")
        f.write("Chr num\t")
        f.write("Pos\t")
        f.write("Original AA\t")
        f.write("Mutant AA\t")
        f.write("Deletion Y/N\t")
        f.write("Mutation Type\t")
        f.write("Clinical significance\t")
        f.write("Phenotype\n")
        for entry in full_var_output:
            f.write(entry[7])
            f.write("\t")
            f.write(entry[1])
            f.write("\t")
            f.write(entry[0])
            f.write("\t")
            f.write(entry[2])
            f.write("\t")
            f.write(entry[3])
            f.write("\t")
            f.write(entry[4])
            f.write("\t")
            f.write(entry[5])
            f.write("\t")
            f.write(entry[6])
            f.write("\n")
        f.write("\nClinical Variants in Strongly Conserved Regions\n")
        f.write("Chr num\t")
        f.write("Pos\t")
        f.write("Original AA\t")
        f.write("Mutant AA\t")
        f.write("Deletion Y/N\t")
        f.write("Mutation Type\t")
        f.write("Clinical significance\t")
        f.write("Phenotype\n")
        for entry in strong_var_output:
            f.write(entry[7])
            f.write("\t")
            f.write(entry[1])
            f.write("\t")
            f.write(entry[0])
            f.write("\t")
            f.write(entry[2])
            f.write("\t")
            f.write(entry[3])
            f.write("\t")
            f.write(entry[4])
            f.write("\t")
            f.write(entry[5])
            f.write("\t")
            f.write(entry[6])
            f.write("\n")
        f.write("\nClinical Variants in Weakly Conserved Regions\n")
        f.write("Chr num\t")
        f.write("Pos\t")
        f.write("Original AA\t")
        f.write("Mutant AA\t")
        f.write("Deletion Y/N\t")
        f.write("Mutation Type\t")
        f.write("Clinical significance\t")
        f.write("Phenotype\n")
        for entry in weak_var_output:
            f.write(entry[7])
            f.write("\t")
            f.write(entry[1])
            f.write("\t")
            f.write(entry[0])
            f.write("\t")
            f.write(entry[2])
            f.write("\t")
            f.write(entry[3])
            f.write("\t")
            f.write(entry[4])
            f.write("\t")
            f.write(entry[5])
            f.write("\t")
            f.write(entry[6])
            f.write("\n")

def main():
    combine_fasta("NR5A_orthologs.fasta", "human.fasta")
    #Creates a new files where the ClustalW alignment will be written to
    out_file = "aligned_combined.aln"
    #Takes in the combined file of human and non-human sequences and aligns them
    muscle_cline=MuscleCommandline(input="combined.fasta", out=out_file, clw=True)
    stdout, stderr = muscle_cline()
    align=AlignIO.read(out_file, "clustal")
    hum_access = get_hum_access("human.fasta")
    clin_var_data_dict = import_variant_summary()
    file_list = read_in_file()
    var_interest = gene_of_interest(clin_var_data_dict,hum_access)
    symbols = convert_to_string(file_list,hum_access)
    full_var_output = full_conserve_CVs(var_interest, symbols, hum_access)
    strong_var_output = strong_conserve_CVs(var_interest, symbols, hum_access)
    weak_var_output = weak_conserve_CVs(var_interest, symbols, hum_access)
    export_file(full_var_output,strong_var_output, weak_var_output)
if __name__ == '__main__':
    main()
