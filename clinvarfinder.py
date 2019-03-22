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
    Finds the human gene accession number

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
                hum_access = hum_access + line[0]
    return hum_access

def import_variant_summary():
    '''
    Imports the parsed clinical variant summary file

    Reads in the clin_var_data_dict from the file parse_clin_var.json. This file
    is the parsed version of the variant summary text.

    Returns:
    clin_var_data_dict - dictionary with accession numbers and data
    '''
    with open("parse_clin_var.json") as var_sum:
        clin_var_data_dict = json.load(var_sum)
    return clin_var_data_dict

def gene_of_interest(clin_var_data_dict,hum_access):
    '''
    Finds the gene that we are interested in from the larger clin_var_data_dict

    This function uses the inputted human accession number to find only the
    clinical variants that occur within the gene of var_interest

    Arg:
    clin_var_data_dict - dictionary with accession numbers and data
    hum_access - human accession number parsed from the sequences the user input

    Returns:
    var_interest - a list of all the clinical variants that occur within the
    gene of interest
    '''
    var_interest = []
    for gene in clin_var_data_dict.keys():
        if gene[:9] == hum_access:
            var_interest.append(clin_var_data_dict[gene])
    return var_interest

def read_in_file():
    '''
    Reads in the multiple sequence alignment (MSA) one line at a time, stripping
    off the new line characters and adds the, to a list

    Returns:
    file_list - list of the aligned sequences
    '''
    with open('aligned_combined.aln', 'r') as aligned_seq:
        file_list = []
        for line in aligned_seq:
            line = line.strip('\n')
            file_list.append(line)
    return file_list

def convert_to_string(file_list,hum_access):
    '''
    Converts the MSA to a string Type
    
    MUSCLE outputs the alignment in type Multipe Sequence Alignment (MSA), which is
    not mutable and so not very helpful for our needs. This function converts
    the alignment into a string.

    Args:
    file_list - list of the aligned sequences
    hum_access - human accession number parsed from the sequences the user input

    Returns:
    symbols - a tuple with a list of aligned symbols (*,:, or .) at position 0,
    list of human sequence at position 1, and a list of non-human sequences at
    position 2
    '''
    header = ''
    hum_sequence = ''
    symbols = ''
    non_hum_sequence = ''
    for pos,line in enumerate(file_list):
        if line.startswith('NM'):
            if line.startswith(hum_access):
                hum_sequence += line[20:]
            if not line.startswith(hum_access):
                non_hum_sequence += line[20:]
        else:
            if file_list[pos-1].startswith('NM'):
                symbols += line[20:]

    return symbols, hum_sequence, non_hum_sequence

def full_conserve_CVs(var_interest, symbols, hum_access):
    '''
    Identify Fully Conserved Clinical Variants Between Protein Sequences

    This method identifies the clinical variants fully conserved between the
    non-human protein sequence and the human protein sequence. Creates a list of
    the cv index positions. Determines if matches are fully conserved clinical
    variants and formats a new alignment string. Identify the fully conserved
    clinical variants and their parameters to generate data for a conserved
    clinical variant table.

    Args:
    var_interest - a list of all the clinical variants that occur within the
    gene of interest
    symbols - a tuple with alignment symbols, human sequence, and non-human
    sequences
    hum_access - the accession number of the human gene of interest

    Returns:
    full_var_output - nested list of fully conserved clinical variants and
    parameters
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

def strong_conserve_CVs(var_interest, symbols, hum_access):
    '''
    Identify Strongly Conserved Clinical Variants Between Protein Sequences

    This method identifies the clinical variants strongly conserved between the
    non-human protein sequence and the human protein sequence. Creates a list of
    the cv index positions. Determines if matches are strongly conserved clinical
    variants and formats a new alignment string. Identify the strongly conserved
    clinical variants and their parameters to generate data for a conserved
    clinical variant table.

    Args:
    var_interest - a list of all the clinical variants that occur within the
    gene of interest
    symbols - a tuple with alignment symbols, human sequence, and non-human
    sequences
    hum_access - the accession number of the human gene of interest

    Returns:
    strong_var_output - nested list of strongly conserved clinical variants and
    parameters
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

def weak_conserve_CVs(var_interest, symbols, hum_access):
    '''
    Identify Weakly Conserved Clinical Variants Between Protein Sequences

    This method identifies the clinical variants weakly conserved between the
    non-human protein sequence and the human protein sequence. Creates a list of
    the cv index positions. Determines if matches are weakly conserved clinical
    variants and formats a new alignment string. Identify the weakly conserved
    clinical variants and their parameters to generate data for a conserved
    clinical variant table.

    Args:
    var_interest - a list of all the clinical variants that occur within the
    gene of interest
    symbols - a tuple with alignment symbols, human sequence, and non-human
    sequences
    hum_access - the accession number of the human gene of interest

    Returns:
    weak_var_output - nested list of weakly conserved clinical variants and
    parameters
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

def not_conserve_CVs(var_interest, symbols, hum_access):
    '''
    Identify Non-Conserved Clinical Variants Between Protein Sequences

    This method identifies the clinical variants non-conserved between the
    non-human protein sequence and the human protein sequence. Creates a list of
    the cv index positions. Determines if matches are non-conserved clinical
    variants and formats a new alignment string. Identify the non -conserved
    clinical variants and their parameters to generate data for a conserved
    clinical variant table.

    Args:
    var_interest - a list of all the clinical variants that occur within the
    gene of interest
    symbols - a tuple with alignment symbols, human sequence, and non-human
    sequences
    hum_access - the accession number of the human gene of interest

    Returns:
    not_var_output - nested list of non-conserved clinical variants and
    parameters
    '''
    index = 0
    not_index = []
    symbol_list = list(symbols[0])
    for character in symbol_list:
        if character == " ":
            position = symbol_list.index(character) + 1
            not_index.append(position)
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

    not_hum_index = []
    for pos in not_index:
        for loc in hum_positions:
            if pos == hum_positions.index(loc):
                not_hum_index.append(loc)

    not_var_output = []
    for index in not_hum_index:
        for line in var_interest:
            for cv in line:
                if index == int(cv[1]):
                    not_var_output.append(cv)
    return not_var_output

def export_file(full_var_output, strong_var_output, weak_var_output,
not_var_output):
    '''
    Creates a Tab Delimited File with the Clinical Variants data

    This function exports a tab delimited file with clinical variants found in
    fully, strongly, weakly, and non conserved regions between the inputted
    human and non-human sequences. The file provides additional information on
    the chromosome numner, position, original amino acid, mutant amino acid,
    whether a deletion occurs, mutaiton type, clinical significance and the
    phenotype of each variant.

    Args:
    full_var_output - nested list of fully conserved clinical variants and
    parameters
    strong_var_output - nested list of strongly conserved clinical variants and
    parameters
    weak_var_output - nested list of weakly conserved clinical variants and
    parameters
    not_var_output - nested list of non-conserved clinical variants and
    parameters
    '''
    with open("clinvar.txt", "w") as f:
        f.write("Chr num" + "\t" + "Pos" + "\t" + "Orig AA" + "\t" + "Mut AA" +
        "\t" + "Del Y/N" + "\t" + "Mutation Type" + "\t" + "Clin Sig" + "\t" +
        "Phenotype" + "\t" + "Conservation"+ "\n")
        for entry in full_var_output:
            f.write(entry[7] + "\t" + entry[1] + "\t" + entry[0] + "\t" +
            entry[2] + "\t" + entry[3] + "\t" + entry[4] + "\t" + entry[5] +
            "\t" + entry[6]+ "\t" + "Full" + "\n")
        for entry in strong_var_output:
            f.write(entry[7] + "\t" + entry[1] + "\t" + entry[0] + "\t" +
            entry[2] + "\t" + entry[3] + "\t" + entry[4] + "\t" + entry[5] +
            "\t" + entry[6] + "\t" + "Strong" + "\n")
        for entry in weak_var_output:
            f.write(entry[7] + "\t" + entry[1] + "\t" + entry[0] + "\t" +
            entry[2] + "\t" + entry[3] + "\t" + entry[4] + "\t" + entry[5] +
            "\t" + entry[6] + "\t" + "Weak" + "\n")
        for entry in not_var_output:
            f.write(entry[7] + "\t" + entry[1] + "\t" + entry[0] + "\t" +
            entry[2] + "\t" + entry[3] + "\t" + entry[4] + "\t" + entry[5] +
            "\t" + entry[6] + "\t" + "No" + "\n")

def list_var_pos(full_var_output, strong_var_output, weak_var_output, not_var_output):
    '''
    This function is not finished yet

    Aims to add a line to the aligned sequence file denoting where the Clinical
    variants occur within the alignment

    Args:
    full_var_output - nested list of fully conserved clinical variants and
    parameters
    strong_var_output - nested list of strongly conserved clinical variants and
    parameters
    weak_var_output - nested list of weakly conserved clinical variants and
    parameters
    not_var_output - nested list of non-conserved clinical variants and
    parameters

    Returns:
    var_list - list of all the clinical variants
    '''
    var_list = []
    for entry in full_var_output:
        var_list.append(entry[1])
    for entry in strong_var_output:
        var_list.append(entry[1])
    for entry in weak_var_output:
        var_list.append(entry[1])
    for entry in not_var_output:
        var_list.append(entry[1])
    return var_list

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
    not_var_output = not_conserve_CVs(var_interest, symbols, hum_access)
    export_file(full_var_output,strong_var_output, weak_var_output, not_var_output)
    var_list = list_var_pos(full_var_output,strong_var_output, weak_var_output, not_var_output)

if __name__ == '__main__':
    main()
