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

def convert_to_string(align):
    '''
    MUSCLE outputs the alignment in type Multipe Sequence Alignment (MSA), which is
    not mutable and so not very helpful for our needs. This function converts
    the alignment into a string.

    Args:
    align - the MSA output from MUSCLE

    Returns:
    symbols - a list of the alignment symbols (*,:, or .)
    '''
    with open('aligned_combined.aln', 'r') as aligned_seq:
        header = ''
        symbols = ''
        for line in aligned_seq:
            #takes the header of the file and sets it to its own variable name
            if line[0] == 'M':
                header = header + line
            #if the line does not contain the header or an accession number,
            #it appends the contents of the line to the symbols string
            if line[0] == ' ':
                line = line.strip('\n')
                symbols = symbols + line
    return symbols

def conserve_CVs(clin_var_data_dict, symbols, hum_access): #KG&AP&AC
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


    #Create List of All Human CV Index Positions
    cv_aa_positions = []
    for nest in clin_var_data_dict.values():
        for inner in nest:
            cv_aa_positions.append(inner[1])
    print(cv_aa_positions)
    #Identify Conserved Clinical Variants & Format Alignment String
    index = 0
    cv_alignment_str = ''
    for symbol in cv_alignment_str:
        if symbol == "*":
            index_pos = symbols.find("*", index, len(symbols))
            # Find the first aa match from index position
            aa_pos = index_pos + 1 #Convert index position to aa position
            aa_pos_str = str(aa_pos) #Convert aa to string to search var_output
            if aa_pos_str in cv_aa_positions: #If aa position in cv index list
                cv_alignment_str = cv_alignment_str + "|"
            else: #Else, conserve the alignment symbol from align_ortho
                cv_alignment_str = cv_alignment_str + symbol
            index = index + 1
        else:  #Else, conserve the alignment symbol from align_ortho
            cv_alignment_str = cv_alignment_str + symbol
            index = index + 1

    #Identifying Conserved Clinical Variant and Parameters (Table)
    conserved_cv_index = [] #List of Conserved CV Index positions
    # start = 0
    # #add in conserved indexes for strongly and weakly conserved nucleotides
    # for symbol in cv_alignment_str:
    #     if symbol == "|":
    #         index = cv_alignment_str.find("|", start) #Find index position of conserved clinical variant
    #         start = index + 1
    #         index = str(index+1)
    #         conserved_cv_index.append(index) #Append conserved CV index position to list
    #     else:
    #         start = start + 1 #Iterate through cv_alignment_str

    conserved_var_output = [] # # AP, KG, AC: List of Conserved CVs w/ Parameters (for Table)
    for index in conserved_cv_index: # For every conserved index...
        for key in clin_var_data_dict:
            key_short = key[:len(hum_access)]
            if key_short == hum_access:
                print('key_short')
        #for cv in clin_var_data_dict[hum_access]: # and for every cv in var_output from id_cv...\
                for item in clin_var_data_dict[key]:
                    #print(item)
                    if index == item[1]: # if conserved cv index is present in the table
                        print("index")
                        conserved_var_output.append(item[1]) # append conserved cvs to table list
                    else:
                        pass
    return conserved_var_output

    #return [align, conserved_var_output] # AP, KG, AC: list containing String of cosnerved Alignment symbols and conserved clincial variants

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
    symbols = convert_to_string(align)
    conserved_var_output = conserve_CVs(clin_var_data_dict, symbols, hum_access)
    print('output',conserved_var_output)
if __name__ == '__main__':
    main()
