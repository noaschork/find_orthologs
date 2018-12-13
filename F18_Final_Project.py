import subprocess
import sys
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO


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

#Creates a new files where the ClustalW alignment will be written to
out_file = "aligned_combined.aln"
#Takes in the combined file of human and non-human sequences and aligns them
muscle_cline=MuscleCommandline(input="combined.fasta", out=out_file, clw=True)
stdout, stderr = muscle_cline()
align=AlignIO.read(out_file, "clustal")


def cleanup_dict(clin_var_data_dict): #NK&NS
    '''
    Converts the amino acids from three to one letter codes.

    Takes in the clin_var_data_dict and accesses the positions containing the
    original and mutant amino acids. Depending on length, the function will
    either convert from three to one letter code or include an arrow (->) to
    indicate a frame shift mutation. In the conversion dictionary, terminating
    amino acids are converted into stars (*) as well as deletions. If the amino
    acid data is not in the dictionary, the data is passed over.

    Args:
    clin_var_data_dict - dictionary with accession numbers and data

    Returns:
    nothing is returned. Only modifies clin_var_data_dict.
    '''
    amino_conversion= {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C',
                       'Glu':'E','Gln':'Q','Gly':'G','His':'H','Ile':'I',
                       'Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P',
                       'Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V',
                       'Ter':'*', 'del':'*'}

    for accession in clin_var_data_dict:
        for entry in clin_var_data_dict[accession]:
            #Converts original amino acid to one letter code
            if len(entry[0]) == 0:
                pass
            elif len(entry[0]) == 3:
                if amino_conversion.get(entry[0]) == None:
                    pass
                else:
                    entry[0] = amino_conversion[entry[0]]
            else:
                pass
            #Converts mutant amino acid to one letter code
            if len(entry[2]) == 0:
                pass
            elif len(entry[2]) == 3:
                if amino_conversion.get(entry[2]) == None:
                    pass
                else:
                    entry[2] = amino_conversion[entry[2]]
            #If mutant results in a frameshift, arrow indicates it
            elif len(entry[2]) == 5:
                if amino_conversion.get(entry[2][0:3]) == None:
                    pass
                else:
                    entry[2] = amino_conversion[entry[2][0:3]] + '->'
            else:
                pass
def clin_var_parser(): #NK&NS
    '''
    Parses variant_summary.txt into accession number and data.

    Creates a dictionary which holds accession numbers as keys and each value is
    a list containing each entry under that accession number as a nested list.
    Each entry list contains: the original amino, the mutation position, the
    mutant amino, whether a deletion occurred ("yes"/"no"), mutation type, its
    clinical significance, the phenotype, and the chromosome the mutation occurs
    in. It excludes entries that lack protein information or if a large deletion
    and an insertion occurred.

    Returns:
    clin_var_data_dict - dictionary with accession numbers and data
    '''
    #Initializes the clin_var_data_dict
    clin_var_data_dict = {}

    with open('variant_summary2.txt', 'r') as fh:
        #Variable containing all column labels in fh
        column_labels = fh.readline()

        #List to contain each line in fh
        lines = []

        #Parses out the mRNA accession number of each line and creates a
        #dictionary entry for each accession number
        for line in fh:
            line = line.strip('\n')
            line = line.split('\t')
            lines.append(line)

            #Variable containing the accession number of the entry
            accession_num = line[2].split('(')[0]
            other_Name_data = line[2].split('(')[-1]

            #Creates an entry in dictionary for mRNA accession number
            if accession_num.find('NM') >= 0:
                if other_Name_data.find('p.') >= 0 and \
                   other_Name_data.find('delins') == -1 and \
                   other_Name_data.find('_') == -1:
                    clin_var_data_dict[accession_num] = []

        for line in lines:
            accession_num = line[2].split('(')[0]
            if accession_num.find('NM') >= 0:
                #Replacing first parentheses with star to find amino info
                line[2] = line[2].replace('(','*',1)

                current_line = line[2].split('(')

                amino_info = current_line[-1]
                original_amino = ''
                position = ''
                mutant_amino = ''
                deletion_y_n = ''

                #CASE: if no protein data is given
                if amino_info.find('p.') == -1:
                    #print 'noprotdata', amino_info
                    pass
                #CASE: if the mutation is a large deletion and insertion
                elif amino_info.find('delins') >= 0:
                    pass
                #CASE: if the nucleotide mutation is the same amino acid
                elif amino_info.find('=') >= 0:
                    original_amino = amino_info[2:5]
                    for char in amino_info:
                        if char.isdigit():
                            position += char
                    position = position
                    mutant_amino = original_amino
                    deletion_y_n = 'no'
                #CASE: if the mutation is a deletion
                elif amino_info.find('del') >= 0:
                    original_amino = amino_info[2:5]
                    for char in amino_info:
                        if char.isdigit():
                            position += char
                    position = position
                    mutant_amino = amino_info[(len(amino_info)-4):-1]
                    deletion_y_n = 'yes'
                #CASE: if the mutation is a frameshift
                elif amino_info.find('fs') >= 0:
                    original_amino = amino_info[2:5]
                    for char in amino_info:
                        if char.isdigit():
                            position += char
                    position = position
                    mutant_amino = amino_info[(len(amino_info)-6):-1]
                    deletion_y_n = 'no'
                #CASE: if the mutation leads to a termination (Ter)
                elif amino_info.find('Ter') >= 0:
                    original_amino = amino_info[2:5]
                    for char in amino_info:
                        if char.isdigit():
                            position += char
                    position = position
                    mutant_amino = amino_info[(len(amino_info)-4):-1]
                    deletion_y_n = 'no'
                #CASE: if the mutation is a simple substitution
                else:
                    original_amino = amino_info[2:5]
                    for char in amino_info:
                        if char.isdigit():
                            position += char
                    position = position
                    mutant_amino = amino_info[(len(amino_info)-4):-1]
                    deletion_y_n = 'no'

                #Appends the relevant information to the entry for each
                #accession number
                other_Name_data = line[2].split('(')[-1]
                if other_Name_data.find('p.') >= 0 and \
                   other_Name_data.find('delins') == -1 and \
                   other_Name_data.find('_') == -1:
                    clin_var_data_dict[accession_num].append([original_amino,\
                        position, mutant_amino, deletion_y_n, line[1],line[6],\
                        line[13], line[18]])


        #Removes duplicates
        for accession_num in clin_var_data_dict:
            copy = []
            for entry in clin_var_data_dict[accession_num]:
                if entry not in copy:
                    copy.append(entry)
            clin_var_data_dict[accession_num] = copy

        cleanup_dict(clin_var_data_dict)

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
    with open('aligned.aln', 'r') as aligned_seq:
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

    #Identify Conserved Clinical Variants & Format Alignment String
    index = 0
    cv_alignment_str = ''
    for symbol in symbols:
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
    start = 0
    #add in conserved indexes for strongly and weakly conserved nucleotides
    for symbol in cv_alignment_str:
        if symbol == "|":
            index = cv_alignment_str.find("|", start) #Find index position of conserved clinical variant
            start = index + 1
            index = str(index+1)
            conserved_cv_index.append(index) #Append conserved CV index position to list
        else:
            start = start + 1 #Iterate through cv_alignment_str

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
    hum_access = get_hum_access("human.fasta")
    clin_var_data_dict = clin_var_parser()
    symbols = convert_to_string(align)
    conserved_var_output = conserve_CVs(clin_var_data_dict, symbols, hum_access)
    print('output',conserved_var_output)
if __name__ == '__main__':
    main()
