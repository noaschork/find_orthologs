import subprocess
import sys
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO

#takes in the test FASTA file and aligns all of the sequences within the file.
#This includes both the human and nonhuman proteins.
out_file = "aligned.aln"
#uses * to indicate a single, fully conserved nucleotide, : to indicate
#conservation between groups of strongly similar properties, and . to indicates
#conservation between groups of weakly similar properties
muscle_cline=MuscleCommandline(input="NR5A_orthologs.fasta", out=out_file, clw=True)
stdout, stderr = muscle_cline()
#align=AlignIO.read(StringIO(stdout), "fasta")
align=AlignIO.read(out_file, "clustal")



def cleanup_dict(clin_var_data_dict):
    '''
    ANS-NK: Converts the amino acids from three to one letter codes.

    Takes in the clin_var_data_dict and accesses the positions containing the
    original and mutant amino acids. Depending on length, the function will
    either convert from three to one letter code or include an arrow (->) to
    indicate a frame shift mutation. In the conversion dictionary, terminating
    amino acids are converted into stars (*) as well as deletions. If the amino
    acid data is not in the dictionary, the data is passed over.

    Args: clin_var_data_dict - dictionary with accession numbers and data
    ...
    Returns: nothing is returned. Only modifies clin_var_data_dict.
    '''
    amino_conversion= {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C',
                       'Glu':'E','Gln':'Q','Gly':'G','His':'H','Ile':'I',
                       'Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P',
                       'Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V',
                       'Ter':'*', 'del':'*'}

    for accession in clin_var_data_dict:
        for entry in clin_var_data_dict[accession]:
            # ANS-NK:Converts original amino acid to one letter code
            if len(entry[0]) == 0:
                pass
            elif len(entry[0]) == 3:
                if amino_conversion.get(entry[0]) == None:
                    pass
                else:
                    entry[0] = amino_conversion[entry[0]]
            else:
                pass
            # ANS-NK:Converts mutant amino acid to one letter code
            if len(entry[2]) == 0:
                pass
            elif len(entry[2]) == 3:
                if amino_conversion.get(entry[2]) == None:
                    pass
                else:
                    entry[2] = amino_conversion[entry[2]]
            #ANS-NK:CASE: if mutant results in a frameshift, arrow indicates it
            elif len(entry[2]) == 5:
                if amino_conversion.get(entry[2][0:3]) == None:
                    pass
                else:
                    entry[2] = amino_conversion[entry[2][0:3]] + '->'
            else:
                pass
def clin_var_parser():
    '''
    ANS-NK: Parses variant_summary.txt into accession number and data.

    Creates a dictionary which holds accession numbers as keys and each value is
    a list containing each entry under that accession number as a nested list.
    Each entry list contains: the original amino, the mutation position, the
    mutant amino, whether a deletion occurred ("yes"/"no"), mutation type, its
    clinical significance, the phenotype, and the chromosome the mutation occurs
    in. It excludes entries that lack protein information or if a large deletion
    and an insertion occurred.

    Returns: clin_var_data_dict - dictionary with accession numbers and data
    '''
    # ANS-NK:Initializes the clin_var_data_dict
    clin_var_data_dict = {}

    with open('variant_summary2.txt', 'r') as fh:
        # ANS-NK:Variable containing all column labels in fh
        column_labels = fh.readline()

        # ANS-NK:List to contain each line in fh
        lines = []

        # ANS-NK:Parses out the mRNA accession number of each line and creates a
        # dictionary entry for each accession number
        for line in fh:
            line = line.strip('\n')
            line = line.split('\t')
            lines.append(line)

            # ANS-NK:Variable containing the accession number of the entry
            accession_num = line[2].split('(')[0]
            other_Name_data = line[2].split('(')[-1]

            # ANS-NK:Creates an entry in dictionary for mRNA accession number
            if accession_num.find('NM') >= 0:
                if other_Name_data.find('p.') >= 0 and \
                   other_Name_data.find('delins') == -1 and \
                   other_Name_data.find('_') == -1:
                    clin_var_data_dict[accession_num] = []

        for line in lines:
            accession_num = line[2].split('(')[0]
            if accession_num.find('NM') >= 0:
                # ANS-NK:Replacing first parenthese with star to find amino info
                line[2] = line[2].replace('(','*',1)

                current_line = line[2].split('(')

                amino_info = current_line[-1]
                original_amino = ''
                position = ''
                mutant_amino = ''
                deletion_y_n = ''

                # ANS-NK:CASES that distinguish how mutation data is parsed
                # ANS-NK:Adjusts original_amino, position, mutant_amino,
		# deletion_y_n values

                # ANS-NK:CASE: if no protein data is given
                if amino_info.find('p.') == -1:
                    #print 'noprotdata', amino_info
                    pass
                # ANS-NK:CASE: if the mutation is a large deletion and insertion
                elif amino_info.find('delins') >= 0:
                    pass
                # ANS-NK:CASE: if the nucleotide mutation is the same amino acid
                elif amino_info.find('=') >= 0:
                    original_amino = amino_info[2:5]
                    for char in amino_info:
                        if char.isdigit():
                            position += char
                    position = position
                    mutant_amino = original_amino
                    deletion_y_n = 'no'
                # ANS-NK:CASE: if the mutation is a deletion
                elif amino_info.find('del') >= 0:
                    original_amino = amino_info[2:5]
                    for char in amino_info:
                        if char.isdigit():
                            position += char
                    position = position
                    mutant_amino = amino_info[(len(amino_info)-4):-1]
                    deletion_y_n = 'yes'
                # ANS-NK:CASE: if the mutation is a frameshift
                elif amino_info.find('fs') >= 0:
                    original_amino = amino_info[2:5]
                    for char in amino_info:
                        if char.isdigit():
                            position += char
                    position = position
                    mutant_amino = amino_info[(len(amino_info)-6):-1]
                    deletion_y_n = 'no'
                # ANS-NK:CASE: if the mutation leads to a termination (Ter)
                elif amino_info.find('Ter') >= 0:
                    original_amino = amino_info[2:5]
                    for char in amino_info:
                        if char.isdigit():
                            position += char
                    position = position
                    mutant_amino = amino_info[(len(amino_info)-4):-1]
                    deletion_y_n = 'no'
                # ANS-NK:CASE: if the mutation is a simple substitution
                else:
                    original_amino = amino_info[2:5]
                    for char in amino_info:
                        if char.isdigit():
                            position += char
                    position = position
                    mutant_amino = amino_info[(len(amino_info)-4):-1]
                    deletion_y_n = 'no'

                # ANS-NK:Appends the relevant information to the entry for each
                # accession number
                other_Name_data = line[2].split('(')[-1]
                if other_Name_data.find('p.') >= 0 and \
                   other_Name_data.find('delins') == -1 and \
                   other_Name_data.find('_') == -1:
                    clin_var_data_dict[accession_num].append([original_amino,\
                        position, mutant_amino, deletion_y_n, line[1],line[6],\
                        line[13], line[18]])


        # ANS-NK:Removes duplicates
        for accession_num in clin_var_data_dict:
            copy = []
            for entry in clin_var_data_dict[accession_num]:
                if entry not in copy:
                    copy.append(entry)
            clin_var_data_dict[accession_num] = copy

        cleanup_dict(clin_var_data_dict)

        # ANS-NK:Prints all keys and values in dict
        #keys = clin_var_data_dict.keys()
        #values = clin_var_data_dict.values()
        #for i in range(len(keys)):
                #print keys[i], '\t', values[i]

    return clin_var_data_dict
def convert_to_string(align):
    with open('aligned.aln', 'r') as aligned_seq:
        header = ''
        symbols = ''
        for line in aligned_seq:
            if line[0] == 'M':
                header = header + line
            if line[0] == ' ':
                line = line.strip('\n')
                symbols = symbols + line
    return symbols

def conserve_CVs(clin_var_data_dict, symbols): # Katie Greene, Andrew Peterson, Anna Catalanotto
    '''
    Identify Conserved Clinical Variants Between Protein Sequences

    This method identifies the clinical variants conserved between the non-human
    protein sequence and the human protein sequence. The method removes the gaps
    from the align_ortholog function's alignment string. Creates a list of the cv
    index positions. Determines if matches are conserved clinical variants and formats
    a new alignment string. Identify the conserved clinical variants and their parameters
    to generate data for a conserved clinical variant table.

    Args:
    var_output = nested list of clinical variants and parameters
    alignment = list with non-human seq, alignment string from align_ortholog, human_seq,
                alignment score, and start/end aa index

    Return:
    Conserved CV List-
    alignment = list with non-human seq, CONSERVED alignment string, human_seq,
                alignment score, and start/end aa index
    conserved_var_output = nested list of CONSERVED clinical variants and parameters


    '''
    # AP, KG, AC: Remove Spaces from Alignment String
    pure_alignment = ''

    for symbol in symbols:
        if symbol is not " ":
            pure_alignment = pure_alignment + symbol
        else:
            pass # AP, KG, AC: If symbol is space, do not add to alignment string



    # AP, KG, AC: Create List of All Human CV Index Positions
    cv_aa_positions = [] # AP, KG, AC: list of cv aa positions (strings)
    for nest in clin_var_data_dict.values():
        for inner in nest:
            cv_aa_positions.append(inner[1])


    # AP, KG, AC: Identify Conserved Clinical Variants & Format Alignment String
    index = 0 # AP, KG, AC: Counter for index position in pure_alignment
    cv_alignment_str = ''
    print(pure_alignment)
    for symbol in pure_alignment:
        if symbol == ":":
            index_pos = pure_alignment.find(":", index, len(pure_alignment))
            # AP, KG, AC: find the first aa match from index position
            aa_pos = index_pos + 1 # AP, KG, AC: convert index position to aa position
            aa_pos_str = str(aa_pos) # AP, KG, AC: convert aa to string to search var_output
            if aa_pos_str in cv_aa_positions: # AP, KG, AC: If aa position in cv index list
                cv_alignment_str = cv_alignment_str + "|"
            else: # AP, KG, AC: Else, conserve the alignment symbol from align_ortho
                cv_alignment_str = cv_alignment_str + symbol
            index = index + 1
        else:  # AP, KG, AC: Else, conserve the alignment symbol from align_ortho
            cv_alignment_str = cv_alignment_str + symbol
            index = index + 1

    # AP, KG, AC: Identifying Conserved Clinical Variant and Parameters (Table)
    conserved_cv_index = [] # AP, KG, AC: List of Conserved CV Index positions
    start = 0
    # print('cas', cv_alignment_str)
    for symbol in cv_alignment_str:
        if symbol == "|":
            #print('i am a |')
            index = cv_alignment_str.find("|", start) # AP, KG, AC: Find index position of conserved clinical variant
            start = index + 1
            index = str(index+1)
            conserved_cv_index.append(index) # AP, KG, AC: Append conserved CV index position to list
        if symbol == ':':
            #print('i am a :')
            index = cv_alignment_str.find(':', start)
            start = index + 1
            index = str(index+1)
        else:
            #print('i am nothing')
            start = start + 1 # # AP, KG, AC: Iterate through cv_alignment_str

    # print('ccv', conserved_cv_index)
    conserved_var_output = [] # # AP, KG, AC: List of Conserved CVs w/ Parameters (for Table)
    for index in conserved_cv_index: # For every conserved index...
        for cv in clin_var_data_dict: # and for every cv in var_output from id_cv...
            #print(cv)
            if index == cv[0]: # if conserved cv index is present in the table
                conserved_var_output.append(cv) # append conserved cvs to table list
            else:
                pass


    # AP, KG, AC: Append Convserved CVs in Alignment String
    complete_cv_alignment_str = ''
    # for symbol in align_seqs()[1]:
    #     if symbol == " ":
    #         complete_cv_alignment_str = complete_cv_alignment_str + ' '
    #     else:
    #         complete_cv_alignment_str = complete_cv_alignment_str + cv_alignment_str[0]
    #         cv_alignment_str = cv_alignment_str[1:]
    # align_seqs()[1] = align_seqs()[1] + complete_cv_alignment_str

    #print('align', align_seqs())
    #print('cvo', conserved_var_output)
    #return [align, conserved_var_output] # AP, KG, AC: list containing String of cosnerved Alignment symbols and conserved clincial variants

def main():
    # print (align)
    clin_var_data_dict = clin_var_parser()
    symbols = convert_to_string(align)
    conserve_CVs(clin_var_data_dict, symbols)

if __name__ == '__main__':
    main()
