import subprocess
import sys
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO

#takes in the test FASTA file and aligns all of the sequences within the file. This includes both the human and nonhuman proteins.
out_file = "aligned.fasta"
muscle_cline=MuscleCommandline(input="NR5A_orthologs.fasta", out=out_file)
stdout, stderr = muscle_cline()
#align=AlignIO.read(StringIO(stdout), "fasta")
align=AlignIO.read(out_file, "fasta")
print(align)


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

    print clin_var_data_dict

def main():
    print (align)
    clin_var_parser()

if __name__ == '__main__':
    main()
