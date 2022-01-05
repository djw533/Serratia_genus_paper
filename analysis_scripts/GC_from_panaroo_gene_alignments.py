#! /usr/env python


from Bio.SeqUtils import GC
import multiprocessing as mp
from multiprocessing import Process
import multiprocessing
import tqdm




def parseArgs():
    """Parse command line arguments"""

    import argparse
    try:
        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter ,
            description="Get GC content for codon positions 1/2/3 and whole gene sequence for multiple alignments of all genes in panaroo")
        parser.add_argument('-i',
    		    '--input',
    		    action='store',
                nargs='+',
                #help='<Required> Set flag',
                required=True,
    		    help='fasta alignment files as input')
        parser.add_argument('-o',
    		    '--output',
    		    action='store',
                #help='<Required> Set flag',
                required=True,
    		    help='output csv file')
        parser.add_argument('-t',
    		    '--threads',
    		    action='store',
                #help='<Required> Set flag',
                default = 1,
    		    help='number of threads')


    except:
        print("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()

    return parser.parse_args()


def __read_alignment_and_remove_gaps__(aln_file):

    # open the data file for the alignment:
    data = open(aln_file)
    lines = data.readlines()
    data.close()

    sequences = {}

    counter = 0
    for line in lines:
        new_sequence = False # set false setting for whether there is a new sequence in the file:
        if line.startswith(">"):
            new_sequence = True
            header_name = line.strip().split('>')[-1].split(';')[0] # only want to get the strain name, not the gene number! - as the output in panaroo is formatted

            if counter == 0:
                counter += 1
                temp_sequence = []
                continue
            else:
                #now write out the temp sequence as a string
                counter += 1
                sequences[header_name] = ''.join(temp_sequence)

            ##now reset the temp_sequence list

            temp_sequence = []


        # if new_sequence == True:
        #     sequences[header_name] = [] # create new list for this sequence:
        #     continue # now continue to next line

        if new_sequence == False:
            temp_sequence.append(line.strip().replace("-","")) #append sequence, getting rid of the new line character and any gaps:


    #return the dictionary:
    return sequences


def __get_GC_content_from_panaroo_alignment__(sequences):


    GC_content_details = {}

    for strain, sequence in sequences.items():

        #get GC content for all:
        GC_full = GC(sequence)

        ### now create strings with the bases in positions 1/2/3

        number_of_codons = int(len(sequence)/3)

        GC1_seq = []
        GC2_seq = []
        GC3_seq = []

        for i in range(number_of_codons):
            GC1_seq.append(sequence[(i*3) + 0])
            GC2_seq.append(sequence[(i*3) + 1])
            GC3_seq.append(sequence[(i*3) + 2])

        ### now get the GC contents:
        GC1 = GC(''.join(GC1_seq))
        GC2 = GC(''.join(GC2_seq))
        GC3 = GC(''.join(GC3_seq))

        ### add this information to the dictionary
        GC_content_details[strain] = {}
        GC_content_details[strain]["GC1"] = GC1
        GC_content_details[strain]["GC2"] = GC2
        GC_content_details[strain]["GC3"] = GC3
        GC_content_details[strain]["GC_full"] = GC_full



    return(GC_content_details)


def main():

    args = parseArgs()


    # get names of all alignments:
    alignments = args.input

    # create dictionary to write all into
    GC_results = {}

    # A - loop through all alignment
    for aln in alignments:

        gene_name = aln.split('/')[-1].replace(".aln.fas","")

        # get all sequences and remove gaps:
        sequences = __read_alignment_and_remove_gaps__(aln) # sequences is a directory - keys are strains, values are

        ## get GC content values for each group:
        GC_results[gene_name] = __get_GC_content_from_panaroo_alignment__(sequences)


    #write out:

    with open(args.output,"w") as output:
        output.write("strain,gene,GC_gene,GC1,GC2,GC3\n")

    ### now go through the GC results:

    for gene_name, sequences in GC_results.items():
        for key in sequences:
            with open(args.output,"a") as output:
                output.write("{strain},{gene},{GC_gene},{GC1},{GC2},{GC3}\n".format(strain = key, gene = gene_name, GC_gene = sequences[key]["GC_full"], GC1 = sequences[key]["GC1"], GC2 = sequences[key]["GC2"], GC3 = sequences[key]["GC3"]))


if __name__ == '__main__':
    main()
