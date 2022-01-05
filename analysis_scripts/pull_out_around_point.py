#! /usr/env python


### argparse

import sys
import os
import traceback
import Bio
import shutil
import re
from Bio import SeqIO

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import translate
from Bio.Seq import reverse_complement
from shutil import copyfile
import numpy as np
from datetime import datetime
from Bio.SeqUtils import GC
import time
from datetime import timedelta

def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter ,
            description="Parse assmebly file for stats and fasta sequence")
        parser.add_argument('-i',
                '--input',
                action='store',
                required=True,
                help='assembly and locus_tag file - with a header for strain and locus_tag - in that order: - dont include .gff extension in the assembly name')
        parser.add_argument('-g',
                '--gffs',
                action='store',
                required=True,
                help='gff directory:')

        parser.add_argument('-o',
			    '--outdir',
			    action='store',
                default="extracted_sequences",
			    help='ourput directory')

        parser.add_argument('-u',
			    '--upstream',
                type = int,
			    action='store',
                default=0,
			    help='Number of nucleotides to include upstream/"left" of gene cluster, default = 0')

        parser.add_argument('-d',
			    '--downstream',
                type = int,
			    action='store',
                default=0,
			    help='Number of nucleotides to include downstream/"right" of gene cluster, default = 0')

    except:
        print("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()

    return parser.parse_args()



##functions

#### extracting from 2 indexes function   -  editing as if there is no upstream of downstream, then the start region is always missing one base...
def extract_a2b(start_index, stop_index, gff_input, strain, contig, upstream, downstream, outname):
    """ Take sequence and annotation from gff3 files between genomic point A and B, given both base indices """

    #### get cwd
    currentdir = os.getcwd()+'/'



    #### first check whether the contig is long enough

    record = SeqIO.read(currentdir+"only_fasta/singlefastas/"+contig+'.fasta', "fasta")
    contig_length = len(record)

    #print 'Checking to see if hit is on the positive or negative strand: -----'
    ### N.B this should always be positive - this is a bit that has been taken over from a previous script that equated blast hits (which had been reversed by blast)
    ### to gff subsections from the original annotated genome
    result = stop_index - start_index
    if result > 0:
        #print 'Hit is on the positive strand, proceeding : ------'
        start = int(start_index) - int(upstream)
        stop = int(stop_index) + int(downstream)

        #check if the contig is long enough:
        if start < 1 :
            start = 1
        if stop > contig_length:
            stop = contig_length

    elif result < 0:
        #print 'Hit is on the negative strand, redefining the start and stop indices'
        start = int(stop_index) - int(downstream)
        stop = int(start_index) + int(upstream)

        #again check if the contig is long enough: (and correct if not):
        if start < 1 :
            start = 1
        if stop > contig_length:
            stop = contig_length

    elif result == 0:
        print('Start and Stop are the same, please check input')

    ## -- loop through each gff file in the list made above and start reading from the 2nd line of the gff file
    ## -- create new file and insert header for a new gff file

    #for gff in gff_list:

    gff_file = open(gff_input)
    gff_file_lines = gff_file.readlines()
    gff_file.close

    with open(outname, "w") as output:
        output.write("##gff-version 3")
        output.write("\n")

   # print 'Beginning extraction from ', strain


    ## -- index the line where the fasta sequence begins, and make a list of all gff features above that line
    ## -- i.e. all annotation

    for line in gff_file_lines:
        if not line.startswith("#"):
            line_to_index = line
            break

    #print line_to_index

    fasta_target = gff_file_lines.index('##FASTA\n')
    hash_index = gff_file_lines.index(line_to_index)
    gff_features = gff_file_lines[hash_index:fasta_target]

    #print "fasta line is: ", fasta_target


    ## -- Next line of gff file written detailing the sequence region, how long the region is, and what the name of the chromosome is in the associated fasta sequence
    with open(outname, "a") as output:
        output.write('##sequence-region'+' '+contig+' 1 '+str(stop-(start+1)))
        output.write("\n")

    ## -- Loop through annotation - define each feature start and stop indices
    ## -- If feature start and stop sites are between vgrG1_1 and vgrG1_2 (start and stop),
    ## -- annotation line is appended to new_list as a new item
    ## -- Feature start and stop indices are corrected so as to refer to the new sequence length

    for gff_line in gff_file_lines[hash_index:(fasta_target)]:
        if gff_line.startswith('#'):
            continue
        else:
            feature_start = gff_line.split('\t')[3]
            feature_stop = gff_line.split('\t')[4]

        if int(feature_start) >= start and int(feature_stop) <= stop and contig in gff_line:
            new_list = gff_line.split("\t")[0:3] + gff_line.split("\t")[5:8]
            new_list.append(" ".join(gff_line.strip().split("\t")[8:]))
            new_list.insert(3, str(int(feature_start) - start))
            new_list.insert(4, str(int(feature_stop) - start))
            linegff3 = "\t".join(new_list)
            with open(outname, "a") as output:
                output.write(linegff3)
                output.write("\n") ####### have removed the newline character "\n"

        else:
            pass

    ## -- Append ##fasta header to the annotation to indicate that the fasta sequence will be below

    with open(outname, "a") as output:
        output.write('##FASTA\n')


    ## -- Append fasta sequence to the end of the gff file

    record = SeqIO.read(currentdir+"only_fasta/singlefastas/"+contig+'.fasta', "fasta")
    with open(outname, "a") as out:
        SeqIO.write(record[start:stop], out, "fasta")

    # with open(outname+".fasta", "w") as out:
    #     SeqIO.write(record[start:stop], out, "fasta")


    ## -- Print counter to terminal screen to show progress

   # print 'Sequence and annotation extracted from', strain

def __get_protein_info__(assembly, gene, assembly_dir):

    data = open(assembly_dir + "/" + assembly + ".gff")
    lines = data.readlines()
    data.close()

    print(gene)


    fasta = False
    for line in lines:
        if line.startswith("##FASTA"):
            print(line) # if got to here - then the gene isn't present in the annotation
            break
        if line.startswith("#"):
            continue

        ### split the line

        toks = line.strip().split('\t')
        if "ID=" in line:

            gene_name = toks[8].split('ID=')[1].split(';')[0]
            if gene_name == gene:
                contig = toks[0]
                start = toks[3]
                stop = toks[4]
                strand = toks[6]
                break

    output_tuple = (gene_name,contig,start,stop,strand)
    return output_tuple



#script

# take input gfa and get the fasta sequence:
#
# also get the stats - legnth , sequence depth - circular or not - links

def __get_fastas_from_gffs(gff,strain):
    ### extract fasta sequence and annotation from gff file
    currentdir = os.getcwd()+'/'

    fasta_data = open(gff)
    gff_lines_for_fasta = fasta_data.readlines()
    number_for_parsing_from = gff_lines_for_fasta.index("##FASTA\n")
    fasta_lines = gff_lines_for_fasta[number_for_parsing_from+1:]
    annotation_lines = gff_lines_for_fasta[:number_for_parsing_from]

    ##fasta
    if not os.path.exists(currentdir+"only_fasta"):
        os.makedirs(currentdir+"only_fasta")
        open(currentdir+"only_fasta/"+strain+".fasta","w")
    with open(currentdir+"only_fasta/"+strain+".fasta","a") as output:
        for line in fasta_lines:
            output.write(line)

    ##annotation
    # if not os.path.exists(currentdir+"only_annotation"):
    #     os.makedirs(currentdir+"only_annotation")
    #     open(currentdir+"only_annotation/"+strain+".gff","w")
    # with open(currentdir+"only_annotation/"+strain+".gff","a") as output:
    #     for line in annotation_lines:
    #         output.write(line)


    ###should also make this a function!!!!!
    ### create singlefastas to be read for each operon
    if not os.path.exists(currentdir+"only_fasta/singlefastas"):
      os.makedirs(currentdir+"only_fasta/singlefastas")

    data = open("only_fasta/"+strain+".fasta")
    lines = data.readlines()

    lines_string = ''.join(lines)
    #print lines_string
    #print GC(lines_string)
    no_entries = [x for x in lines if not x.startswith('>')]
    no_entries_string = ''.join(no_entries).replace('\n','')
    GC_genome = GC(no_entries_string)

    entries_list = (''.join(lines)).split('>')

    ### make dictionary for the contig lengths:
    contig_lengths = {}

    ##make dictionary for protein hits results:
    hits_results = {}

    for entry in entries_list:
         if entry != "":
             name = entry.split('\n')[0].split()[0]
             contig_lengths[name] = len(''.join(entry.split('\n')[1:]))
             #print name
             with open(currentdir+"only_fasta/singlefastas/"+name+'.fasta', "w") as output:
                  output.write('>'+entry)


##

def main():

    args = parseArgs()


    #create outdir :
    if os.path.isdir(args.outdir):
        sys.exit(args.outdir+" is already a directory! Exiting script")
    else:
        os.makedirs(args.outdir)


    assemblys_and_genes_data = open(args.input)
    assemblys_and_genes = assemblys_and_genes_data.readlines()
    assemblys_and_genes_data.close()

    assemblies_parsed = {}

    for input_data in assemblys_and_genes[1:]:
        names_toks = input_data.strip().split(',')
        assembly_name = names_toks[0]
        print(assembly_name)
        gene_name = names_toks[1]

        ##extracted something from this assembly already? - check with dict
        if assembly_name not in assemblies_parsed:
            assemblies_parsed[assembly_name] = 0
        elif assembly_name in assemblies_parsed:
            assemblies_parsed[assembly_name] += 1

        ### extract fasta and singlefastas:
        __get_fastas_from_gffs(args.gffs+'/'+assembly_name+".gff",assembly_name)

        gene_position_data = __get_protein_info__(assembly_name, gene_name, args.gffs)

        gene_name_extracted = gene_position_data[0]

        if gene_name_extracted != gene_name: # check that the right gene was found:
            print("names not the same")

        contig = gene_position_data[1]
        start = gene_position_data[2]
        stop = gene_position_data[3]
        strand = gene_position_data[4]


        # create the output name:

        outname = "{outdir}/{assembly_name}_cluster_{num}.gff".format(outdir = args.outdir, assembly_name = assembly_name, num = str(assemblies_parsed[assembly_name] + 1))


        ### now extract this! - with 30kb upstream and downstream
        extract_a2b(int(start), int(stop), args.gffs+'/'+assembly_name+".gff", assembly_name, contig, args.upstream, args.downstream, outname)

        ###clean up:
        currentdir = os.getcwd()+'/'
        shutil.rmtree(currentdir+"only_fasta")




if __name__ == '__main__':
    main()
