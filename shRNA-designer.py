#
# shRNA Design tool for mir30 based backbones
# by Maximilian Pfisterer (c) 2022
# based on data by Matveeva et al., 2012
#
# This script processes an mRNA sequence by dividing it into chunks of
# specified lengths. The chunks are scored by correlation values for
# each nucleotide position with the efficacy of the sequence as an
# shRNA.
#
# This script should give similar results to the miR_Scan program of
# the authors of the paper but allows for more customization and easier
# export of shRNA as ready to use oligos.
#


import os
import random
import time
import pickle
import argparse as ap
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqUtils import gc_fraction


# Open files needed for the scoring
fd = os.path.dirname(os.path.realpath(__file__))+"/"
coeffs = pd.read_pickle(fd+"coeffs.pkl")
motifs_to_avoid = ["TTTT","AAAA","GGGG","UGGC","CGGA"]
deltaG_bins = pd.read_pickle(fd+"dG_bins.pkl")

# From Xia et al., 1998
with open(fd+"deltaG.pkl", "rb") as handle:
    deltaG = pickle.load(handle)

# Configure all arguments
argparser = ap.ArgumentParser(
    description="shRNA Design Tool for mir30 based shRNA backbones.",
    epilog="""
    Based on work by Matveeva et al., 2012. No guarantees.
    
    This program processes a FASTA file with mRNA sequences and ranks all
    possible shRNAs for their efficiency in a mir30a based backbone. shRNAs are
    scored by the effect of nucleotides at specific positions on the targeting
    efficiency, observed in empirical studies. The shRNAs are ranked from
    -1.306 to 1.336. Additionally duplex stability for each dinucleotide is
    calculated and also scored for effects on silencing efficiency observed in
    empirical studies.
    """)
argparser.add_argument(
    "inputfile",
    metavar = "PATH", nargs = "?",
    default = "-",
    help = "Path to FASTA File")
argparser.add_argument(
    "-l", "--length",
    metavar = "nt", dest = "length",
    type = int, default = 22,
    help = "Length of the guide sequence, (default 22)")
argparser.add_argument(
    "--gc",
    metavar = "%", dest = "gc",
    type = int, default = [0,100], nargs = 2,
    help = "Minimum and maximum GC content, (default 0, 100)")
argparser.add_argument(
    "-s", "--min-score",
    metavar = "", dest = "min_score",
    type = float, default = -9999,
    help = "Minimum score")
argparser.add_argument(
    "--sort-dg", dest = "score_dg",
    action = "store_true",
    help = "Sort shRNAs by ∆G-Score")
argparser.add_argument(
    "--remove-motifs", dest = "removemotifs",
    action = "store_true",
    help = "Remove shRNAs with low complexity motifs (e.g. TTTT)")
argparser.add_argument(
    "--multitarget", dest = "multitarget",
    type = int, default = 1,
    help = "Rquires shRNAs to target at least n targets from input file")
argparser.add_argument(
    "--guide-on-5p", dest = "guide_on_5p",
    action = "store_true",
    help = "Put the guide strand on the 5' arm of the hairpin")
argparser.add_argument(
    "--oligo",
    metavar = "NNN", dest = "oligo",
    nargs = 3,
    help = "Sequences for creating Oligo (5p-flank, loop, 3p-flank")
argparser.add_argument(
    "--annealing",
    metavar = "Nt", dest = "annealing",
    nargs = 2, type = int,
    # [0] = left (5') overhang [1] = right (3') overhang
    help = "Create forward and reverse oligo for  annealing, specify  5' and \
            3' overhang size")
argparser.add_argument(
    "--format",
    default = ["xlsx"], dest = "format",
    nargs = "*",
    help = "Output format (default xlsx)")
argparser.add_argument(
    "-o", "--outfile",
    metavar = "PATH",
    default = "~/shRNAs",
    help = "Output filename, without suffix")
argparser.add_argument(
    "-gff", "--gff3",
    dest = "verbose", action = "store_true",
    help = "Output shRNA positions as GFF3 file")
argparser.add_argument(
    "--rna",
    dest = "rna", action = "store_true",
    help = "Output all sequences as RNA")
argparser.add_argument(
    "--sense",
    dest = "sense", action = "store_true",
    help = "Output targeting sequences as sense strand")
argparser.add_argument(
    "--stemmismatch",
    default = 1, type = int, nargs = "+",
    help = "Introduce mismatches at positions on passenger strand N [N N N] (default 1, requires --oligo)")
argparser.add_argument(
    "-v", "--verbose",
    dest = "verbose", action = "store_true",
    help = "Be verbose, for debugging")

args = argparser.parse_args()

def process_arguments():
    if args.annealing:
        if args.annealing[0] > len(args.oligo[0]) \
        or args.annealing[1] > len(args.oligo[2]):
            print("WARNING: The overhang size you specified is greater than"
                  "the oligo flanking \nsequences")
    if args.multitarget > len(mrnas):
        args.multitarget = len(mrnas)
        print("WARNING: Multitarget was set higher than the number of \n"
              "sequences provided. multitarget was set to the maxmimum",
                len(mrnas))

def verbose_out(list_of_strings):
    if args.verbose:
        string_to_print = " ".join(str(e) for e in list_of_strings)
        print(string_to_print)

# Creates a list of all possible shRNA sequences with a rolling frame
# over mRNA of args.length.
def make_shrnas(list_of_mrnas):
    start = time.time()
    pre_shrnas = {}
    for mrna in list_of_mrnas:
        verbose_out(["Processing:", mrna.id])
        pre_shrnas[mrna.id] = list(mrna.seq.upper()
                              .reverse_complement()[0+i:args.length+i]
                              for i in range(0, len(mrna)-args.length, 1))
        verbose_out([len(pre_shrnas[mrna.id]), "possible shRNA sequences for",
                     mrna.id])
    end = time.time()
    verbose_out(["make_shrna took [s]", (end-start)])
    return pre_shrnas

def scoreBins(position, dG_NN):
    if dG_NN > -1:
        return deltaG_bins.iloc[position, 0]
    elif dG_NN > -2:
        return deltaG_bins.iloc[position, 1]
    elif dG_NN > -3:
        return deltaG_bins.iloc[position, 2]
    else:
        return deltaG_bins.iloc[position, 3]

# Score possible shRNAs with the coefficients of correlation from the
# publication. Positive correlation coeff correspond to a positive
# effect on shRNA effectiveness, and vice versa. The values are added
# for each position to calculate a score.
def process_shrnas(list_of_pre_shrnas):
    start = time.time()
    scored_shrnas = {}
    
    for mrna in list_of_pre_shrnas:
        new_scored_shrnas = []
        i = 0
        for shrna in list_of_pre_shrnas[mrna]:
            i += 1
            # Filter out sequences by their GC content
            if not (args.gc[0]< gc_fraction(shrna) < args.gc[1]):
                print("filtered by GC content")
                continue
            

            # Calculate Score
            score = 0
            for position in range(0, len(shrna)):
                score += coeffs.loc[position, shrna[position]]
            
            # Score after deltaG bins
            list_of_dinucs = list(Seq(shrna).transcribe()[0+i:2+i] for i in range(0, len(shrna)-2, 1))
            dG = 0.0
            dG_score = 0.0
            pos = 0
            # If the dinucleotide is not found in dict, use the reverse

            for dinuc in list_of_dinucs:
                if dinuc in deltaG:
                    dG_NN = deltaG[dinuc]
                    dG += dG_NN
                    dG_score += scoreBins(pos, dG_NN)
                else:
                    dG_NN = deltaG[dinuc.reverse_complement().transcribe()]
                    dG += dG_NN
                    dG_score += scoreBins(pos, dG_NN)
                pos += 1
            dG += 4.09
            
            # Filter out sequences with too low score
            if score < args.min_score:
                print("filtered by score")
                continue

            # Filter out sequences with low complexity motifs
            if args.removemotifs:
                if any(motif in shrna for motif in motifs_to_avoid):
                    print("filtered by motif")
                    continue

            # Filter out sequences by their number of targets
            targets = 0
            for mrna_ in mrnas:
                if mrna_.seq.reverse_complement().count(shrna) > 0:
                    targets += 1
            if targets < args.multitarget:
                print("filtered by target")
                continue

            fwd_strand = ""
            rev_strand = ""
            full_stem = ""
            if args.oligo:
                # Creates a hairpin sequence
                guide = shrna
                passenger = shrna.reverse_complement()
                if args.stemmismatch:
                    if args.stemmismatch != 0:
                        if type(args.stemmismatch) == int:
                            args.stemmismatch = [args.stemmismatch]
                        passenger = MutableSeq(passenger)
                        for mm_position in args.stemmismatch:
                            #print(passenger[mm_position-1])
                            if passenger[mm_position-1] == "G": #C on guide
                                list_of_nucleotides = ["T","C","A"]
                            if passenger[mm_position-1] == "C": #G on guide
                                list_of_nucleotides = ["C","A","G"]
                            if passenger[mm_position-1] == "T": #A on guide
                                list_of_nucleotides = ["A","C", "G"]
                            if passenger[mm_position-1] == "A": #U on guide
                                list_of_nucleotides = ["C","T"]
                            #list_of_nucleotides.remove(passenger[mm_position-1])

                            passenger[mm_position-1] = "".join(random.choices(list_of_nucleotides))
                full_stem = args.oligo[0] \
                             + (guide if args.guide_on_5p else passenger) \
                             + args.oligo[1] \
                             + (passenger if args.guide_on_5p else guide) \
                             + args.oligo[2]

                if args.annealing:
                    # Create the second oligo for annealing
                    rev_strand = Seq(full_stem).reverse_complement()[:-args.annealing[0]]
                    fwd_strand = full_stem[:-args.annealing[1]]
            new_scored_shrnas.append({
                "name":os.path.splitext(args.inputfile.split(os.sep)[-1])[0]+"-sh"+str(i) if args.multitarget == len(mrnas) else mrna+"-sh"+str(i),
                "start_pos":len(list_of_pre_shrnas[mrna])-i,
                "antisense": shrna.transcribe() if args.rna else shrna,
                "sense": shrna.reverse_complement().transcribe() if args.rna else shrna.reverse_complement(),
                "Score":round(score,4),
                "dG_Score":round(dG_score,1),
                "GC%":round(gc_fraction(shrna),1),
                "dG":round(dG,1),
                "Targets":"/".join(str(e) for e in [targets, len(mrnas)]),
                "insert":full_stem,
                "fwd_strand":fwd_strand,
                "rev_strand":rev_strand})
        if len(new_scored_shrnas) > 0:
            scored_shrnas[mrna] = pd.DataFrame(new_scored_shrnas) \
                                  .sort_values(by='Score' if not args.score_dg else 'dG_Score', ascending=False)
        else:
            scored_shrnas[mrna] = pd.DataFrame()

    end = time.time()
    verbose_out(["score_shrnas took [s]", (end-start)])
    return scored_shrnas

def export_shrnas(list_of_shrnas):
    for format_type in args.format:
        if format_type == "csv":
            for mrna in list_of_shrnas:
                outfilepath = "".join([os.path.splitext(args.outfile)[0],".csv"]) if args.multitarget >= len(mrnas) else "".join([os.path.dirname(args.outfile), mrna,".csv"])
                list_of_shrnas[mrna].to_csv(outfilepath)
    
        if format_type == "xlsx":
            # Creates an excel file with 1 sheet for each mRNA
            excel_writer = pd.ExcelWriter("".join([os.path.splitext(args.outfile)[0], ".xlsx"]),
                                         engine='xlsxwriter')
            for mrna in list_of_shrnas:
                # Sheet names have to be limited to 31 chars
                list_of_shrnas[mrna].to_excel(excel_writer, sheet_name=mrna[0:31])
                worksheet =  excel_writer.sheets[mrna[0:31]]
                # Sets column widths to optimal length
                for idx, col in enumerate(list_of_shrnas[mrna]):
                    series = list_of_shrnas[mrna][col]
                    
                    # Determine needed width of column
                    # Size of longest content or header row, whichever is greater
                    max_len = max((
                        series.astype(str).map(len).max(),  # len of largest item
                        len(str(series.name))  # len of column name/header
                        )) + 1
                    
                    # Set column width, +1 since the first row in excel
                    # sheet is not contained in dataframe as it is the row index.
                    worksheet.set_column(idx+1, idx+1, max_len)

                worksheet.conditional_format(
                    "F2:F"+str(len(list_of_shrnas[mrna])+1),
                    {"type":"3_color_scale", "min_value":-1.306, "max_value":1.336, "mid_value":0, "min_type":"num", "mid_type":"num", "max_type":"num"})
                worksheet.conditional_format(
                    "G2:G"+str(len(list_of_shrnas[mrna])+1),
                    {"type":"3_color_scale", "min_value":17.28, "max_value":46.5, "mid_value":31.89, "min_type":"num", "mid_type":"num", "max_type":"num"})
                worksheet.conditional_format(
                    "H2:H"+str(len(list_of_shrnas[mrna])+1),
                    {"type":"3_color_scale", "min_color":"#63BE7B", "max_color":"#F8696B", "min_value":0, "max_value":100, "mid_value":50, "min_type":"num", "mid_type":"num", "max_type":"num"})
                worksheet.conditional_format(
                    "I2:I"+str(len(list_of_shrnas[mrna])+1),
                    {"type":"3_color_scale", "min_color":"#F8696B", "max_color":"#F8696B", "mid_color":"#63BE7B", "min_value":-40, "mid_value":-33, "max_value":-20, "min_type":"num", "mid_type":"num", "max_type":"num"})
                
                # When only shRNAs that target all sequences are
                # exported, all sheets would be identical, so only one
                # is saved.
                if args.multitarget == len(mrnas):
                    break
                    
            excel_writer.close()

with open(args.inputfile) as filehandle:
    mrnas = list(SeqIO.parse(filehandle, "fasta"))
process_arguments()
processed_shrnas = process_shrnas(make_shrnas(mrnas))
export_shrnas(processed_shrnas)
