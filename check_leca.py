#!/usr/bin/env python3
import sys
from eukarya import *
from scrollsaw import *
import argparse

# ------

parser = argparse.ArgumentParser(description = "This script checks if a set of sequences fulfills the LECA criterion by using the corresponding BLAST hits.")
parser.add_argument("pfam", metavar = "Pfam", help = "the Pfam for which you want to check if it was in LECA")
parser.add_argument("sequence_ids", metavar = "sequence_ID", nargs = '+', help = 'the sequences from the BBHs')
parser.add_argument("-l", metavar = "0.xx", help = "coverage threshold for LECA calling (DEFAULT: 0.15)", type = float, default = 0.15)
parser.add_argument("-b", metavar = "BLASTdir", help = "directory containing the BLAST output files and the file containing all sequence IDs (DEFAULT: current)", default = '.')
parser.add_argument('-o', metavar = 'outdir', help = 'directory for output files (DEFAULT: current)', default = '.')
args = parser.parse_args()

prefix = args.pfam
bbh_seqs = args.sequence_ids
coverage_criterion = args.l
blastdir = args.b
outdir = args.o

with open(f'{blastdir}/{prefix}_seqids.list') as all_seqs_file:
    all_seqs = [line.rstrip() for line in all_seqs_file]
representing = assign_all_seqs(bbh_seqs, all_seqs, euk_only = True, prefix = prefix, blast_path = blastdir)
coverage, copy_no, species = infer_coverage_redundancy(bbh_seqs, representing, supergroups5, tree = False)
if coverage <= coverage_criterion:
    sys.exit(f'No LECA in this Pfam. Coverage is {coverage}.')
else:
    with open(f'{outdir}/{prefix}_lecas.tsv', 'w') as lecas_out:
        lecas_out.write('Pfam\tFECA\tAncestry\tLECA\tSupport\tCoverage\tCopy number\tSeqs\n')
        print(prefix, 'NA', 'Eukaryotic', 'OG1.1', 'NA', coverage, copy_no, ','.join(bbh_seqs), sep = '\t', file = lecas_out)
    print('Pfam\tFECAs (normal)\tFECAs (strict)\tFECAs (after merging)\tLECAs\tUnknowns\tNon-FECAs')
    print(prefix, 'NA', 'NA', 'NA', 1, 0, 'NA', sep = '\t')
