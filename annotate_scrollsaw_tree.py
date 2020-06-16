#!/usr/bin/env python

# Load modules
from ete3 import TreeStyle
from feca2leca import *
import os
import argparse
ts = TreeStyle()

#-----------------------------------------------------------------------------------------------------------------------------------------------

# Parse arguments
parser = argparse.ArgumentParser(description = "This script identifies FECA-2-LECA duplications, determines the best prokaryotic outgroup (if any), and performs a branch length analysis.")
parser.add_argument("tree", help = 'tree file (either .treefile or .contree)')
parser.add_argument("-p", metavar = "prefix", help = "prefix for output files (DEFAULT: basename tree file), also used for finding the BLAST file")
parser.add_argument("-b", metavar = "BLASTdir", help = 'directory containing the BLAST output files and the file containing all sequence IDs (DEFAULT: current)')
parser.add_argument('-o', metavar = 'outdir', help = 'directory for output files (DEFAULT: current)')
group = parser.add_mutually_exclusive_group()
group.add_argument('-e', help = 'only eukaryotes (DEFAULT: off)', action = 'store_true')
group.add_argument('-i', help = 'filter interspersing prokaryotes (DEFAULT: off)', action = 'store_true')
parser.add_argument('-f', help = 'use only farthest leaf for rooting (DEFAULT: off)', action = 'store_true')
parser.add_argument('-d', metavar = 'xx|0.xx', help = 'threshold for duplication consistency (< 1) or species overlap (> 1) for duplications calling (DEFAULT: 0.2)', default = 0.2, type = float)
parser.add_argument("-l", metavar = "0.xx", help = "coverage threshold for LECA calling (DEFAULT: 0.15)", type = float, default = 0.15)
parser.add_argument('-m', help = 'mode for calculating the branch lengths in case of duplications (minimum (DEFAULT), maximum, median or mean)', choices = ('median', 'minimum', 'maximum', 'mean'), default = 'minimum')
parser.add_argument('-c', help = 'if mode is minimum or maximum, find the branch length corresponding to the minimal/maximal raw branch length (DEFAULT: off)', action = 'store_true')
parser.add_argument("-r", metavar = "root", help = "position of eukaryotic root (DEFAULT: Opimoda-Diphoda)", default = "OD")
parser.add_argument("-s", metavar = "supergroups", help = "supergroups definition used", type = int, choices = (4, 5, 6), default = 5)
args = parser.parse_args()

if not args.p:
    prefix = os.path.basename(args.tree)
    prefix = prefix[:prefix.find('.')]
else:
    prefix = args.p
if args.b:
    blastdir = args.b
else:
    blastdir = '.'
if args.o:
    outdir = args.o
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    outdir_prefix = outdir + '/' + prefix
else:
    outdir_prefix = prefix
euk_only = False; farthest = False; filtering = False; corresponding = False # Default values
if args.e:
    euk_only = True
elif args.i:
    filtering = True
if args.f:
    filtering = True
    sys.stderr.write('Warning: only farthest not fully implemented yet\n')
if args.d < 1:
    duplication_criterion = float(args.d)
    consistency = True
else:
    duplication_criterion = int(args.d)
    consistency = False
coverage_criterion = float(args.l)
if coverage_criterion > 1 or coverage_criterion < 0:
    sys.exit('Error: coverage criterion should be between 0 and 1')
mode = args.m
if args.c:
    if mode in ('minimum', 'maximum'):
        corresponding = True
    else:
        sys.exit('Error: -c can only be used in the minimum or maximum mode')
root_name = args.r
if args.s == 5:
    supergroups = supergroups5
elif args.s == 4:
    supergroups = supergroups4
else:
    supergroups = supergroups6

# Open tree, get supergroups and annotate leaves
tree = open_tree(args.tree)
root_daughters = get_root_daughters(root_name, supergroups4)
annotate_prokaryotic_eukaryotic_leaves(tree, euk_only, root_daughters, supergroups)
leca_dict = {}
dupl_dict = {}
unknown_dict = {}

# Assign all original sequences to representing sequence
with open(f'{blastdir}/{prefix}_seqids.list') as all_seqs_file:
    all_seqs = [line.rstrip() for line in all_seqs_file]
if euk_only:
    tree_seqs = [leaf.name for leaf in tree]
else:
    tree_seqs = [leaf.name for leaf in tree.iter_search_nodes(prok_euk = 'Eukaryote')]
representing = assign_all_seqs(tree_seqs, all_seqs, euk_only, prefix, blastdir)

# Filter interspersing prokaryotes
if filtering:
    if len(tree.search_nodes(prok_euk = 'Prokaryote')) > 1: # At least 2 prokaryotic sequences in the tree
        to_prune = True
        while to_prune:
            to_prune = prok_filter(tree)
        if len(tree.search_nodes(prok_euk = 'Prokaryote')) == 0: # All prokaryotic sequences filtered
            sys.stderr.write('All prokaryotic sequences removed!\n')
            if len(tree.search_nodes(prok_euk = 'Eukaryote')) == 2: # Only 2 eukaryotic sequences remaining
                sys.exit('Error: All prokaryotic sequences removed and only 2 eukaryotic sequences remaining; no further tree analysis possible!')
            euk_only = True
    else:
        if 1935183 in ncbi.get_lineage(tree.search_nodes(prok_euk = 'Prokaryote')[0].taxid): # Asgard archaeon
            sys.stderr.write('Single prokaryotic sequence, but from Asgard archaeon\n')
        else: # Most likely a transfer from eukaryotes to prokaryotes, so removed
            sys.stderr.write('Single prokaryotic sequence removed\n')
            tree.prune(tree.search_nodes(prok_euk = 'Eukaryote'), preserve_branch_length = True)
            euk_only = True

# Perform analysis
if euk_only:
    # First check if the tree fulfills the FECA-2-LECA criterion (always the case for 2 supergroups trees)
    if not feca2leca(tree.get_leaves()):
        sys.exit('No LECAs in this tree')
    # Midpoint rooting to have a root to work with
    root = tree.get_midpoint_outgroup()
    tree.set_outgroup(root)

    # Call duplications in unrooted mode and reroot tree
    tree = annotate_and_reroot_euk_only(tree, supergroups, representing, coverage_criterion, duplication_criterion, consistency)
    if not annotate_overlap_all_assigned(tree, representing, supergroups, coverage_criterion, duplication_criterion, consistency): # Annotate nodes (new duplications and LECAs may be found, because now only one consistency value is considered) and get relevant coverage and redundancy values
        sys.exit(f'No confident LECA in this tree: coverage = {tree.coverage}')
        # Actually, if there are new LECAs now, the position of the root should be recalculated, and then again, the duplications and LECAs should be inferred, and again...

    # Collect information on LECA nodes, duplications (including duplication lengths) and unknown nodes
    collect_leca_information(tree, 1, leca_dict, representing)
    collect_dupl_information(tree, 1, dupl_dict, mode = mode, corresponding = corresponding)
    collect_unknown_information(tree, 1, unknown_dict, representing)

else: # prok + euk
    ancestry_dict = {}
    euk_clades = get_euk_clades(tree)
    feca_clades = []
    non_feca_clades = []
    feca_count = 0
    for euk_clade in euk_clades:
        if feca2leca(euk_clade.get_leaves()):
            feca_clades.append(euk_clade) # To prevent reuse of same euk clades in the loop
        else:
            non_feca_clades.append(euk_clade)
    failed_fecas = []
    for euk_clade in feca_clades:
        feca_confidence = annotate_overlap_all_assigned(euk_clade, representing, supergroups, coverage_criterion, duplication_criterion, consistency)
        if not feca_confidence:
            failed_fecas.append(euk_clade)
            non_feca_clades.append(euk_clade)
    for failed_feca in failed_fecas:
        feca_clades.remove(failed_feca)
    bl_information = {}
    untrusted_fecas = {}
    for euk_clade in feca_clades:
        feca_count += 1
        feca_no = 'FECA' + str(feca_count)
        euk_clade.add_features(feca = True, feca_no = feca_no)
        euk_clade.add_face(TextFace(feca_no), column = 0, position = 'branch-bottom')
    ancestry_out = open(outdir_prefix + '_ancestry.tsv', 'w')
    ancestry_out.write('Pfam\tFECA\tFECA support\tLECAs\tUnknowns\tSister1 name\tSister1 support\tSister2 name\tSister2 support\tLCA sister\tAncestry\tMajority\tSupport\n')
    for feca_count, euk_clade in enumerate(feca_clades): # Again, now FECA nodes assigned
        feca_no = euk_clade.feca_no
        feca_count += 1
        feca_support = euk_clade.support
        tree_copy = tree.copy(method = 'cpickle') # Copy to prevent errors due to rerooting on sequence of other FECA
        feca_clade = tree_copy.search_nodes(feca_no = feca_no)[0]
        lecas = feca_clade.search_nodes(identity = 'LECA')
        unknowns = feca_clade.search_nodes(identity = 'unknown')
        support = None
        ancestry = ''
        majority = False
        if len(tree_copy.search_nodes(prok_euk = 'Prokaryote')) == 1:
            prok = tree_copy.children[0] # As the tree is rooted on this single prokaryote
            ancestry = classify_sister(prok.taxid)
            lca_sel = ncbi.get_taxid_translator([prok.taxid])[prok.taxid]
            lcas = [(prok.taxid, lca_sel, 'NA'),('NA', 'NA', 'NA')]
            support = 'NA'
            feca_support = 'NA'
            sys.stderr.write(f'Warning: one prokaryotic sequence for {prefix}, so branch length analysis may not be accurate and is therefore not performed\n')
        else:
            lcas = get_prokaryotic_sister(feca_clade, tree_copy, farthest)
            ancestries = []
            for lca in lcas:
                pot_ancestry = classify_sister(lca[0])
                ancestries.append(pot_ancestry)
            if ancestries[0] == ancestries[1]:
                ancestry = ancestries[0]
                if 'NA' in (lcas[0][2], lcas[1][2]):
                    support = 'NA'
                else:
                    support = max(lcas[0][2], lcas[1][2])
                farthest_leaf = reroot(feca_clade, tree_copy) # Root on the farthest leaf
                sys.stderr.write(f'Warning: for {feca_no} in {prefix} not enough information to choose a good outgroup and therefore rooted on farthest leaf {farthest_leaf}\n')
                lca_sel = get_prokaryotic_sister(feca_clade, tree_copy, farthest = True)[1]
            else:
                order = ('Alphaproteobacteria', 'Asgard archaea', 'ABG proteobacteria', 'Asgard+TACK group', 'Betaproteobacteria', 'Gammaproteobacteria', 'TACK archaea')
                for ancestry in order:
                    if ancestry in ancestries:
                        index = ancestries.index(ancestry)
                        support = lcas[index][2]
                        lca_sel = lcas[index][1]
                        if index == 1:
                            root = feca_clade.get_sisters()[0] # Root on the old sister
                            tree_copy.set_outgroup(root)
                        break
                else:
                    old_sister_leaves = [leaf.name for leaf in feca_clade.get_sisters()[0]]
                    farthest_leaf = reroot(feca_clade, tree_copy) # Root on the farthest leaf
                    sys.stderr.write(f'Warning: for {feca_no} in {prefix} not enough information to choose a good outgroup and therefore rooted on farthest leaf {farthest_leaf}\n')
                    if farthest_leaf in old_sister_leaves: # To check in which of the two possible sisters the farthest leaf is --> other is the sister group
                        index = 1
                    else:
                        index = 0
                    ancestry = ancestries[index]
                    support = lcas[index][2]
                    lca_sel = lcas[index][1]
            if ancestry in ('Bacteria', 'Archaea', 'cellular organisms'):
                reclassified = reclassify_majority_sister(feca_clade)
                if reclassified:
                    majority = True
                    ancestry = reclassified
        ancestry_dict[feca_no] = ancestry
        print(prefix, feca_no, feca_support, len(lecas), len(unknowns), lcas[0][1], lcas[0][2], lcas[1][1], lcas[1][2], lca_sel, ancestry, majority, support, sep = '\t', file = ancestry_out)

        # Check if there are FECAs in the sistergroup and if not, calculate branch lengths
        sister = feca_clade.get_sisters()[0]
        sister_fecas = sister.search_nodes(feca = True)
        if len(sister_fecas) > 0:
            sys.stderr.write(f'Warning: for {feca_no} in {prefix} there is another FECA in the sister group and therefore the branch lengths may not be accurate\n')
            untrusted_fecas[feca_no] = '+'.join([sister_feca.feca_no for sister_feca in sister_fecas])
            bl_information[feca_no] = [ancestry, str(len(lecas))] + ['NA'] * 4
        else:
            lepca = feca_clade.up
            psbl, rsl, sl, ebl = calculate_stem_lengths(tree_copy, lecas, lepca, sister, mode = mode, corresponding = corresponding)
            bl_information[feca_no] = [ancestry, str(len(lecas)), str(psbl), str(rsl), str(sl), str(ebl)]

        # Collect information on LECA nodes, duplications (including duplication lengths) and unknown nodes
        collect_leca_information(euk_clade, feca_count, leca_dict, representing)
        collect_dupl_information(euk_clade, feca_count, dupl_dict, mode = mode, corresponding = corresponding)
        collect_unknown_information(euk_clade, feca_count, unknown_dict, representing)

    ancestry_out.close()
    # Annotate non-FECA clades
    with open(outdir_prefix + '_non_feca.tsv', 'w') as non_feca_out:
        non_feca_out.write('Pfam\tNon-FECA\tDonor taxon\tAncestry\tSupport\tCoverage\tSpecies\tSequence IDs\tRepresenting species\n')
        for count,euk_clade in enumerate(non_feca_clades):
            count += 1
            euk_clade.add_features(non_feca_no = count)
            euk_clade.add_face(TextFace('Non-FECA' + str(count), fgcolor = 'red'), column = 0, position = 'branch-top')
            tree_copy = tree.copy(method = 'cpickle')
            non_feca_clade = tree_copy.search_nodes(non_feca_no = count)[0]
            lca, lca_name, support = get_non_feca_sister(non_feca_clade, tree_copy)
            ancestry = classify_sister(lca)
            species = []
            sequences = []
            repr_species = []
            for leaf in non_feca_clade:
                species.append(leaf.taxid)
                sequences.append(leaf.name)
            coverage, copy_no, list_species = infer_coverage_redundancy(euk_clade, representing, supergroups)
            print(prefix, 'Non-FECA' + str(count), lca_name, ancestry, support, coverage, ','.join(species), ','.join(sequences), ','.join(list_species), sep = '\t', file = non_feca_out)

    # Print branch lengths and FECA information
    if feca_count != 0:
        with open(outdir_prefix + '_fecas.tsv', 'w') as fecas_out, open(outdir_prefix + '_branch_lengths.tsv', 'w') as branch_lengths_out:
            fecas_out.write('Pfam\tFECA\tTrusted\tAncestry\tLECAs\tSister-FECAs\n')
            branch_lengths_out.write('Pfam\tFECA\tAncestry\tLECAs\tProkaryotic sister branch length\tRaw stem lengths\tStem lengths\tEukaryotic branch lengths\n')
            if len(tree.search_nodes(prok_euk = 'Prokaryote')) == 1: # No accurate branch lengths
                print(prefix, 'FECA1', ancestry, len(tree.search_nodes(identity = 'LECA')), '\t'.join(['NA']*4), sep = '\t', file = branch_lengths_out)
                print(prefix, 'FECA1', 'Yes', ancestry, len(tree.search_nodes(identity = 'LECA')), 'NA', sep = '\t', file = fecas_out)
            else:
                for feca_no in bl_information:
                    print(prefix, feca_no, '\t'.join(bl_information[feca_no]), sep = '\t', file = branch_lengths_out)
                    if feca_no in untrusted_fecas:
                        print(prefix, feca_no, 'No', '\t'.join(bl_information[feca_no][0:2]), untrusted_fecas[feca_no], sep = '\t', file = fecas_out)
                    else:
                        print(prefix, feca_no, 'Yes', '\t'.join(bl_information[feca_no][0:2]), 'NA', sep = '\t', file = fecas_out)

    # Merge untrusted FECAs with (un)trusted FECAs in case of aunt/niece (non-FECAs in between are ignored and if merged included in the merged FECA)
    merged_fecas = []
    strict_fecas = {}
    # Check if tree not rooted within this FECA?
    for feca in untrusted_fecas:
        if feca in merged_fecas:
            continue
        sister_feca = untrusted_fecas[feca] # Get the/a FECA in the sistergroup to root the tree
        if '+' in sister_feca:
            sister_feca = sister_feca[:sister_feca.find('+')]
        feca_node = tree.search_nodes(feca_no = feca)[0]
        sister_feca_node = tree.search_nodes(feca_no = sister_feca)[0]
        sister_node = feca_node.get_sisters()[0] # What if not bifurcating?!
        if sister_feca_node not in sister_node: # First reroot
            tree.set_outgroup(sister_node)
            sister_node = feca_node.get_sisters()[0] # What if not bifurcating?!
        prok = False
        merged = False
        mfeca_count = len(tree.search_nodes(identity = 'mFECA')) + 1
        while(not prok):
            try:
                if sister_node.identity == 'mFECA':
                    mfeca = sister_node.mfeca_no
                    try:
                        strict_fecas[mfeca_count] += strict_fecas[mfeca]
                    except KeyError:
                        strict_fecas[mfeca_count] = [feca] + strict_fecas[mfeca]
                    del strict_fecas[mfeca]
                    merged = True
                    merged_fecas.append(feca)
                    break
            except AttributeError:
                pass
            for niece in sister_node.get_children():
                if niece in feca_clades:
                    incl_feca = niece.feca_no
                    merged_fecas.append(incl_feca)
                    try:
                        strict_fecas[mfeca_count].append(incl_feca)
                    except KeyError:
                        strict_fecas[mfeca_count] = [feca, incl_feca]
                        merged_fecas.append(feca)
                    sister_node = niece.get_sisters()[0] # What if not bifurcating?!
                    merged = True
                    break
                elif niece in non_feca_clades:
                    incl_non_feca = 'Non-FECA' + str(niece.non_feca_no)
                    try:
                        strict_fecas[mfeca_count].append(incl_non_feca)
                    except KeyError:
                        strict_fecas[mfeca_count] = [feca, incl_non_feca]
                        merged_fecas.append(feca)
                    sister_node = niece.get_sisters()[0]
                    merged = True
                    break
            else:
                prok = True
        if merged:
            feca_node.up.add_features(identity = 'mFECA', mfeca_no = mfeca_count)
            feca_node.up.add_face(TextFace('mFECA' + str(mfeca_count)), column = 0, position = 'branch-top')
    with open(outdir_prefix + '_fecas_strict.tsv', 'w') as fecas_strict_out:
        fecas_strict_out.write('Pfam\tFECA\tTrusted\tAncestry\tLECAs\tContent\n')
        for feca_no in bl_information:
            if feca_no in merged_fecas:
                continue
            if feca_no in untrusted_fecas:
                trusted = 'No'
            else:
                trusted = 'Yes'
            print(prefix, feca_no, trusted, '\t'.join(bl_information[feca_no][0:2]), 'NA', sep = '\t', file = fecas_strict_out)
        for mfeca,content in strict_fecas.items():
            feca_no = 'FECA' + '_'.join([feca[4:] for feca in content if 'Non' not in feca])
            ancestry = bl_information[content[0]][0] # Should be all the same
            lecas = 0
            for feca in content:
                if 'Non' not in feca:
                    lecas += int(bl_information[feca][1])
            print(prefix, feca_no, 'No', ancestry, lecas, '+'.join(content), sep = '\t', file = fecas_strict_out)

no_lecas = len(tree.search_nodes(identity = 'LECA'))
no_unknowns = len(tree.search_nodes(identity = 'unknown'))
print('Pfam\tFECAs (normal)\tFECAs (strict)\tFECAs (after merging)\tLECAs\tUnknowns\tNon-FECAs')
if euk_only:
    print(prefix, 'NA', 'NA', 'NA', no_lecas, no_unknowns, 'NA', sep = '\t')
else:
    print(prefix, feca_count, feca_count - len(untrusted_fecas), feca_count - len(set(merged_fecas)) + len(strict_fecas), no_lecas, no_unknowns, len(non_feca_clades), sep = '\t')

# Write output of LECA, duplication and unknown nodes
with open(outdir_prefix + '_lecas.tsv', 'w') as lecas_out, open(outdir_prefix + '_lecas_all_seqs.tsv', 'w') as lecas_all_seqs_out:
    lecas_all_seqs_out.write('Pfam\tLECA\tSupport\tTree seqs\tRepresenting seqs\n')
    lecas_out.write('Pfam\tFECA\tAncestry\tLECA\tSupport\tCoverage\tCopy number\tSeqs\n')
    for leca, info in leca_dict.items():
        if euk_only:
            feca_no = 'NA'
            ancestry = 'Eukaryotic'
        else:
            feca_no = 'FECA' + leca[2:leca.find('.')]
            ancestry = ancestry_dict[feca_no]
        print(prefix, feca_no, ancestry, leca, '\t'.join(info[:4]), sep = '\t', file = lecas_out)
        print(prefix, leca, info[0], info[3], info[4], sep = '\t', file = lecas_all_seqs_out)

with open(outdir_prefix + '_duplication_lengths.tsv', 'w') as duplication_lengths_out:
    duplication_lengths_out.write('Pfam\tFECA\tAncestry\tDuplication\tSupport\tSpecies overlap\tDuplication consistency\tOGs\tRaw duplication lengths\tDuplication lengths\tEukaryotic branch lengths\n')
    for dupl, info in dupl_dict.items():
        if euk_only:
            feca_no = 'NA'
            ancestry = 'Eukaryotic'
        else:
            feca_no = 'FECA' + dupl[1:dupl.find('.')]
            ancestry = ancestry_dict[feca_no]
        print(prefix, feca_no, ancestry, dupl, '\t'.join([str(x) for x in info]), sep = '\t', file = duplication_lengths_out)

with open(outdir_prefix + '_unknowns.tsv', 'w') as unknowns_out:
    unknowns_out.write('Pfam\tFECA\tAncestry\tUnknowns\tSupport\tCoverage\tSeqs\n')
    for unknown, info in unknown_dict.items():
        if euk_only:
            feca_no = 'NA'
            ancestry = 'Eukaryotic'
        else:
            feca_no = 'FECA' + unknown[1:unknown.find('.')]
            ancestry = ancestry_dict[feca_no]
        print(prefix, feca_no, ancestry, unknown, '\t'.join(info), sep = '\t', file = unknowns_out)

# Output annotated tree
ts.show_branch_support = True
ts.show_leaf_name = False
tree.render(outdir_prefix + '_annotated_tree.pdf', tree_style = ts)
tree.write(format = 1, outfile = outdir_prefix + '_annotated_tree.nw')
