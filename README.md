# Analysis of duplications during eukaryogenesis (from FECA to LECA)

This repository contains scripts to annotate phylogenetic trees containing duplications during eukaryogenesis, as used in "Timing the origin of eukaryotic cellular complexity with ancient duplications" (Vosseberg, Van Hooff, et al.; preprint at bioRxiv: https://www.biorxiv.org/content/10.1101/823484v1).

## Dependencies
- Python 3.6 or higher
- ETE3
- NumPy

## Pfam-ScrollSaw trees

### Input
- Newick tree file
- BLAST files: an all-versus-all (`<pfam>\_blastp.txt`) and own-versus-own file (`<pfam>\_blastp\_own.txt`)
- File containing a list of all sequence IDs prior to ScrollSaw downsampling (`<pfam>_seqids.list`)

### Usage
```
usage: annotate_scrollsaw_tree.py [-h] [-p prefix] [-b BLASTdir] [-o outdir]
                                  [-e | -i] [-f] [-d xx|0.xx] [-l 0.xx]
                                  [-m {median,minimum,maximum,mean}] [-c]
                                  [-r root] [-s supergroups]
                                  tree

This script identifies FECA-2-LECA duplications, determines the best
prokaryotic outgroup (if any), and performs a branch length analysis.

positional arguments:
  tree                  tree file (either .treefile or .contree)

optional arguments:
  -h, --help            show this help message and exit
  -p prefix             prefix for output files (DEFAULT: basename tree file),
                        also used for finding the BLAST file
  -b BLASTdir           directory containing the BLAST output files and the
                        file containing all sequence IDs (DEFAULT: current)
  -o outdir             directory for output files (DEFAULT: current)
  -e                    only eukaryotes (DEFAULT: off)
  -i                    filter interspersing prokaryotes (DEFAULT: off)
  -f                    use only farthest leaf for rooting (DEFAULT: off)
  -d xx|0.xx            threshold for duplication consistency (< 1) or species
                        overlap (> 1) for duplications calling (DEFAULT: 0.2)
  -l 0.xx               coverage threshold for LECA calling (DEFAULT: 0.15)
  -m {median,minimum,maximum,mean}
                        mode for calculating the branch lengths in case of
                        duplications (minimum (DEFAULT), maximum, median or
                        mean)
  -c                    if mode is minimum or maximum, find the branch length
                        corresponding to the minimal/maximal raw branch length
                        (DEFAULT: off)
  -r root               position of eukaryotic root (DEFAULT: Opimoda-Diphoda)
  -s supergroups        supergroups definition used

```
_Typical usage_
```
annotate_scrollsaw_tree.py -e -b examples -o test examples/PF03028.iqtree.LG4X.treefile
annotate_scrollsaw_tree.py -i -b examples -o test examples/PF04157.iqtree.LG4x.contree
```

## KOG-to-COG trees

### Input
- Newick tree file

### Usage
```
usage: annotate_kogcog_tree.py [-h] [-p prefix] [-o outdir] [-e | -i] [-f]
                               [-d xx|0.xx] [-m {median,minimum,maximum,mean}]
                               [-c] [-r root] [-s supergroups]
                               tree

This script identifies FECA-2-LECA duplications in KOG-COG trees, determines
the best prokaryotic outgroup (if any), and performs a branch length analysis.

positional arguments:
  tree                  tree file (either .treefile or .contree)

optional arguments:
  -h, --help            show this help message and exit
  -p prefix             prefix for output files (DEFAULT: basename tree file)
  -o outdir             directory for output files (DEFAULT: current)
  -e                    only eukaryotes (DEFAULT: off)
  -i                    filter interspersing prokaryotes (DEFAULT: off)
  -f                    use only farthest leaf for rooting (DEFAULT: off)
  -d xx|0.xx            threshold for duplication consistency (< 1) or species
                        overlap (> 1) for duplications calling (DEFAULT: 0.2)
  -m {median,minimum,maximum,mean}
                        mode for calculating the branch lengths in case of
                        duplications (minimum (DEFAULT), maximum, median or
                        mean)
  -c                    if mode is minimum or maximum, find the branch length
                        corresponding to the minimal/maximal raw branch length
                        (DEFAULT: off)
  -r root               position of eukaryotic root (DEFAULT: Opimoda-Diphoda)
  -s supergroups        supergroups definition used
```

_Typical usage_
```
annotate_kogcog_tree.py -e -o test examples/Efamily_100.tree.fasttree.ml.wag.contree
annotate_kogcog_tree.py -i -o test examples/COG0002.tree.fasttree.ml.wag.contree
```

## Output files
- ancestry.tsv (if tree includes prokaryotes): per acquisition the two candidate sister groups and which one was chosen
- fecas.tsv (if tree includes prokaryotes): per acquisition if there was another acquisition in its sister clade
- fecas_strict.tsv (if tree includes prokaryotes): acquisitions after merging nested acquisitions
- non_feca.tsv (if tree includes prokaryotes): eukaryotic clades not fulfilling the acquisition criteria
- branch_lengths.tsv (if tree includes prokaryotes): stem lengths for each acquisition
- duplication_lengths.tsv: information about each duplication including its duplication length
- lecas.tsv: information about each LECA family
- lecas_all_seqs.tsv (only ScrollSaw): per LECA family all sequences it is comprised of
- unknowns.tsv: eukaryotic clades in an acquisition not fulfilling the LECA criteria
- annotated_tree.nw: annotated Newick tree file
- annotated_tree.pdf: PDF with the annotated tree
