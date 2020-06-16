#!/usr/bin/env python

# Functions for annotating FECA-to-LECA trees

import sys
from ete3 import PhyloTree
from ete3 import NCBITaxa
from ete3 import TextFace
import collections
from numpy import mean, median
from eukarya import *
from scrollsaw import *
ncbi = NCBITaxa()

def open_tree(tree_file_path):
    """Opens tree (contree or treefile) and assigns support values to nodes in case of a standard tree file"""
    if 'contree' in tree_file_path:
        tree = PhyloTree(tree_file_path, sp_naming_function = None)
    elif 'treefile' in tree_file_path: # Branch supports in SH-aLRT support (%) / ultrafast bootstrap support (%)
        tree = PhyloTree(tree_file_path, sp_naming_function = None, format = 1)
        for node in tree.iter_descendants():
            if not node.is_leaf():
                support_values = node.name.split('/')
                try:
                    node.support = float(support_values[1])
                except IndexError: # No support values when sequences were identical --> set support artifically to 100.0
                    node.support = 100.0
                #node.add_features(shalrt = float(support_values[0])) # Not necessary...
    else:
        sys.exit('Error: tree format not recognised')
    return tree

def annotate_prokaryotic_eukaryotic_leaves(tree, euk_only, root_daughters = supergroups2, supergroups = supergroups5):
    """Distinguishes prokaryotic (NCBI taxid) and eukaryotic leave names and annotates them"""
    if euk_only:
        for leaf in tree:
            taxid = leaf.name[0:4]
            try:
                leaf.add_features(taxid = taxid, root_daughter = root_daughters[taxid], supergroup = supergroups[taxid])
            except KeyError:
                sys.exit(f'Error: species {taxid} (sequence ID: {leaf.name}) not recognised')
    else:
        for leaf in tree:
            if leaf.name[0].isdigit():
                taxid = int(leaf.name[:leaf.name.find('.')])
                leaf.add_features(taxid = taxid, prok_euk = 'Prokaryote')
            else:
                taxid = leaf.name[0:4]
                leaf.add_features(taxid = taxid, root_daughter = root_daughters[taxid], supergroup = supergroups[taxid], prok_euk = 'Eukaryote')

def prok_filter(tree, level = 'genus'):
    """Filters interspersing prokaryotes (either single or same species or other level).
Note: if there is only one prokaryotic leaf, this one will be removed."""
    euk_leave_names = [leaf.name for leaf in tree.iter_search_nodes(prok_euk = 'Eukaryote')]
    tree.set_outgroup(euk_leave_names[0])
    filtered_leaves = []
    for clade in tree.get_monophyletic(values = ['Prokaryote'], target_attr = 'prok_euk'):
        if level == 'single' and len(clade) != 1: # So, only singletons
            continue
        elif len(clade) != 1: # Check if same level if not singleton
            previous = ''
            try:
                if clade.diverse:
                    continue
            except AttributeError:
                diverse = False
                for leaf in clade:
                    ranks = ncbi.get_rank(ncbi.get_lineage(leaf.taxid))
                    for name, rank in ranks.items():
                        if rank == level:
                            if previous == '':
                                previous = name
                            elif previous != name:
                                diverse = True
                            break
                    else: # Level not detected
                        sys.stderr.write(f'Warning: {level} level not detected for {leaf.name}\n')
                        diverse = True
                    if diverse:
                        break
                if diverse:
                    clade.add_features(diverse = True)
                    continue
        parent = clade.up
        if len(clade) == 1:
            prok_name = clade.name
            prok_leaves = [clade]
        else:
            prok_leaves = [leaf for leaf in clade]
            prok_name = ','.join([leaf.name for leaf in prok_leaves])
            clade.name = prok_name
        sister = clade.get_sisters()
        if len(sister) == 1:
            sister = sister[0]
        else:
            continue
        if len(parent.get_sisters()) != 1:
            continue
        if parent.up.is_root(): # Immediate sister is where the tree is rooted on
            for child in sister.get_children():
                if len(child.search_nodes(prok_euk = 'Prokaryote')) == 0:
                    sys.stderr.write(f'Interspersing prokaryote type 1: {prok_name}\n')
                    filtered_leaves.extend(prok_leaves)
                    break
        elif len(sister.search_nodes(prok_euk = 'Prokaryote')) != 0:
            tree.set_outgroup(parent)
            clade = tree&prok_name
            parent = clade.up
            sister = clade.get_sisters()[0]
            new_sister = parent.get_sisters()[0]
            if len(new_sister.search_nodes(prok_euk = 'Prokaryote')) == 0:
                for child in sister.get_children():
                    if len(child.search_nodes(prok_euk = 'Prokaryote')) == 0:
                        sys.stderr.write(f'Interspersing prokaryote type 2: {prok_name}\n')
                        filtered_leaves.extend(prok_leaves)
                        tree.set_outgroup(tree&euk_leave_names[0])
                        break
                else:
                    tree.set_outgroup(tree&euk_leave_names[0])
            else:
                tree.set_outgroup(tree&euk_leave_names[0])
        else:
            p_sister = parent.get_sisters()[0]
            if len(p_sister.search_nodes(prok_euk = 'Prokaryote')) == 0:
                sys.stderr.write(f'Interspersing prokaryote type 3: {prok_name}\n')
                filtered_leaves.extend(prok_leaves)
            else:
                tree.set_outgroup(p_sister)
                clade = tree&prok_name
                parent = clade.up
                new_sister = parent.get_sisters()[0]
                if len(new_sister.search_nodes(prok_euk = 'Prokaryote')) == 0:
                    sys.stderr.write(f'Interspersing prokaryote type 4: {prok_name}\n')
                    filtered_leaves.extend(prok_leaves)
                tree.set_outgroup(tree&euk_leave_names[0])
    if len(filtered_leaves) == 0:
        pruned = False
        return(pruned)
    leaves_kept = [leaf for leaf in tree if leaf not in filtered_leaves]
    tree.prune(leaves_kept, preserve_branch_length = True)
    pruned = True
    return pruned

def feca2leca(leaves):
    """Determines if clade fulfills FECA-2-LECA criteria: both sides of the root present"""
    root_daughter_groups = set()
    for leaf in leaves:
        root_daughter_groups.add(leaf.root_daughter)
        if len(root_daughter_groups) == 2:
            break
    else:
        return False
    return True

def duplication_check(leaves1, leaves2):
    """Determines if a node fulfills the FECA-2-LECA duplication criteria: both daughters have both Opimoda and Diphoda"""
    if feca2leca(leaves1) and feca2leca(leaves2):
        return True
    else:
        return False

def duplication_check_unrooted(node, tree, tips, supergroups, representing = None, coverage_criterion = 0.15, duplication_criterion = 0.2, consistency = True, scrollsaw = True):
    """Duplication check in euk only tree, see annotate_and_reroot_euk_only"""
    if len(node.get_children()) > 2:
        node.resolve_polytomy(recursive = False)
    leaves1, leaves2 = [child.get_leaves() for child in node.get_children()]
    leaves3 = tips - (set(leaves1) | set(leaves2))
    if not feca2leca(leaves1) or not feca2leca(leaves2) or not feca2leca(leaves3):
        return False
    species = []
    for leaves in [leaves1, leaves2, leaves3]:
        if scrollsaw:
            coverage, copies, repr_species = infer_coverage_redundancy(leaves, representing, supergroups)
            if coverage < coverage_criterion:
                return False
            species.append(set(repr_species))
        else:
            species.append(set([leaf.taxid for leaf in leaves]))
    overlaps = [len(species[0] & species[1]), len(species[0] & species[2]), len(species[1] & species[2])]
    if consistency:
        dupl_cons = [overlaps[0] / len(species[0] | species[1]), overlaps[1] / len(species[0] | species[2]), overlaps[2] / len(species[1] | species[2])]
        if dupl_cons[0] >= duplication_criterion and dupl_cons[1] >= duplication_criterion and dupl_cons[2] >= duplication_criterion:
            return True
        else:
            return False
    else:
        if overlaps[0] >= duplication_criterion and overlaps[1] >= duplication_criterion and overlaps[2] >= duplication_criterion:
            return True
        else:
            return False

def duplication_check_rooted(node, supergroups, representing = None, coverage_criterion = 0.15, duplication_criterion = 0.2, consistency = True, scrollsaw = True):
    """Duplication check in euk only tree, see annotate_and_reroot_euk_only"""
    daughter1, daughter2 = node.get_children()
    if not duplication_check(daughter1.get_leaves(), daughter2.get_leaves()):
        return False
    if scrollsaw:
        repr_daughter = {}
        repr_daughter[1] = infer_coverage_redundancy(daughter1, representing, supergroups)
        if repr_daughter[1][0] < coverage_criterion:
            return False
        repr_daughter[2] = infer_coverage_redundancy(daughter2, representing, supergroups)
        if repr_daughter[2][0] < coverage_criterion:
            return False
        overlap = len(set(repr_daughter[1][2]) & set(repr_daughter[2][2]))
        dupl_consistency = overlap / len(set(repr_daughter[1][2]) | set(repr_daughter[2][2]))
    else:
        species1 = set([leaf.taxid for leaf in daughter1.get_leaves()])
        species2 = set([leaf.taxid for leaf in daughter2.get_leaves()])
        overlap = len(species1 & species2)
        dupl_consistency = overlap / len(species1 | species2)
    if consistency and dupl_consistency < duplication_criterion or not consistency and overlap < duplication_criterion:
        return False
    else:
        if consistency:
            return dupl_consistency
        else:
            return overlap

def annotate_and_reroot_euk_only(tree, supergroups, representing = None, coverage_criterion = 0.15, duplication_criterion = 0.2, consistency = True, scrollsaw = True):
    """Root tree on mid of longest possible distance between the LECAs. Written by Jolien, adapted for ScrollSaw trees."""
    for node in tree.traverse("preorder"):
        node.add_features(identity = "?")
    tips = set(tree.get_leaves())
    dup_counter=0
    for node in tree.iter_descendants("preorder"): ##don't visit root: only two directions (only descendants)
        if not node.is_leaf():
            if duplication_check_unrooted(node, tree, tips, supergroups, representing, coverage_criterion, duplication_criterion, consistency, scrollsaw):
                # Duplications only called if there are at least 2 duplications, otherwise no internal node fulfilling this check
                dup_counter += 1
                node.add_features(identity = "duplication", name = "D"+str(dup_counter))
    if dup_counter > 0: # So, at least 2 duplications in the tree
        duplications = tree.search_nodes(identity = 'duplication')
        duplication_paths = {}
        for dup in duplications:
            if dup.get_children()[0] in duplications or dup.get_children()[1] in duplications:
                continue
            name = dup.name
            duplication_paths[name] = [dup.name]
            tree_copy = tree.copy(method = "cpickle")
            search_node = tree_copy&name
            while search_node:
                parent = search_node.up
                if parent.identity == "duplication":
                    duplication_paths[name].append(parent.name)
                    search_node = parent
                elif parent.is_root():
                    duplication_paths[name].append("root")
                    del tree_copy
                    break
                else:
                    del tree_copy
                    break
        convergence_point = None
        # Check if there is a continuous path
        for i, dupl in enumerate(duplication_paths):
            if i == 0:
                convergence_point = set(duplication_paths[dupl])
            else:
                convergence_point = convergence_point & set(duplication_paths[dupl])
        if len(convergence_point) == 0:
            sys.stderr.write(f"Tree contains more than one duplication path with {dup_counter} duplications in total\n")
        duplications = tree.search_nodes(identity = "duplication")
        tree.set_outgroup(duplications[0]) ## tmp rooting on duplication node
        lecas = []
        leca_counter = 0
        for d in duplications:
            daughters = d.get_children()
            if d.up.is_root():
                daughters.append(d.get_sisters()[0])
            else:
                daughters.append(d.up)
            for daughter in daughters:
                if daughter not in duplications and len(daughter.search_nodes(identity = 'duplication')) == 0:
                    leca_counter += 1
                    name = "OG" + str(leca_counter)
                    daughter.add_features(identity = "LECA", name = name)
                    lecas.append(name)
        ##calculate all distances between LECAs, find the longest distance, root in the middle of this distance
        ##root = duplication
        longest_distance = 0
        longest_distance_node1 = ""
        longest_distance_node2 = ""
        for x in range(0, len(lecas)-1):
            for y in range(1, len(lecas)):
                if x < y:
                    distance = tree.get_distance(lecas[x], lecas[y])
                    if distance > longest_distance:
                        longest_distance = distance
                        longest_distance_node1 = lecas[x]
                        longest_distance_node2 = lecas[y]
        mid_longest_distance = 0.5 * longest_distance
        ##find path between node1 and node2 via their last common ancestor
        node1 = tree&longest_distance_node1
        node2 = tree&longest_distance_node2
        lca = tree.get_common_ancestor(node1, node2)
        path = collections.OrderedDict() # Ordered dictionary in which the distance between nodes are stored
        while node1:
            if node1 != lca:
                path[node1] = tree.get_distance(node1, node1.up)
                node1 = node1.up
            else:
                break
        path2 = collections.OrderedDict()
        while node2:
            if node2 != lca:
                path2[node2.up] = tree.get_distance(node2, node2.up)
                node2 = node2.up
            else:
                break
        path2_reverse = collections.OrderedDict(reversed(path2.items()))
        path.update(path2_reverse)
        ##make a list of the path
        node_length_path = []
        for key, value in path.items():
            node_length_path.extend([key, value])
        node_length_path.append(tree&longest_distance_node2) ##complete the path list
        ##find root position
        nodes_on_path = node_length_path[0::2] ##take nodes only
        parent_node_root = ""
        child_node_root = ""
        parent_node_distance = 0
        child_node_distance = 0
        dist_sum = 0
        for x in range(len(nodes_on_path)):
            n = nodes_on_path[x]
            dist = node_length_path[x * 2 + 1]
            dist_sum += dist
            if dist_sum > mid_longest_distance:
                n1 = nodes_on_path[x + 1]
                if n1 == n.up:
                    child_node_root = n
                    parent_node_root = n1
                    parent_node_distance = dist_sum - mid_longest_distance
                else:
                    child_node_root = n1
                    parent_node_root = n
                    parent_node_distance = mid_longest_distance - (dist_sum - dist)
                child_node_distance = dist - parent_node_distance
                break
        ##set root position
        child_node_root_detached = child_node_root.detach()
        child_node_root_detached.add_features(dist = child_node_distance)
        parent_node_root.add_child(name = "R", dist = parent_node_distance)
        R = tree&"R"
        R.add_child(name = "O", dist = 0) ##outgroup
        R.add_child(child_node_root_detached)
        tree.set_outgroup("O")
        tree = R.detach()
        R.add_features(identity = "duplication", name = "D" + str(dup_counter + 1), support = 101)
        supports = [child.support for child in R.get_children()]
        if supports[0] != supports[1]:
            if supports[0] == 1.0:
                R.get_children()[0].support = supports[1]
            elif supports[1] == 1.0:
                R.get_children()[1].support = supports[0]
            else:
                sys.exit('Error with duplicating support values at both sides of the root')
    else:
        ## Either no duplications or 1 duplication
        ## Try rooting on each internal node
        ## Would the root be a duplication node?
        ## What is the number of species overlap or duplication consistency?
        ## Define LECAs for maximal species overlap or duplication consistency
        tmp_outgroup = ""
        tmp_overlap = 0
        tree_copy = tree.copy(method = 'cpickle') ##keeps original rooting position and lengths
        default_outgroup = tree.get_children()[0]
        tree.set_outgroup(default_outgroup)
        for node in tree.iter_descendants("preorder"):
            if node.is_leaf():
                continue
            tree.set_outgroup(node)
            check = duplication_check_rooted(tree, supergroups, representing, coverage_criterion, duplication_criterion, consistency, scrollsaw) # If fulfilling criteria length overlap, else True/False
            if check:
                if check > tmp_overlap:
                    tmp_outgroup = node
                    tmp_overlap = check
            ##reset the root
            ##branch lengths to root likely will be altered, but that doesn't matter for now.
            tree.set_outgroup(default_outgroup)
        if tmp_outgroup != "":
            outgroup = tmp_outgroup
            tree.set_outgroup(outgroup)
            sister_outgroup = outgroup.get_sisters()[0]
            outgroup.add_features(identity = "LECA", name = "OG1")
            sister_outgroup.add_features(identity = "LECA", name = "OG2")
            root = tree.get_tree_root()
            root.add_features(identity = "duplication", name = "D1")
        ##no duplications found: root = LECA
        else:
            tree = tree_copy ##take back the original tree
            tree.add_features(identity = "LECA", name = 'OG1') # For ebl would be nice to have a good root here as well
    return tree

def get_euk_clades(tree):
    """Get monophyletic eukaryotic clades"""
    root = tree.search_nodes(prok_euk = 'Prokaryote')[0] # Root on first prokaryotic sequence
    tree.set_outgroup(root)
    clades = tree.get_monophyletic(values = ['Eukaryote'], target_attr = 'prok_euk')
    return clades

def annotate_overlap_all_assigned(feca, representing, supergroups, coverage_criterion = 0.15, duplication_criterion = 0.2, consistency = True):
    """Annotate eukaryotic nodes, including the represented sequences"""
    for node in feca.traverse('preorder'):
        if not node.is_leaf():
            if len(node.get_children()) > 2:
                node.resolve_polytomy(recursive = False)
        coverage, copies, repr_species = infer_coverage_redundancy(node, representing, supergroups)
        node.add_features(coverage = coverage, redundancy = copies, repr_species = repr_species, identity = '?')
    for node in feca.traverse('preorder'):
        if node.is_leaf():
            continue
        leaves = node.get_leaves()
        if not feca2leca(leaves):
            continue
        daughter1, daughter2 = node.get_children()
        species_overlap = set(daughter1.repr_species) & set(daughter2.repr_species)
        all_species = set(daughter1.repr_species) | set(daughter2.repr_species)
        dupl_consistency = len(species_overlap) / len(all_species)
        if consistency and dupl_consistency >= duplication_criterion or not consistency and len(species_overlap) >= duplication_criterion:
            if daughter1.coverage >= coverage_criterion and daughter2.coverage >= coverage_criterion:
                if duplication_check(daughter1.get_leaves(), daughter2.get_leaves()):
                    node.add_features(identity = 'duplication', overlap = len(species_overlap), consistency = dupl_consistency)
                    continue
        if node == feca:
            if node.coverage >= coverage_criterion:
                node.add_features(identity = 'LECA')
            else:
                node.add_features(coverage = coverage)
                return False
        elif node.up.identity == 'duplication':
            node.add_features(identity = 'LECA')
        else:
            node.add_features(identity = 'post-LECA')
    for leca in feca.iter_search_nodes(identity = 'LECA'):
        if len(leca.search_nodes(identity = 'duplication')) > 0:
            leca.identity = 'unknown'
        else:
            lecas_down = leca.search_nodes(identity = 'LECA')
            if len(lecas_down) > 1:
                for leca_down in lecas_down[1:]:
                    leca_down.identity = '?'
    for postleca in feca.iter_search_nodes(identity = 'post-LECA'):
        if len(postleca.search_nodes(identity = 'LECA')) > 0:
            postleca.identity = 'unknown'
    # Change rare unknown nodes that don't fulfill the duplication criterion but are duplications because there are duplications in both their children
    for unknown in feca.iter_search_nodes(identity = "unknown"):
        daughter1, daughter2 = unknown.get_children()
        if len(daughter1.search_nodes(identity = 'duplication')) > 0 and len(daughter2.search_nodes(identity = 'duplication')) > 0:
            unknown.identity = 'duplication'
            species_overlap = set(daughter1.repr_species) & set(daughter2.repr_species)
            all_species = set(daughter1.repr_species) | set(daughter2.repr_species)
            dupl_consistency = len(species_overlap) / len(all_species)
            unknown.add_features(overlap = len(species_overlap), consistency = dupl_consistency)
    return True

def annotate_non_scrollsaw(feca, duplication_criterion = 0.2, consistency = True):
    """Annotate eukaryotic nodes for non-ScrollSaw trees"""
    for node in feca.traverse('preorder'):
        if not node.is_leaf():
            if len(node.get_children()) > 2:
                node.resolve_polytomy(recursive = False)
        node.add_features(identity = '?')
    for node in feca.traverse('preorder'):
        if node.is_leaf():
            continue
        leaves = node.get_leaves()
        if not feca2leca(leaves):
            continue
        daughters = node.get_children()
        if duplication_check(daughters[0].get_leaves(), daughters[1].get_leaves()):
            daughter1_sp = set([leaf.taxid for leaf in daughters[0]])
            daughter2_sp = set([leaf.taxid for leaf in daughters[1]])
            species_overlap = daughter1_sp & daughter2_sp
            all_species = daughter1_sp | daughter2_sp
            dupl_consistency = len(species_overlap) / len(all_species)
            if consistency and dupl_consistency >= duplication_criterion or not consistency and len(species_overlap) >= duplication_criterion:
                node.add_features(identity = 'duplication', overlap = len(species_overlap), consistency = dupl_consistency)
                continue
        if node == feca:
            node.add_features(identity = 'LECA')
        elif node.up.identity == 'duplication':
            node.add_features(identity = 'LECA')
        else:
            node.add_features(identity = 'post-LECA')
    for leca in feca.iter_search_nodes(identity = 'LECA'):
        if len(leca.search_nodes(identity = 'duplication')) > 0:
            leca.identity = 'unknown'
        else:
            lecas_down = leca.search_nodes(identity = 'LECA')
            if len(lecas_down) > 1:
                for leca_down in lecas_down[1:]:
                    leca_down.identity = '?'
    for postleca in feca.iter_search_nodes(identity = 'post-LECA'):
        if len(postleca.search_nodes(identity = 'LECA')) > 0:
            postleca.identity = 'unknown'
    # Change rare unknown nodes that don't fulfill the duplication criterion but are duplications because there are duplications in both their children
    for unknown in feca.iter_search_nodes(identity = "unknown"):
        daughter1, daughter2 = unknown.get_children()
        if len(daughter1.search_nodes(identity = 'duplication')) > 0 and len(daughter2.search_nodes(identity = 'duplication')) > 0:
            unknown.identity = 'duplication'
            daughter1_sp = set([leaf.taxid for leaf in daughter1])
            daughter2_sp = set([leaf.taxid for leaf in daughter2])
            species_overlap = daughter1_sp & daughter2_sp
            all_species = daughter1_sp | daughter2_sp
            dupl_consistency = len(species_overlap) / len(all_species)
            unknown.add_features(overlap = len(species_overlap), consistency = dupl_consistency)

def get_prokaryotic_sister(euk_clade, tree, farthest):
    """Determines both possible prokaryotic sister groups in an unrooted way or a rooted way using the rooting on the farthest leaf"""
    sister = euk_clade.get_sisters() # Should be checked if there are any eukaryotic sequences in the sister group
    if len(sister) == 1:
        prok_leaves_sister = sister[0].search_nodes(prok_euk = 'Prokaryote')
    else: # In case of multifurcation: take all sisters (written out for clarity, but does not have to be split between bifurcating and multifurcating)
        prok_leaves_sister = []
        for sis in sister:
            prok_leaves_sister.extend(sis.search_nodes(prok_euk = 'Prokaryote'))
    if farthest:
        if len(tree) - (len(euk_clade) + len(prok_leaves_sister)) == 1: # So only 1 other non-sister sequence
            support = 'NA'
        else:
            support = euk_clade.up.support
        prok_taxids = [prok.taxid for prok in prok_leaves_sister]
        sp_tree = ncbi.get_topology(prok_taxids) # Get NCBI species tree, to get the identity of the LCA
        lca = sp_tree.taxid
        lca_name = ncbi.translate_to_names([lca])[0]
        return lca, lca_name, support
    else:
        other_prok_leaves = set(tree.search_nodes(prok_euk = 'Prokaryote')) - set(prok_leaves_sister)
        lcas = []
        if len(other_prok_leaves) == 1: # So only 1 other non-sister sequence
            supports = ['NA']
        else:
            supports = [euk_clade.up.support]
        if len(prok_leaves_sister) == 1:
                supports.append('NA')
        else:
            if len(sister) == 1:
                supports.append(sister[0].support)
            else:
                supports.append('NA')
        for i, group in enumerate([prok_leaves_sister, other_prok_leaves]):
            prok_taxids = []
            for prok in group: # Collect tax ids of prokaryotic sister leaves
                prok_taxids.append(prok.taxid)
            sp_tree = ncbi.get_topology(prok_taxids) # Get NCBI species tree, to get the identity of the LCA
            lca = sp_tree.taxid
            if lca == 1224: # Proteobacteria
                for proteo in sp_tree:
                    lineage = ncbi.get_lineage(proteo.name)
                    if not 28211 in lineage and not 1236 in lineage and not 28216 in lineage: # So, not an alpha/gamma/beta
                        lca_name = 'Proteobacteria'
                        break
                else: # So, only alpha/beta/gamma proteobacteria
                    lca = 'abgprot'
                    lca_name = 'ABG proteobacteria'
            elif lca == 2157: # Archaea
                for arch in sp_tree:
                    lineage = ncbi.get_lineage(arch.name)
                    if not 1935183 in lineage and not 1783275 in lineage: # So, not an Asgard or TACK
                        lca_name = 'Archaea'
                        break
                else: # So, only Asgards + TACK
                    lca = 'asgtack'
                    lca_name = 'Asgard+TACK group'
            else:
                lca_name = ncbi.translate_to_names([lca])[0]
            support = supports[i]
            lcas.append((lca, lca_name, support))
        return lcas

def classify_sister(lca): # Added alpha/beta/gamma proteo superclass and TACK+Asgard supersuperphylum
    """Classifies the prokaryotic sister-group"""
    if lca == 'abgprot':
        return 'ABG proteobacteria'
    elif lca == 'asgtack':
        return 'Asgard+TACK group'
    ancestors = ncbi.get_lineage(lca)
    if 28211 in ancestors:
        return 'Alphaproteobacteria'
    elif 1935183 in ancestors:
        return 'Asgard archaea'
    elif 28216 in ancestors:
        return 'Betaproteobacteria'
    elif 1236 in ancestors:
        return 'Gammaproteobacteria'
    elif 1783275 in ancestors:
        return 'TACK archaea'
    else:
        desired = 'phylum'
        if 1224 in ancestors:
            desired = 'class'
        ranks = ncbi.get_rank(ancestors)
        names = ncbi.get_taxid_translator(ancestors)
        for taxon in ranks: # Return phylum and if that is not present, then lowest group
            if ranks[taxon] == desired:
                return names[taxon]
        else:
            return names[ancestors[-1]]

def reroot(euk_clade, tree):
    """Reroots the tree on the farthest leaf from the eukaryotic clade"""
    tree.set_outgroup(euk_clade) # Root on this eukaryotic clade
    sister = euk_clade.get_sisters()[0]
    farthest = sister.get_farthest_leaf()[0]
    tree.set_outgroup(farthest) # Root on the leaf farthest from this eukaryotic clade (can be a false positive for example)
    return farthest.name

def reclassify_majority_sister(feca_clade):
    """Tries to reclassify the prokaryotic sister-group in case it was classified broadly as bacterial, archaeal or cellular"""
    sister_taxids = [leaf.taxid for sister in feca_clade.get_sisters() for leaf in sister if leaf.prok_euk == 'Prokaryote']
    sister_counts = {}
    for taxid in sister_taxids:
        ancestors = ncbi.get_lineage(taxid)
        if 1935183 in ancestors:
            sister_taxon = 'Asgard archaea'
        elif 1783275 in ancestors:
            sister_taxon = 'TACK archaea'
        else: # Use phylum
            ranks = ncbi.get_rank(ancestors)
            names = ncbi.get_taxid_translator(ancestors)
            desired = 'phylum'
            if 1224 in ancestors:
                desired = 'class'
            for taxon in ranks:
                if ranks[taxon] == desired:
                    sister_taxon = names[taxon]
                    break
            else: # In case no phylum annotation known
                sister_taxon = 'unclassified'
        sister_counts[sister_taxon] = sister_counts.get(sister_taxon, 0) + 1
    majority = len(sister_taxids) / 2
    for taxon in sister_counts:
        if sister_counts[taxon] > majority: # In case simple majority
            return taxon
    asgtack = sister_counts.get('Asgard archaea', 0) + sister_counts.get('TACK archaea', 0)
    if asgtack > majority: # Majority TACK or Asgard
        return 'Asgard+TACK group'
    abgprot_counts = 0
    for abgprot in ('Alphaproteobacteria', 'Betaproteobacteria', 'Gammaproteobacteria'):
        abgprot_counts += sister_counts.get(abgprot, 0)
    if abgprot_counts > majority: # Majority alhpa/beta/gamma proteobacteria
        return 'ABG proteobacteria'
    proteo_counts = abgprot_counts
    for proteo in ('Deltaproteobacteria', 'Epsilonproteobacteria', 'Zetaproteobacteria', 'Acidithiobacillia', 'Oligoflexia'):
        proteo_counts += sister_counts.get(proteo, 0)
    if proteo_counts > majority: # Majority proteobacteria
        return 'Proteobacteria'
    else:
        return False

def pick_branch_length(raw_branch_lengths, branch_lengths, ebls, mode = 'minimum', corresponding = False):
    """Picks the single branch length that corresponds to the minimum (default) or other value.
    If corresponding, the raw branch length is selected and the corresponding branch length and ebl."""
    if mode == 'median':
        rbl = median(raw_branch_lengths)
        bl = median(branch_lengths)
        ebl = median(ebls)
    elif mode == 'minimum':
        rbl = min(raw_branch_lengths)
        if not corresponding:
            bl = min(branch_lengths)
            ebl = min(ebls)
    elif mode == 'maximum':
        rbl = max(raw_branch_lengths)
        if not corresponding:
            bl = max(branch_lengths)
            ebl = max(ebls)
    elif mode == 'mean':
        rbl = mean(raw_branch_lengths)
        bl = mean(branch_lengths)
        ebl = mean(ebls)
    else:
        sys.stderr.write(f'Warning: mode "{mode}" for picking the branch lengths not detected!\n')
        return 'NA', 'NA', 'NA'
    if corresponding and mode in ('minimum', 'maximum'):
        bl = branch_lengths[raw_branch_lengths.index(rbl)]
        ebl = ebls[raw_branch_lengths.index(rbl)]
    return rbl, bl, ebl

def calculate_stem_lengths(tree, lecas, lepca, sister, mode = 'minimum', corresponding = False):
    """Calculates the branch lengths for a single FECA"""
    prok_branch_lengths = []
    for prok in sister:
        prok_branch_lengths.append(tree.get_distance(lepca, prok))
    prok_bl_med = median(prok_branch_lengths)
    raw_stem_lengths = []
    stem_lengths = []
    ebls = []
    for leca in lecas:
        euk_branch_lengths = []
        for euk in leca:
            euk_branch_lengths.append(tree.get_distance(leca, euk))
        euk_bl_med = median(euk_branch_lengths)
        rsl = tree.get_distance(lepca, leca)
        sl = rsl / euk_bl_med
        ebls.append(euk_bl_med)
        raw_stem_lengths.append(rsl)
        stem_lengths.append(sl)
    rsl, sl, ebl = pick_branch_length(raw_stem_lengths, stem_lengths, ebls, mode, corresponding)
    return prok_bl_med, rsl, sl, ebl

def calculate_duplication_lengths(tree, duplication, lecas, mode = "minimum", corresponding = False):
    """Calculates the duplications lengths. Note: these can be inconsistent!"""
    raw_dupl_lengths = []
    dupl_lengths = []
    ebls = []
    for leca in lecas:
        euk_branch_lengths = []
        for euk in leca:
            euk_branch_lengths.append(tree.get_distance(leca, euk))
        ebl_med = median(euk_branch_lengths)
        raw_dupl_length = tree.get_distance(leca, duplication)
        raw_dupl_lengths.append(raw_dupl_length)
        dupl_length = raw_dupl_length / ebl_med
        dupl_lengths.append(dupl_length)
        ebls.append(ebl_med)
    return pick_branch_length(raw_dupl_lengths, dupl_lengths, ebls, mode, corresponding)

def get_non_feca_sister(non_feca_node, tree): # In a rooted way (on farthest leaf)
    """Identifies the donor of the non-FECA clade"""
    tree.set_outgroup(non_feca_node)
    farthest_leaf = non_feca_node.get_sisters()[0].get_farthest_leaf()[0] # Might impact the results in case of multifurcation
    tree.set_outgroup(farthest_leaf)
    sister = non_feca_node.get_sisters()
    if len(sister) == 1:
        prok_leaves_sister = sister[0].search_nodes(prok_euk = 'Prokaryote')
    else:
        prok_leaves_sister = []
        for sis in sister:
            prok_leaves_sister.extend(sis.search_nodes(prok_euk = 'Prokaryote'))
    support = non_feca_node.up.support
    prok_taxids = []
    for prok in prok_leaves_sister:
        prok_taxids.append(prok.taxid)
    sp_tree = ncbi.get_topology(prok_taxids) # Get NCBI species tree, to get the identity of the LCA
    lca = sp_tree.taxid
    lca_name = ncbi.translate_to_names([lca])[0]
    fecas = []
    for sis in sister:
        for feca_clade in sis.iter_search_nodes(feca = True):
            fecas.append(feca_clade.feca_no)
    if len(fecas) > 0:
        sys.stderr.write('Warning: for non-FECA eukaryotic sequences the sister group does contain a FECA\n')
        lca_name += '+' + '+'.join(fecas)
    return lca, lca_name, support

def collect_leca_information(euk_clade, feca_count, leca_dict, representing = None, scrollsaw = True):
    for i, leca in enumerate(euk_clade.iter_search_nodes(identity = 'LECA')):
        leca_id = f'OG{feca_count}.{i+1}'
        leca.name = leca_id
        leca.add_face(TextFace(leca_id, fgcolor = 'green'), column = 0, position = 'branch-top')
        seqs = [leaf.name for leaf in leca.iter_leaves()]
        if scrollsaw:
            repr_seqs = [repr_seq for seq in seqs for repr_seq in representing.get(seq,'')]
            leca_dict[leca_id] = (str(leca.support), str(leca.coverage), str(leca.redundancy), ','.join(seqs), ','.join(repr_seqs))
        else:
            leca_dict[leca_id] = (str(leca.support), ','.join(seqs))

def collect_dupl_information(euk_clade, feca_count, dupl_dict, mode = 'minimum', corresponding = False):
    for i, duplication in enumerate(euk_clade.iter_search_nodes(identity = 'duplication')):
        dupl_id = f'D{feca_count}.{i+1}'
        duplication.name = dupl_id
        duplication.add_face(TextFace(dupl_id, fgcolor = 'green'), column = 0, position = 'branch-top')
        dupl_lecas = duplication.search_nodes(identity = 'LECA')
        rdl, dl, ebl = calculate_duplication_lengths(euk_clade, duplication, dupl_lecas, mode = mode, corresponding = corresponding) # Tree --> euk_clade
        support = duplication.support
        if support == 101.0:
            support = 'NA'
        children = duplication.get_children()
        ogs1 = [leca.name for leca in children[0].iter_search_nodes(identity = 'LECA')]
        ogs2 = [leca.name for leca in children[1].iter_search_nodes(identity = 'LECA')]
        dupl_dict[dupl_id] = (support, duplication.overlap, duplication.consistency, ','.join(ogs1) + ' - ' + ','.join(ogs2), rdl, dl, ebl)

def collect_unknown_information(euk_clade, feca_count, unknown_dict, representing, scrollsaw = True):
    for i, unknown in enumerate(euk_clade.iter_search_nodes(identity = 'unknown')):
        unknown_id = f'U{feca_count}.{i+1}'
        unknown.name = unknown_id
        for child in unknown.get_children():
            if child.identity == 'duplication' or child.identity == 'unknown': # Only want the child that makes this node 'unknown'
                continue
            seqs = [leaf.name for leaf in child.iter_leaves()]
            if scrollsaw:
                unknown_dict[unknown_id] = (str(unknown.support), str(child.coverage), ','.join(seqs))
            else:
                unknown_dict[unknown_id] = (str(unknown.support), ','.join(seqs))
