# This script compares k-mers from sample read to reference database
# Requires sample file in fastq format
# Outputs percent identity to known reference genomes

from Bio import SeqIO
from itertools import islice
import json
import os
import numpy as np
import pydot
# from get_test_short_reads import get_test_reads, get_test_tax
from get_test_reads_sra import get_test_reads, get_test_tax
from build_ref_dict import build_ref_tax_dict, build_ref_kmer_dict_no_overlap, write_to_file
from operator import itemgetter


# a class containing functions to retrieve unclassified test reads
class GetFiles(object):

    def __init__(self):
        pass

    @staticmethod
    # reads must be in the working directory
    def get_with_fasta_filename(filename):
        unclassified_reads = SeqIO.parse(filename, "fastq")
        return unclassified_reads

    @staticmethod
    # reads must be in the SRA_Test_Sequences directory under the working directory
    def get_from_directory(num_test_seqs, file):
        unclassified_reads = get_test_reads(file, 150, num_test_seqs)
        true_tax = get_test_tax(file)
        return unclassified_reads, true_tax


# this function returns the complement of a DNA sequence
def complement_DNA_seq(k_mer):
    opposite_k_mer = ""
    for base in range(0, len(k_mer)):
        if k_mer[base] == 'A':
            opposite_k_mer = opposite_k_mer + 'T'
        elif k_mer[base] == 'T':
            opposite_k_mer = opposite_k_mer + 'A'
        elif k_mer[base] == 'C':
            opposite_k_mer = opposite_k_mer + 'G'
        elif k_mer[base] == 'G':
            opposite_k_mer = opposite_k_mer + 'C'
        else:
            opposite_k_mer = opposite_k_mer + 'N'
    return opposite_k_mer


# this function counts the number of occurrences of a k-mer from a test read and its reverse, complement,
# and reverse complement in each reference organism
def count_k_mers(k, seq, ref_k_mer_dict):
    common_k_mer_dict = {}
    # for every k-mer in sample seq, check for matches in reference dictionary
    for s in range(0, len(seq) - k + 1):
        k_mer = str(seq[s:s + k])
        reverse_k_mer = k_mer[::-1]
        opposite_k_mer = complement_DNA_seq(k_mer)
        opposite_rev_k_mer = complement_DNA_seq(reverse_k_mer)
        if k_mer in ref_k_mer_dict:  # check for k-mer presence in reference genome dictionary
            organism_count_dict = ref_k_mer_dict[k_mer]
            for seq_id in organism_count_dict:
                if seq_id not in common_k_mer_dict:
                    common_k_mer_dict[seq_id] = 0
                common_k_mer_dict[seq_id] += organism_count_dict[seq_id]
        if reverse_k_mer in ref_k_mer_dict:  # check for reverse of k-mer too
            organism_count_dict = ref_k_mer_dict[reverse_k_mer]
            for seq_id in organism_count_dict:
                if seq_id not in common_k_mer_dict:
                    common_k_mer_dict[seq_id] = 0
                common_k_mer_dict[seq_id] += organism_count_dict[seq_id]
        if opposite_k_mer in ref_k_mer_dict:  # check for reverse of k-mer too
            organism_count_dict = ref_k_mer_dict[opposite_k_mer]
            for seq_id in organism_count_dict:
                if seq_id not in common_k_mer_dict:
                    common_k_mer_dict[seq_id] = 0
                common_k_mer_dict[seq_id] += organism_count_dict[seq_id]
        if opposite_rev_k_mer in ref_k_mer_dict:  # check for reverse of k-mer too
            organism_count_dict = ref_k_mer_dict[opposite_rev_k_mer]
            for seq_id in organism_count_dict:
                if seq_id not in common_k_mer_dict:
                    common_k_mer_dict[seq_id] = 0
                common_k_mer_dict[seq_id] += organism_count_dict[seq_id]
    return common_k_mer_dict


# this function determines what taxonomic label to assign to each unclassified read
def classify(seq_num, common_k_mer_dict, ref_tax_dict, tax_class_dict, unclassified, req_hits):

    # remove ref genomes with fewer than the specified # of k-mer matches
    acc_ids = list(common_k_mer_dict)
    for acc_id in acc_ids:
        if common_k_mer_dict[acc_id] < req_hits:
            del common_k_mer_dict[acc_id]
    acc_ids = list(common_k_mer_dict)

    # if hits found only in one reference genome:
    if len(acc_ids) == 1:
        # this is currently in memory & can get very large depending on number of reference & test sequences
        tax_class_dict[seq_num] = ref_tax_dict[acc_ids[0]]
        # to classify only if matches are found in multiple reference genomes:
        # unclassified += 1

    # if hits found in multiple reference genomes:
    elif len(acc_ids) > 1:
        # 1) To classify by LCA:
        # taxonomy = []
        # for acc_id in acc_ids:
        #     taxonomy.append(ref_tax_dict[acc_id])
        # n = 6  # species level
        # while len(taxonomy) > 1:
        #     for i in range(0, len(taxonomy)):
        #         taxonomy[i] = taxonomy[i][:n]
        #     j = 0
        #     while j < len(taxonomy) - 1:
        #         if taxonomy[j] == taxonomy[j+1]:
        #             del(taxonomy[j+1])
        #         else:
        #             j += 1
        #     n -= 1
        #     tax_class_dict[seq_num] = taxonomy

        # 2) To classify by most k-mers shared
        # key_value_tuples = common_k_mer_dict.items()
        # max_match_id = max(key_value_tuples, key=itemgetter(1))[0]
        # tax_class_dict[seq_num] = ref_tax_dict[max_match_id]

        # 3) To classify based on LCA from within 0.5 SD of max # k-mers shared
        half_std = np.std(list(common_k_mer_dict.values()))/2
        key_value_tuples = common_k_mer_dict.items()
        sorted_tuples = sorted(key_value_tuples, key=lambda x: x[1])
        tuples_above_std = filter(lambda x: x[1] > half_std, sorted_tuples)
        filtered_common_k_mer_dict = dict(list(tuples_above_std))

        taxonomy = []
        acc_ids = list(filtered_common_k_mer_dict)
        for acc_id in acc_ids:
            taxonomy.append(ref_tax_dict[acc_id])
        n = 6  # genus level
        while len(taxonomy) > 1:
            for i in range(0, len(taxonomy)):
                taxonomy[i] = taxonomy[i][:n]
            j = 0
            while j < len(taxonomy) - 1:
                if taxonomy[j] == taxonomy[j+1]:
                    del(taxonomy[j+1])
                else:
                    j += 1
            n -= 1
        tax_class_dict[seq_num] = taxonomy[0]

    # count reads with no hits above required threshold as unclassified
    else:
        unclassified += 1

    return unclassified


# this function writes the precision and specificity of classifications to a file
def write_results(true_tax, tax_class_dict, total_reads, k, req_hits, num_test_seqs):
    # check accuracy
    class_correct = {'Domain': 0, 'Phylum': 0, 'Order': 0, 'Family': 0, 'Genus': 0, 'Species': 0}
    class_incorrect = {'Domain': 0, 'Phylum': 0, 'Order': 0, 'Family': 0, 'Genus': 0, 'Species': 0}
    total_correct = 0
    total_incorrect = 0
    total_specific_correct = 0
    # this assumes annotations('taxonomy') is consistently the same taxonomic levels:
    index_to_taxa = {0: 'Domain', 1: 'Phylum', 2: 'Order', 3: 'Family', 4: 'Genus', 5: 'Species'}

    for key in tax_class_dict:
        # for correctly classified reads
        if tax_class_dict[key][-1] in true_tax:
            # print(tax_class_dict[key][-1])
            tax_index = true_tax.index(tax_class_dict[key][-1])
            try:
                class_correct[index_to_taxa[tax_index]] += 1
            except KeyError:
                print(tax_class_dict[key])
                print("does not fit taxonomy format 'Domain, Phylum, Order, Family, Genus, Species'")
            total_correct += 1
            if tax_class_dict[key][-1] != true_tax:
                total_specific_correct += 1
        # for incorrectly classified reads
        else:
            tax_index = tax_class_dict[key].index(tax_class_dict[key][-1])
            try:
                class_incorrect[index_to_taxa[tax_index]] += 1
            except KeyError:
                print(tax_class_dict[key])
                print("does not fit taxonomy format 'Domain, Phylum, Order, Family, Genus, Species'")
            total_incorrect += 1

    try:
        precision = total_correct*100/(total_correct + total_incorrect)
    except ZeroDivisionError:
        precision = 0
    sensitivity_domain = (class_correct['Domain'] + class_correct['Phylum'] + class_correct['Order'] + class_correct[
        'Family'] + class_correct['Genus'] + class_correct['Species'])*100/total_reads
    sensitivity_phylum = (class_correct['Phylum'] + class_correct['Order'] + class_correct['Family'] + class_correct[
        'Genus'] + class_correct['Species']) * 100 / total_reads
    sensitivity_order = (class_correct['Order'] + class_correct['Family'] + class_correct['Genus'] + class_correct[
        'Species']) * 100 / total_reads
    sensitivity_family = (class_correct['Family'] + class_correct['Genus'] + class_correct[
        'Species']) * 100 / total_reads
    sensitivity_genus = (class_correct['Genus'] + class_correct['Species']) * 100 / total_reads
    sensitivity_species = (class_correct['Species']) * 100 / total_reads

    # write results to file "classification_results.txt"
    file = open("classification_results.txt", 'a')
    file.write('organism: {3}, k = {0}, required k-mer hits = {1}, # of sequences tested = {2}\n\n'.format(k, req_hits, num_test_seqs, true_tax[-1]))
    file.write("{0:.2f}% precision\n".format(precision))
    file.write('{0:.2f}% sensitivity at the domain level\n'.format(sensitivity_domain))
    file.write('{0:.2f}% sensitivity at the phylum level\n'.format(sensitivity_phylum))
    file.write('{0:.2f}% sensitivity at the order level\n'.format(sensitivity_order))
    file.write('{0:.2f}% sensitivity at the family level\n'.format(sensitivity_family))
    file.write('{0:.2f}% sensitivity at the genus level\n'.format(sensitivity_genus))
    file.write('{0:.2f}% sensitivity at the species level\n\n'.format(sensitivity_species))
    file.close()


# this function adds test organism taxonomy to the tree figure
def add_to_tree_dict(tree, path):
    if len(path) == 0:
        return
    if path[0] not in tree:
        tree[path[0]] = {}
    add_to_tree_dict(tree[path[0]], path[1:])
    return tree


# this function prints the nested dictionary form of the tree figure to the command line
def print_tree(tree, offset=0):
    for k in tree.keys():
        print("*"*offset + k)
        print_tree(tree[k], offset=offset+1)


# this function creates a nested dictionary tree of reference genome taxonomic relationships
def create_tree(ref_tax_dict):
    # create tree based on first reference genome tax_list
    keys = list(ref_tax_dict)
    tax_list = ref_tax_dict[keys[0]]
    value = {tax_list[-1]: {}}
    for i in range(1, len(tax_list)):
        value2 = {tax_list[-1 - i]: value}
        value = value2
    tree = value

    # add other genomes to tree
    for key in range(1, len(keys)):
        tree = add_to_tree_dict(tree, ref_tax_dict[keys[key]])
    return tree


# this function draws the tree figure (see visit function)
def draw(parent_name, child_name, graph):
    edge = pydot.Edge(parent_name, child_name)
    graph.add_edge(edge)


# this function records how many test reads are classified to a certain node in the tree on the figure
# (see add_classified_reads_to_graph function)
def update_node_label(node_name, number, graph, style):
    if number != 'NA':
        node = pydot.Node(node_name, label=node_name+': '+str(number), style=style)
        graph.add_node(node)
    else:
        node = pydot.Node(node_name, style=style)
        graph.add_node(node)
    return graph


# this function draws the tree figure (see draw function)
def visit(tree, graph, parent=None):
    for k in tree.keys():
        if isinstance(tree[k], dict):
            if parent:
                draw(parent, k, graph)
            visit(tree[k], graph, k)
        else:
            draw(parent, k, graph)
            # drawing the label using a distinct name
            draw(k, k+'_'+tree[k], graph)


# this function saves the tree figure to a file
def graph_tree(tree, graph, filename):
    visit(tree, graph)
    graph.write_png(filename)
    return graph


# this function records how many test reads are classified to a certain node in the tree on the figure
# (see update_node_label function)
def add_classified_reads_to_graph(tax_class_dict, graph):
    classification_tally = {}
    for key in list(tax_class_dict):
        node_name = tax_class_dict[key][-1]
        if node_name not in classification_tally:
            classification_tally[node_name] = 1
        else:
            classification_tally[node_name] += 1
    # among all test sequences tally how many classified to each taxonomy (node_name)
    for key in list(tax_class_dict):
        node_name = tax_class_dict[key][-1]
        number = classification_tally[node_name]
        graph = update_node_label(node_name, number, graph, style='filled')
    return graph


def main():
    # define parameters
    try:
        k = int(raw_input("Enter k-mer size, as integer: "))
        req_hits = int(
            raw_input("Enter number of k-mer matches to require, as integer: "))
        num_test_seqs = int(
            raw_input("Enter number of test sequences to generate, as integer: "))

    except NameError:
        k = int(input("Enter k-mer size, as integer: "))
        req_hits = int(input("Enter number of k-mer matches to require, as integer: "))
        num_test_seqs = int(
            input("Enter number of test sequences to generate, as integer: "))

    print("Remember to rebuild k-mer dict if switching training organism set!")
    try:
        print("...loading k-mer dictionary")
        file = open("Ref_kmer_dict" + str(k) + ".txt", 'r')
        dict_string = file.read()
        ref_kmer_dict = json.loads(dict_string)
        file.close()
        # read taxonomic info for reference sequences from file (created using build_ref_dict.py script)
        print("...loading taxonomic dictionary")
        file = open("Ref_tax_dict" + str(k) + ".txt", 'r')
        tax_string = file.read()
        ref_tax_dict = json.loads(tax_string)
        file.close()

    except IOError or FileNotFoundError:
        ref_tax_dict, ref_genome_dict = build_ref_tax_dict()
        ref_kmer_dict = build_ref_kmer_dict_no_overlap(ref_genome_dict, ref_tax_dict, k)
        write_to_file(ref_kmer_dict, ref_tax_dict, k)

    # get test read file names
    print("...getting short test sequence samples")
    file_list = []
    for file_name in os.listdir("SRA_Test_Sequences"):
        if file_name.endswith("fastq"):
            file_name = os.path.join("SRA_Test_Sequences", file_name)
            file_list.append(file_name)

    # create tree and graph with reference organisms
    tree = create_tree(ref_tax_dict)
    graph = pydot.Dot(graph_type='graph')
    graph_tree(tree, graph, filename='Figures/Graph_of_ref_organisms.png')
    # generate list of taxa in reference set
    taxa_in_ref = []
    for key in list(ref_tax_dict):
        for entry in ref_tax_dict[key]:
            taxa_in_ref.append(entry)
    unclassified_taxa = []

    for file in file_list:
        unclassified_reads, true_tax = GetFiles.get_from_directory(num_test_seqs, file)
        total_reads = len(unclassified_reads)

        # generate list of taxa in test set
        for entry in true_tax:
            unclassified_taxa.append(entry)

        # classify taxonomy of sample reads
        print("...classifying reads")
        tax_class_dict = {}
        seq_num = 0
        unclassified = 0
        for read in unclassified_reads:
            # seq = read.seq # for when parsing fastq file
            seq = read # for when input is gbff file
            common_k_mer_dict = count_k_mers(k, seq, ref_kmer_dict)
            classify(seq_num, common_k_mer_dict, ref_tax_dict, tax_class_dict, unclassified, req_hits)
            seq_num += 1
        # print("tax class dict: {}".format(tax_class_dict))

        # check accuracy and write results to file
        print("...checking accuracy")
        write_results(true_tax, tax_class_dict, total_reads, k, req_hits, num_test_seqs)

        # add newly classified reads to graph node labels
        tree = add_to_tree_dict(tree, true_tax)
        graph = add_classified_reads_to_graph(tax_class_dict, graph)

    # add test genomes to figure and denote nodes with dotted lines
    print("...graphing")
    for i in range(0, len(unclassified_taxa)):
        if unclassified_taxa[i] not in taxa_in_ref:
            graph = update_node_label(unclassified_taxa[i], 'NA', graph, style='dotted')
    graph_tree(tree, graph, filename='Figures/Graph_with_classified_reads.png')


if __name__ == "__main__":

    main()
