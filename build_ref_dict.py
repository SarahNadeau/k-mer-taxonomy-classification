# This script builds a dictionary of k_mers from reference organism genomes

from Bio import SeqIO
import os
import json
from math import floor


# This function builds a tree of reference genome k-mers w/
# where each node is the set of k-mers existing in all the members of that taxa
def build_ref_tax_dict():
    print("...building reference taxonomy dictionary")
    file_list = []
    for file_name in os.listdir("Genomes"):
        if file_name.endswith(".gbff") and not file_name.endswith("test.gbff"):
            file_name = os.path.join("Genomes", file_name)
            file_list.append(file_name)
    # make dictionary of reference organism sequence objects labeled by seq_id
    ref_genome_dict = {}
    ref_tax_dict = {}
    for file in file_list:
        for accession in SeqIO.parse(file, 'genbank'):
            ref_genome_dict[accession.id] = accession.seq
            taxonomy = accession.annotations['taxonomy']
            taxonomy.append(accession.annotations['organism'])
            ref_tax_dict[accession.id] = taxonomy
    return ref_tax_dict, ref_genome_dict

def build_ref_kmer_dict_overlapping(ref_genome_dict, ref_tax_dict, k):

    # ref_kmer_dict keys are k-mers and values are NCBI seq-ids
    ref_kmer_dict = {}
    print("...adding genomes to reference k-mer dictionary")
    for seq_id in ref_genome_dict:
        # print("Adding {} to reference genome k-mer dictionary.".format(ref_tax_dict[seq_id]))
        for s in range(0, len(ref_genome_dict[seq_id]) - k + 1):
            k_mer = str(ref_genome_dict[seq_id][s:s+k])
            if k_mer not in ref_kmer_dict:  # if this is 1st occurrence of the k-mer, add empty dict. entry for org.
                ref_kmer_dict[k_mer] = {seq_id: 0}
            if seq_id not in ref_kmer_dict[k_mer]:  # if this is 1st occurrence of org. for the k-mer, add the org.
                ref_kmer_dict[k_mer][seq_id] = 0
            ref_kmer_dict[k_mer][seq_id] += 1  # add this k-mer occurrence to the total count for the org.

    return ref_kmer_dict

def build_ref_kmer_dict_no_overlap(ref_genome_dict, ref_tax_dict, k):
    print("...adding genomes to reference k-mer dictionary")
    # ref_kmer_dict keys are k-mers and values are NCBI seq-ids
    ref_kmer_dict = {}
    for seq_id in ref_genome_dict:
        # print("Adding {} to reference genome k-mer dictionary.".format(ref_tax_dict[seq_id]))
        for s in range(0, int(floor(len(ref_genome_dict[seq_id])/k))):
            k_mer = str(ref_genome_dict[seq_id][k * s : k * s + k])
            if k_mer not in ref_kmer_dict:  # if this is 1st occurrence of the k-mer, add empty dict. entry for org.
                ref_kmer_dict[k_mer] = {seq_id: 0}
            if seq_id not in ref_kmer_dict[k_mer]:  # if this is 1st occurrence of org. for the k-mer, add the org.
                ref_kmer_dict[k_mer][seq_id] = 0
            ref_kmer_dict[k_mer][seq_id] += 1  # add this k-mer occurrence to the total count for the org.

    return ref_kmer_dict


def write_to_file(ref_k_mer_dict, ref_tax_dict, k):
    # make k_mer dictionary into string using json
    print("...writing k_mer dictionary to file")
    dict_string = json.dumps(ref_k_mer_dict)
    # write string representation of dictionary to a file
    file = open("Ref_kmer_dict" + str(k) + ".txt", 'w')
    file.write(dict_string)
    file.close()

    # make taxonomy dictionary into string using json
    print("...writing taxonomy dictionary to file")
    tax_string = json.dumps(ref_tax_dict)
    file = open("Ref_tax_dict" + str(k) + ".txt", 'w')
    file.write(tax_string)
    file.close()


def main():
    try:
        k = int(raw_input("Enter kmer size as integer: "))
    except NameError:
        k = int(input("Enter kmer size as integer: "))
    ref_tax_dict, ref_genome_dict = build_ref_tax_dict()
    ref_k_mer_dict = build_ref_kmer_dict_no_overlap(ref_genome_dict, ref_tax_dict, k)
    # if overlapping k-mers desired:
    # ref_k_mer_dict = build_ref_kmer_dict_overlapping(ref_genome_dict, ref_tax_dict, k)
    write_to_file(ref_k_mer_dict, ref_tax_dict, k)


if __name__ == "__main__":
    main()



