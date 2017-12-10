from Bio import SeqIO
import os
import numpy as np
from numpy.random import randint


def get_test_reads(file_list, kmer_len, num_reads):
    print("test genome file: {}".format(file_list))
    unclassified_reads = []
    seq_count = 0
    accession_num = 0
    for file in file_list:
        for accession in SeqIO.parse(file, 'genbank'):
            seq_count += 1

        seqs_from_each = int(num_reads/seq_count)
        print("seqs from each: {}".format(seqs_from_each))
        seqs_from_last = num_reads - seqs_from_each * seq_count + num_reads%seqs_from_each
        print("seqs from last: {}".format(seqs_from_last))

        for accession in SeqIO.parse(file, 'genbank'):
            if accession_num < seq_count:
                for i in range(0, seqs_from_each):
                    len_seq = len(accession.seq)
                    rand_seq_start = randint(0, len_seq - kmer_len)
                    unclassified_reads.append(str(accession.seq[rand_seq_start: rand_seq_start + kmer_len]))
                accession_num += 1
            elif seqs_from_last > 0:
                for i in range(0, seqs_from_last):
                    len_seq = len(accession.seq)
                    rand_seq_start = randint(0, len_seq - kmer_len)
                    unclassified_reads.append(str(accession.seq[rand_seq_start: rand_seq_start + kmer_len]))
    return unclassified_reads


def get_test_tax(file_list):
    for file in file_list:
        for accession in SeqIO.parse(file, 'genbank'):
            true_tax = accession.annotations['taxonomy']
    return true_tax


def main():
    # get one exact test sequence from each
    kmer_len = 150
    num_reads = 5

    file_list = []
    for file_name in os.listdir("."):
        if file_name.endswith("test_.gbff"):
            file_list.append(file_name)

    unclassified_reads = get_test_reads(file_list, kmer_len, num_reads)
    get_test_tax(file_list)


if __name__ == "__main__":
    main()






