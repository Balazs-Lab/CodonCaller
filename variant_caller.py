import pysam
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
from multiprocessing import Pool
from collections import defaultdict, OrderedDict
import csv

codon_to_amino_acid = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def write_to_csv(frequency_dictionary, csv_file_path='../test_data/aa_frequencies.csv'):
    # Get a set of all unique codon types
    all_aminoacids = set(codon_to_amino_acid.values())

    # Create a CSV file and write headers
    with open(csv_file_path, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)

        # Write headers (position and all codon types)
        headers = ['Position'] + list(all_aminoacids)
        csv_writer.writerow(headers)

        # Write frequency data for each position
        for position, freq_data in frequency_dictionary.items():
            row_data = [position] + [freq_data.get(aminoacid, 0) for aminoacid in all_aminoacids]
            csv_writer.writerow(row_data)

    print(f'CSV file saved at: {csv_file_path}')


class VariantCaller:
    def __init__(self, bed_file, bam_file, reference_fasta, read_quality_threshold=20, base_quality_threshold=35):
        self.reference_fasta = None
        self.bed_file = bed_file
        self.bam_file = bam_file
        self.fasta_file = reference_fasta
        self.read_quality_threshold = read_quality_threshold
        self.base_quality_threshold = base_quality_threshold * 3
        self.reference_sequence = self.read_fasta()
        self.cds_regions = self.load_cds_regions()
        self.cds_aminoacids = self.translate_cds()
        self.reads = self.load_reads()
        self.aa_counts = {}
        self.aa_freq = {}

    def load_cds_regions(self):
        """
        :return: dict of CDS regions - name, start, stop
        """
        cds_regions = dict()
        with open(self.bed_file, 'r') as bed:
            for line in bed:
                chrom, start, end = line.strip().split('\t')
                cds_regions[chrom] = (int(start) - 1, int(end) - 1)
        return cds_regions

    def read_fasta(self):
        """
        :return: SeqRecord of reference sequence
        """
        record = SeqIO.read(self.fasta_file, "fasta")
        return record

    def translate_cds(self):
        """
        :return: dict with translations of all CDS regions / annotations in the reference sequence
        """
        cds_aminoacids = dict()
        for key, value in self.cds_regions.items():
            sequence = self.reference_sequence
            cds_start = self.cds_regions[key][0]
            cds_end = self.cds_regions[key][1]
            cds_sequence = sequence[cds_start:cds_end]  # Assuming 1-based indexing
            cds_aminoacids[key] = cds_sequence.translate()
        return cds_aminoacids

    def load_reads(self):
        """
        load mapped reads
        :return:
        """
        with pysam.AlignmentFile(self.bam_file, "rb") as bam:
            reads = list(bam.fetch())
        return reads

    def filter_amino_acid(self, pairs):
        key, value = pairs
        # amino acid value - codon with an "N" maps to Unknown
        if value[0] == "Unknown":
            return False
        # codon sum score - sum of all 3 nucleotide mapping scores
        elif value[1] < self.read_quality_threshold:
            return False
        else:
            return True

    def parse_read(self, read):
        """

        :param read: perform operation on a specific read
        :return:component sequence, position, and quality
        """

        read_sequence = Seq(read.query_sequence)
        read_positions = read.get_reference_positions()
        quality = read.query_alignment_qualities

        return read_sequence, read_positions, quality

    def read_to_codon(self, cds, read_sequence, read_positions, quality):
        """
        :param quality:
        :param read_positions:
        :param read_sequence:
        :param cds: specify the coding region to define the cds frame
        :param read: perform operation on a specific read
        :return:
        """

        cds_start = self.cds_regions[cds][0]
        # cds_end = self.cds_regions[cds][1]

        if read_positions[0] < cds_start:
            adjust = cds_start - read_positions[0]
            read_positions = read_positions[adjust:]
            read_sequence = read_sequence[adjust:]
            quality = quality[adjust:]

        # set the frame based on current position
        if (read_positions[0] - cds_start) % 3 == 0:
            codonshift: int = 0
        elif (read_positions[0] - cds_start) % 3 == 1:
            codonshift: int = 2
        elif (read_positions[0] - cds_start) % 3 == 2:
            codonshift: int = 1

        # adjust the start of the sequence
        adjusted_positions = [pos + codonshift for pos in read_positions]

        codon_positions = [1 + ((i - cds_start + codonshift) // 3) for i in adjusted_positions[0::3] if
                           i >= cds_start]

        # adjust the end of the sequence
        if len(read_sequence[codonshift:]) % 3 == 0:
            end_shift = 0
            adjusted_sequence = read_sequence[codonshift: len(read_sequence) - end_shift]
            adjusted_quality = quality[codonshift: len(quality) - end_shift]
            codon_positions = codon_positions[:]

        elif len(read_sequence[codonshift:]) % 3 == 1:
            end_shift = 1
            adjusted_sequence = read_sequence[codonshift: len(read_sequence) - end_shift]
            adjusted_quality = quality[codonshift: len(quality) - end_shift]
            codon_positions = codon_positions[:-1]

        elif len(read_sequence[codonshift:]) % 3 == 2:
            end_shift = 2
            adjusted_sequence = read_sequence[codonshift: len(read_sequence) - end_shift]
            adjusted_quality = quality[codonshift: len(quality) - end_shift]
            codon_positions = codon_positions[:-2]

        codons = [adjusted_sequence[i:i + 3] for i in range(0, len(adjusted_sequence), 3)]
        quality_sum = [sum(adjusted_quality[i:i + 3]) for i in range(0, len(adjusted_quality), 3)]

        return codon_positions, codons, quality_sum

    def process_read(self, cds, read):

        read_sequence, read_positions, quality = self.parse_read(read)
        codon_positions, codons, quality_sum = self.read_to_codon(cds, read_sequence, read_positions, quality)

        # convert codons to amino acids
        amino_acids = [codon_to_amino_acid.get(codon, 'Unknown') for codon in codons]

        d = dict(zip(codon_positions, zip(amino_acids, quality_sum)))

        ## filter out Ns and low quality mapping codons
        d_clean = dict(filter(self.filter_amino_acid, d.items()))

        # print(d_clean)
        return d_clean

    def process_sample(self, cds):

        count = 0
        read_data = []
        for read in self.reads:
            if (read.mapping_quality > self.read_quality_threshold) and (read.reference_start > 0):
                result = self.process_read(cds, read)
                read_data.append(result)

        compiled_dataset = defaultdict(lambda: defaultdict(int))
        for dict_data in read_data:
            try:
                for position, codon_data in dict_data.items():
                    if position not in compiled_dataset.keys():
                        compiled_dataset[position] = [codon_data]
                    else:
                        compiled_dataset[position].extend([codon_data])
            except AttributeError:
                pass

        ordered_dict = OrderedDict(sorted(compiled_dataset.items()))

        # Dictionary to store codon frequencies for each position
        # self.codon_frequencies = {}

        # Iterate through each position in the mapped_positions dictionary
        for position, codon_metadata_list in ordered_dict.items():
            # Dictionary to store codon frequencies for the current position
            position_codon_counts = {}

            # Iterate through each codon-metadata tuple in the list
            for codon, metadata in codon_metadata_list:
                # Count the frequency of each codon for the current position
                position_codon_counts[codon] = position_codon_counts.get(codon, 0) + 1

            # Store the codon frequencies for the current position in the result dictionary
            self.aa_counts[position] = position_codon_counts

        # Now, codon_frequencies is a dictionary where the key is the position,
        # and the value is another dictionary containing codon frequencies
        # print(self.codon_frequencies)

    def aa_count_to_freq(self, cds):

        cds_start = self.cds_regions[cds][0]
        cds_end = self.cds_regions[cds][1]
        reference_protein = self.cds_aminoacids[cds]

        reference_dict = {i + 1: reference_protein[i] for i in range(len(reference_protein))}
        common_positions = set(reference_dict.keys()).intersection(self.aa_counts.keys())

        positions = []
        mutant_freq = []

        # Compare amino acids at each position
        for position in common_positions:
            reference_aa = reference_dict[position]
            experimental_aa = self.aa_counts[position]
            wt = 0
            mut = 0
            tmp_dict = dict()
            for key, value in experimental_aa.items():
                if key == reference_aa:
                    wt = value
                else:
                    mut += value
            for key, value in experimental_aa.items():
                if key != reference_aa:
                    tmp_dict[key] = value / (mut + wt)
            tmp_dict["WT"] = wt / (wt + mut)

            mutant_freq.append(tmp_dict)
            positions.append(position)

        self.aa_freq = dict(zip(positions, mutant_freq))


    def write_to_csv(self, csv_file_path='../test_data/aa_frequencies.csv'):

        # Get a set of all unique codon types
        all_aminoacids = set(codon_to_amino_acid.values())

        # Create a CSV file and write headers
        with open(csv_file_path, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)

            # Write headers (position and all codon types)
            headers = ['Position'] + list(all_aminoacids)
            csv_writer.writerow(headers)

            # Write frequency data for each position
            for position, freq_data in self.aa_freq.items():
                row_data = [position] + [freq_data.get(aminoacid, 0) for aminoacid in all_aminoacids]
                csv_writer.writerow(row_data)

        print(f'CSV file saved at: {csv_file_path}')
