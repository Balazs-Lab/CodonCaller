# import pysam
# import pandas as pd
# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqUtils import seq1
# from multiprocessing import Pool
# from collections import defaultdict, OrderedDict
# import csv
#
# codon_to_amino_acid = {
#             'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
#             'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
#             'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
#             'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
#             'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
#             'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
#             'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
#             'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
#             'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
#             'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
#             'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
#             'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
#             'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
#             'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
#             'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
#             'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
#         }
#
# class AminoAcidAnalyzer:
#     def __init__(self, bed_file, bam_file, fasta_file, quality_threshold=20):
#         self.bed_file = bed_file
#         self.bam_file = bam_file
#         self.fasta_file = fasta_file
#         self.quality_threshold = quality_threshold
#         self.cds_regions = self.load_cds_regions()
#         self.reference_sequence = self.read_fasta()
#
#     def load_cds_regions(self):
#         cds_regions = []
#         with open(self.bed_file, 'r') as bed:
#             for line in bed:
#                 chrom, start, end = line.strip().split('\t')
#                 cds_regions.append((chrom, int(start), int(end)))
#         return cds_regions
#
#     def read_fasta(self):
#         record = SeqIO.read(self.fasta_file, "fasta")
#         return record
#
#     def translate_cds(self, sequence, cds_start, cds_end):
#         cds_sequence = sequence[cds_start - 1:cds_end]  # Assuming 1-based indexing
#         amino_acids = cds_sequence.translate()
#         return amino_acids
#
#     def read_to_codon(self, codonstart, refpositions, sequence, quality, codonend):
#
#         if refpositions[0] < codonstart:
#             adjust = codonstart-refpositions[0]
#             refpositions = refpositions[adjust:]
#             sequence = sequence[adjust:]
#             quality = quality[adjust:]
#
#         if (refpositions[0] - codonstart) % 3 == 0:
#             codonshift: int = 0
#         elif (refpositions[0] - codonstart) % 3 == 1:
#             codonshift: int = 2
#         elif (refpositions[0] - codonstart) % 3 == 2:
#             codonshift:int  = 1
#
#         # adjust start position
#         adjusted_positions = [pos + codonshift for pos in refpositions]
#
#         codon_positions = [1 + ((i - codonstart + codonshift) // 3) for i in adjusted_positions[0::3] if
#                            i >= codonstart]
#
#         # fix ends of the sequence
#         if len(sequence[codonshift:]) % 3 == 0:
#             endshift = 0
#             adjusted_sequence = sequence[codonshift: len(sequence) - endshift]
#             adjusted_quality = quality[codonshift: len(quality) - endshift]
#             codon_positions = codon_positions[:]
#
#         elif len(sequence[codonshift:]) % 3 == 1:
#             endshift = 1
#             adjusted_sequence = sequence[codonshift: len(sequence) - endshift]
#             adjusted_quality = quality[codonshift: len(quality) - endshift]
#             codon_positions = codon_positions[:-1]
#
#         elif len(sequence[codonshift:]) % 3 == 2:
#             endshift = 2
#             adjusted_sequence = sequence[codonshift: len(sequence) - endshift]
#             adjusted_quality = quality[codonshift: len(quality) - endshift]
#             codon_positions = codon_positions[:-2]
#
#
#         codons = [adjusted_sequence[i:i + 3] for i in range(0, len(adjusted_sequence), 3)]
#         quality_sum = [sum(adjusted_quality[i:i + 3]) for i in range(0, len(adjusted_quality), 3)]
#
#         # TODO Trim regions beyond the CDS
#
#         # if max(codon_positions) > codonend // 3:
#         #     print(max(codon_positions))
#         #     codons = codons[:codonend // 3]
#         #     codon_positions = codon_positions[:codonend // 3]
#         #     quality_sum = quality_sum[:codonend // 3]
#
#         return codon_positions, codons, quality_sum
#
#     def filter_amino_acid(self, pairs):
#         key, value = pairs
#         if value[0] == "Unknown":
#             return False
#         elif value[1] < 100:
#             return False
#         else:
#             return True
#
#     def process_read(self, read, quality_threshold):
#
#     # Load CDS Regions
#         start = cds_regions[0][1] - 1
#         end = cds_regions[0][2] - 1
#
#         read_seq = Seq(read.query_sequence)
#         ref_positions = read.get_reference_positions()
#         quality = read.query_alignment_qualities
#
#         codon_positions, codons, quality_sum = read_to_codon(start, ref_positions, read_seq, quality,end)
#
#         # convert codons to amino acids
#         amino_acids = [codon_to_amino_acid.get(codon, 'Unknown') for codon in codons]
#
#         d = dict(zip(codon_positions, zip(amino_acids, quality_sum)))
#
#         ## filter out Ns and low quality mapping codons
#         d_clean = dict(filter(filter_amino_acid, d.items()))
#
#         #print(d_clean)
#         return d_clean
#
#     def process_sample(self):
#         with pysam.AlignmentFile(self.bam_file, "rb") as bam:
#             reads = list(bam.fetch())
#
#         count = 0
#         read_data = []
#         for read in reads:
#             if (read.mapping_quality > self.quality_threshold) and (read.reference_start > 0):
#                 result = self.process_read(read, self.quality_threshold)
#                 read_data.append(result)
#
#         compiled_dataset = defaultdict(lambda: defaultdict(int))
#         for dict_data in read_data:
#             try:
#                 for position, codon_data in dict_data.items():
#                     if position not in compiled_dataset.keys():
#                         compiled_dataset[position] = [codon_data]
#                     else:
#                         compiled_dataset[position].extend([codon_data])
#             except AttributeError:
#                 pass
#
#         ordered_dict = OrderedDict(sorted(compiled_dataset.items()))
#
#         codon_frequencies = {}
#
#         for position, codon_metadata_list in ordered_dict.items():
#             position_codon_counts = {}
#             for codon, metadata in codon_metadata_list:
#                 position_codon_counts[codon] = position_codon_counts.get(codon, 0) + 1
#             codon_frequencies[position] = position_codon_counts
#
#         return codon_frequencies
#
#     def call_amino_acid_variants(self, experimental_dict):
#         cds_start = self.cds_regions[0][1]
#         cds_end = self.cds_regions[0][2]
#         reference_protein = self.translate_cds(self.reference_sequence.seq, cds_start, cds_end)
#         reference_dict = {i + 1: reference_protein[i] for i in range(len(reference_protein))}
#         common_positions = set(reference_dict.keys()).intersection(experimental_dict.keys())
#
#         positions = []
#         mutant_freqs = []
#         wt_freqs = []
#
#         for position in common_positions:
#             reference_aa = reference_dict[position]
#             experimental_aa = experimental_dict[position]
#             WT = 0
#             MUT = 0
#             tmp_dict = dict()
#             for key, value in experimental_aa.items():
#                 if key == reference_aa:
#                     WT = value
#                 else:
#                     MUT += value
#             for key, value in experimental_aa.items():
#                 if key != reference_aa:
#                     tmp_dict[key] = value / (MUT + WT)
#             tmp_dict["WT"] = WT / (WT + MUT)
#
#             mutant_freqs.append(tmp_dict)
#             positions.append(position)
#
#         return dict(zip(positions, mutant_freqs))
#
#     def write_to_csv(self, frequency_dictionary, csv_file_path='../test_data/aa_frequencies.csv'):
#         all_aminoacids = set(self.codon_to_amino_acid.values())
#
#         with open(csv_file_path, 'w', newline='') as csvfile:
#             csv_writer = csv.writer(csvfile)
#
#             headers = ['Position'] + list(all_aminoacids)
#             csv_writer.writerow(headers)
#
#             for position, freq_data in frequency_dictionary.items():
#                 row_data = [position] + [freq_data.get(aminoacid, 0) for aminoacid in all_aminoacids]
#                 csv_writer.writerow(row_data)
#
#         print(f'CSV file saved at: {csv_file_path}')

# Example usage:
# Replace 'your_reference.fasta', 'your_cds.bed', 'your_bam.bam' with actual file paths
# and 'your_codon_to_amino_acid' with the actual codon_to_amino_acid dictionary
# analyzer = AminoAcidAnalyzer(codon_to_amino_acid, 'your_cds.bed', 'your_bam.bam', 'your_reference.fasta')
# frequency_dictionary = analyzer.process_sample()
# analyzer.write_to_csv(frequency_dictionary)



import pysam
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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


def load_cds_regions(bed_file):
    cds_regions = []
    with open(bed_file, 'r') as bed:
        for line in bed:
            chrom, start, end = line.strip().split('\t')
            cds_regions.append((chrom, int(start), int(end)))

    return cds_regions


def is_within_cds(position, cds_regions):
    for chrom, start, end in cds_regions:
        # print(position,start,end)
        if start <= position < end:
            return True
    return False


def read_to_codon(codonstart, refpositions, sequence, quality, codonend):

    #print(refpositions[0])
    if refpositions[0] < codonstart:
        adjust = codonstart-refpositions[0]
        refpositions = refpositions[adjust:]
        sequence = sequence[adjust:]
        quality = quality[adjust:]
        # print(refpositions)
        # print(sequence)
        # print(quality)

    if (refpositions[0] - codonstart) % 3 == 0:
        codonshift: int = 0
    elif (refpositions[0] - codonstart) % 3 == 1:
        codonshift: int = 2
    elif (refpositions[0] - codonstart) % 3 == 2:
        codonshift:int  = 1

    # adjust start position
    adjusted_positions = [pos + codonshift for pos in refpositions]

    # print(adjusted_positions)

    codon_positions = [1 + ((i - codonstart + codonshift) // 3) for i in adjusted_positions[0::3] if
                       i >= codonstart]

    # fix ends of the sequence

    if len(sequence[codonshift:]) % 3 == 0:
        endshift = 0
        adjusted_sequence = sequence[codonshift: len(sequence) - endshift]
        adjusted_quality = quality[codonshift: len(quality) - endshift]
        codon_positions = codon_positions[:]

    elif len(sequence[codonshift:]) % 3 == 1:
        endshift = 1
        adjusted_sequence = sequence[codonshift: len(sequence) - endshift]
        adjusted_quality = quality[codonshift: len(quality) - endshift]
        codon_positions = codon_positions[:-1]

    elif len(sequence[codonshift:]) % 3 == 2:
        endshift = 2
        adjusted_sequence = sequence[codonshift: len(sequence) - endshift]
        adjusted_quality = quality[codonshift: len(quality) - endshift]
        codon_positions = codon_positions[:-2]

    # print(sequence)
    # print(adjusted_sequence)

    codons = [adjusted_sequence[i:i + 3] for i in range(0, len(adjusted_sequence), 3)]
    quality_sum = [sum(adjusted_quality[i:i + 3]) for i in range(0, len(adjusted_quality), 3)]

    # TODO Trim regions beyond the CDS

    # if max(codon_positions) > codonend // 3:
    #     print(max(codon_positions))
    #     codons = codons[:codonend // 3]
    #     codon_positions = codon_positions[:codonend // 3]
    #     quality_sum = quality_sum[:codonend // 3]

    return codon_positions, codons, quality_sum


def filter_amino_acid(pairs):
    key,value = pairs
    #value in codon
    if value[0] == "Unknown":
        return False
    # codon sum score
    elif value[1] < 100:
        return False
    else:
        return True

def process_read(read, quality_threshold, cds_regions):

    # Load CDS Regions
    start = cds_regions[0][1] - 1
    end = cds_regions[0][2] - 1

    read_seq = Seq(read.query_sequence)
    ref_positions = read.get_reference_positions()
    quality = read.query_alignment_qualities

    codon_positions, codons, quality_sum = read_to_codon(start, ref_positions, read_seq, quality,end)

    # convert codons to amino acids
    amino_acids = [codon_to_amino_acid.get(codon, 'Unknown') for codon in codons]

    d = dict(zip(codon_positions, zip(amino_acids, quality_sum)))

    ## filter out Ns and low quality mapping codons
    d_clean = dict(filter(filter_amino_acid, d.items()))

    #print(d_clean)
    return d_clean


def process_sample(bam_file, bed_file, quality_threshold=20):

    # Load reads
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        reads = list(bam.fetch())

    # Load bed regions
    cds_regions = load_cds_regions(bed_file)

    count = 0
    read_data = []
    for read in reads:
        if (read.mapping_quality > quality_threshold) and (read.reference_start > 0):
            result = process_read(read, quality_threshold, cds_regions)
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
    codon_frequencies = {}

    # Iterate through each position in the mapped_positions dictionary
    for position, codon_metadata_list in ordered_dict.items():
        # Dictionary to store codon frequencies for the current position
        position_codon_counts = {}

        # Iterate through each codon-metadata tuple in the list
        for codon, metadata in codon_metadata_list:
            # Count the frequency of each codon for the current position
            position_codon_counts[codon] = position_codon_counts.get(codon, 0) + 1

        # Store the codon frequencies for the current position in the result dictionary
        codon_frequencies[position] = position_codon_counts

    # Now, codon_frequencies is a dictionary where the key is the position,
    # and the value is another dictionary containing codon frequencies
    print(codon_frequencies)
    # df = pd.DataFrame.from_dict(codon_frequencies, orient='index')
    # print(df)

    # Get a set of all unique codon types
    all_codons = set(codon for freq_data in codon_frequencies.values() for codon in freq_data)

    # Create a CSV file and write headers
    csv_file_path = '../test_data/codon_frequencies.csv'
    with open(csv_file_path, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)

        # Write headers (position and all codon types)
        headers = ['Position'] + list(all_codons)
        csv_writer.writerow(headers)

        # Write frequency data for each position
        for position, freq_data in codon_frequencies.items():
            row_data = [position] + [freq_data.get(codon, 0) for codon in all_codons]
            csv_writer.writerow(row_data)

    print(f'CSV file saved at: {csv_file_path}')

    return codon_frequencies

def read_fasta(fasta_file):
    # Read in the FASTA reference sequence
    record = SeqIO.read(fasta_file, "fasta")
    return record

def translate_cds(sequence, cds_start, cds_end):
    # Extract the coding sequence and translate it into amino acids
    cds_sequence = sequence[cds_start - 1:cds_end]  # Assuming 1-based indexing
    amino_acids = cds_sequence.translate()
    return amino_acids


def call_amino_acid_variants(experimental_dict,reference_sequence,cds_regions):

    cds_start = cds_regions[0][1]
    cds_end = cds_regions[0][2]
    reference_protein = translate_cds(reference_sequence,cds_start,cds_end)
    reference_dict = {i + 1: reference_protein[i] for i in range(len(reference_protein))}
    common_positions = set(reference_dict.keys()).intersection(experimental_dict.keys())

    positions = []
    mutant_freqs = []
    wt_freqs = []
    # Compare amino acids at each position
    for position in common_positions:
        reference_aa = reference_dict[position]
        experimental_aa = experimental_dict[position]
        WT = 0
        MUT = 0
        tmp_dict = dict()
        for key, value in experimental_aa.items():
            if key == reference_aa:
                WT = value
            else:
                MUT += value
        for key, value in experimental_aa.items():
            if key != reference_aa:
                tmp_dict[key] = value / (MUT + WT)
        tmp_dict["WT"] = WT / (WT + MUT)



        mutant_freqs.append(tmp_dict)
        positions.append(position)


    return dict(zip(positions,mutant_freqs))




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
