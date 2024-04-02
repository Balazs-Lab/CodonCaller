import pysam
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict, OrderedDict

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
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'NNN': 'DEL'
}


def parse_read(read):
    """

    :param read: perform operation on a specific read
    :return:component sequence, position, quality, and a dict of INDEL positions
    """

    read_sequence = Seq(read.query_sequence)
    read_positions = read.get_reference_positions()
    quality = read.query_alignment_qualities
    indels = {"insertion": [], "deletion": []}
    position = list(read_positions)[0]

    for i in range(len(read.cigartuples)):
        item = read.cigartuples[i]
        if (item[0] != 1) & (item[0] != 2) & (item[0] != 3) & (item[0] != 6):
            # these are the mapped reads
            position += item[1]
        if item[0] == 1:
            indels["insertion"].append((position, item[1]))
            # position += item[1]
        elif (item[0] == 2) or (item[0] == 3):
            indels["deletion"].append((position, item[1]))
            position += item[1]

    return read_sequence, read_positions, quality, indels


def fix_insertions(ordered_dict):
    """

    :param ordered_dict: a dictionary where the key is a base and the value is the amino acid
    :return: integrates insertions in the dictionary by annotating the amino acid to indicate insertion
    then appending those values to the list of amino acids at the site of insertion
    """
    deletion_list = []
    new_ordered_dict = ordered_dict
    for key, value in ordered_dict.items():
        if "i" in key:
            # annotate mutations to be insertions
            insertions = []
            for item in value:
                insertions.append((item[0] + "i", item[1]))

            location = key[:-1]
            if location in new_ordered_dict.keys():
                new_ordered_dict[location].extend(insertions)
            else:
                new_ordered_dict[location] = insertions
            # ordered_dict[location].append(insertions)
            deletion_list.append(key)

    for key in deletion_list:
        del new_ordered_dict[key]
    return new_ordered_dict


class VariantCaller:
    def __init__(self, bed_file, bam_file, reference_fasta, read_quality_threshold=20, base_quality_threshold=35):
        self.aa_coverage = {}
        self.mut_call = None
        self.codon_frequencies = None
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

    def map_correct(self, read):

        # # adjust mis-mapped bases in deletions
        reference = str(read.reference_name)
        start = self.cds_regions[reference][0]
        end = self.cds_regions[reference][1]
        codon_starts = [i for i in range(start, end, 3)]
        codon_one_shift = [i - 1 for i in range(start, end, 3)]
        codon_two_shift = [i - 2 for i in range(start, end, 3)]

        read_sequence, read_positions, quality, indels = parse_read(read)

        corrected_deletions = list()
        for deletion in indels['deletion']:
            if deletion[1] % 3 != 0:
                corrected_deletions.append(deletion)

            elif deletion[0] in codon_starts:
                corrected_deletions.append(deletion)

            elif deletion[0] in codon_one_shift:
                index = read_positions.index(deletion[0]+deletion[1])
                read_positions[index] = deletion[0]
                corrected_deletions.append((deletion[0]+1, deletion[1]))

            elif deletion[0] in codon_two_shift:
                index = read_positions.index(deletion[0]+deletion[1])
                read_positions[index] = deletion[0]
                corrected_deletions.append((deletion[0]+2, deletion[1]))

        indels['deletion'] = corrected_deletions

        return read_sequence, read_positions, quality, indels



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
        # amino acid value - codon with an "N" maps to Unknown except for NNN which is a DEL
        # I still do not like how this works - should be updated to better handel DELs
        if value[0] == "Unknown":
            return False
        # codon sum score - sum of all 3 nucleotide mapping scores
        elif value[1] < self.read_quality_threshold:
            return False
        else:
            return True

    def read_to_codon(self, cds, read_sequence, read_positions, quality, indels):
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

        ## Fix Indels
        read_sequence, read_positions, quality = fix_indels(read_sequence, read_positions, quality, indels)
        # print(read_sequence)
        # print(read_positions)

        # adjust read position based on CDS
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

        codon_positions = [str(1 + ((i - cds_start + codonshift) // 3)) for i in adjusted_positions[0::3] if
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

        # read_sequence, read_positions, quality, indels = parse_read(read)
        read_sequence, read_positions, quality, indels = self.map_correct(read)

        ## TODO turn this into a proper filter
        if len(read_positions) > 100:
            codon_positions, codons, quality_sum = self.read_to_codon(cds, read_sequence, read_positions, quality,
                                                                      indels)
            # print(len(codon_positions), codon_positions)
            # print(len(codons), codons)
            # print(len(quality_sum), quality_sum)

            # convert codons to amino acids
            amino_acids = [codon_to_amino_acid.get(codon, 'Unknown') for codon in codons]
            # print(amino_acids)

            # annotate insertions to new location
            for i, j in enumerate(codon_positions[:-1]):
                if j == codon_positions[i + 1]:
                    codon_positions[i] = str(codon_positions[i]) + 'i'
            # print(codon_positions)
            d = dict(zip(codon_positions, zip(amino_acids, quality_sum)))

            ## filter out Ns and low quality mapping codons
            d_clean = dict(filter(self.filter_amino_acid, d.items()))

            # print(d_clean)
            return d_clean
        else:
            pass

    def compile_reads(self, cds):
        """

        :param cds:
        :return: dictionary where the key is each ammino acid position (including insertions) and the value is a list
        tuples with the amino acid and summed quality score
        """
        count = 0
        read_data = []
        for read in self.reads:
            if (read.mapping_quality > self.read_quality_threshold) and (read.reference_start >= 0):
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

        ordered_dict = compiled_dataset  # OrderedDict(sorted(compiled_dataset.items())) might need to reverse this
        return ordered_dict

    def process_sample(self, cds):
        ordered_dict = self.compile_reads(cds)
        # move insertions to original base
        ordered_dict = fix_insertions(ordered_dict)

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
            self.aa_coverage[position] = sum(position_codon_counts.values())

        # Now, codon_frequencies is a dictionary where the key is the position,
        # and the value is another dictionary containing codon frequencies
        # print(self.codon_frequencies)

    def aa_count_to_freq(self, cds):

        # cds_start = self.cds_regions[cds][0]
        # cds_end = self.cds_regions[cds][1]
        # cds_len = (cds_end - cds_start)//3

        reference_protein = self.cds_aminoacids[cds]

        reference_dict = {str(i + 1): reference_protein[i] for i in range(len(reference_protein))}
        # print(reference_dict)

        # iterate through all positions from mapped data (including insertions)
        # common_positions = list(self.aa_counts.keys())
        # common_positions.sort()
        # print(common_positions)
        common_positions = set(reference_dict.keys()).intersection(self.aa_counts.keys())
        # print(common_positions)

        mutant_call_list = []
        positions = []
        mutant_freq = []
        # print(reference_dict)
        # print(common_positions)

        # Compare amino acids at each position

        for position in common_positions:

            reference_aa = reference_dict[position]
            experimental_aa = self.aa_counts[position]
            # print(reference_aa,experimental_aa)

            wt = 0
            mut = 0
            tmp_dict = dict()
            wt_id = reference_aa
            for key, value in experimental_aa.items():
                if key == reference_aa:
                    wt = value
                else:
                    mut += value
            if mut > 0:
                for key, value in experimental_aa.items():
                    if key != reference_aa:
                        tmp_dict[key] = value / (mut + wt)
                        tmp_out = [position, wt_id, key, value / (mut + wt)]
                        mutant_call_list.append(tmp_out)
                tmp_dict["WT"] = wt / (wt + mut)

            mutant_freq.append(tmp_dict)
            positions.append(position)
        # print(mutant_call_list)
        mut_call = pd.DataFrame(mutant_call_list, columns=['POS_AA', 'REF_AA', "ALT_AA", "ALT_FREQ"])
        self.mut_call = mut_call
        # print(self.mut_call)
        # self.aa_freq = dict(zip(positions, mutant_freq))

    def write_to_csv(self, csv_file_path='../test_data/aa_frequencies.csv'):
        self.mut_call.to_csv(csv_file_path)

    def write_coverage_to_csv(self, csv_file_path='../test_data/aa_coverage.csv'):
        coverage_df = pd.DataFrame.from_dict(self.aa_coverage, orient='index')
        coverage_df = coverage_df.reset_index(0)
        coverage_df.columns = ["POS_AA", "COVERAGE"]
        coverage_df.to_csv(csv_file_path)


def sort_indel_dict(indels):
    sorted_indel = dict()
    for k, v in indels.items():
        type = k
        for item in v:
            start = item[0]
            length = item[1]
            sorted_indel[start] = (type, length)
    sorted_indel = dict(OrderedDict(sorted(sorted_indel.items())))

    return sorted_indel


def fix_indels(read_sequence, read_positions, quality, indels):
    # order the indels

    sorted_indel = sort_indel_dict(indels)

    # need to shift the deletion position for every insertion called
    del_shift = 0
    # insert counter
    total_insert_shift = 0
    # iterate though positions
    # print(indels)
    # print('pre-fix',len(read_positions),len(read_sequence))
    for position in sorted_indel.keys():

        if sorted_indel[position][0] == "insertion":
            start = position
            length = sorted_indel[position][1]
            read_positions = fix_insert(read_positions, start, length, del_shift + total_insert_shift)
            total_insert_shift += length

        elif sorted_indel[position][0] == "deletion":
            start = position
            length = sorted_indel[position][1]
            shift = total_insert_shift
            del_shift += length
            read_sequence, read_positions, quality = fix_deletion(read_sequence, read_positions, quality, start, length,
                                                                  shift)

    # print('adjusted positions',len(read_positions),len(read_sequence))
    return read_sequence, read_positions, quality


def fix_insert(read_positions, start, length, del_shift):
    """
    :param length:
    :param start:
    :param read_positions: a list of positions corresponding to a read sequence
    :param insertion: an insertion tuple with the format (start_position, insertion_length)
    :return: a corrected set of read positions where the inserted positions have a character "a" appended to them
    and the positions following the insertion are adjusted based on the size of the insert to map to their
    original positions
    """

    full_insertion = list(range(start, start + length))
    # corrected_read_positions = []

    corrected_read_positions = read_positions[:start - min(read_positions) + del_shift]
    corrected_read_positions.extend(full_insertion)
    corrected_read_positions.extend(read_positions[start - min(read_positions) + del_shift:])

    return corrected_read_positions


def fix_deletion(read_sequence, read_positions, quality, start, length, shift):
    # full_deletion = list(range(start - shift, start - shift + length))
    full_deletion = list(range(start, start + length))
    print(full_deletion)

    # update sequence
    corrected_read_sequence = read_sequence[:start - min(read_positions) + shift]
    corrected_read_sequence = corrected_read_sequence + length * "N"
    corrected_read_sequence = corrected_read_sequence + read_sequence[start - min(read_positions) + shift:]

    # update positions
    # print(len(read_positions))
    corrected_read_positions = read_positions[: start - min(read_positions) + shift]
    corrected_read_positions.extend(full_deletion)
    corrected_read_positions.extend(read_positions[start - min(read_positions) + shift:])
    # print(len(corrected_read_positions))

    # update quality scores
    corrected_quality = quality[: start - min(read_positions) + shift]
    corrected_quality.extend(length * [99])
    corrected_quality.extend(quality[start - min(read_positions) + shift:])

    return corrected_read_sequence, corrected_read_positions, corrected_quality
