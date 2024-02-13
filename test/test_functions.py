import unittest
from CodonCaller.variant_caller import *
import pysam
from collections import Counter


class TestRead:
    def __init__(self, sequence, position_start, quality):
        self.read = pysam.AlignedSegment()
        self.positions = list(range(position_start + len(sequence)))
        # self.read.reference_start = position_start
        # self.read.reference_end = self.positions[-1]
        self.read.query_sequence = str(sequence)
        self.read.query_qualities = len(sequence) * [quality]


# class TestMockRead(unittest.TestCase):
#     def setUp(self):
#         bed_file = "../test_data/bed_file/REJOc_multiple.bed"
#         bam_file = "../test_data/mapped_reads/CD00-m000-0-REJOc.bam"
#         reference_fasta = "../test_data/reference_sequence/REJOc-reference.fa"
#         read_quality_threshold = 40
#         base_quality_threshold = 30
#         self.test_instance = VariantCaller(bed_file, bam_file, reference_fasta, read_quality_threshold,
#                                            base_quality_threshold)
#
#     def test_load_cds_regions(self):
#         self.test_instance.load_cds_regions()
#         self.assertEqual(self.test_instance.cds_regions['Env'][0], 52)
#         self.assertEqual(self.test_instance.cds_regions['Env'][1], 2619)


class TestYourCode(unittest.TestCase):

    def setUp(self):
        bed_file = "../test_data/bed_file/REJOc_multiple.bed"
        bam_file = "../test_data/mapped_reads/CD5-m265-25-REJOc.bam"
        reference_fasta = "../test_data/reference_sequence/REJOc-reference.fa"
        read_quality_threshold = 40
        base_quality_threshold = 30
        self.test_instance = VariantCaller(bed_file, bam_file, reference_fasta, read_quality_threshold,
                                           base_quality_threshold)

    def test_load_cds_regions(self):
        self.test_instance.load_cds_regions()
        self.assertEqual(self.test_instance.cds_regions['Env'][0], 52)
        self.assertEqual(self.test_instance.cds_regions['Env'][1], 2619)

    def test_read_fasta(self):
        self.test_instance.read_fasta()
        self.assertEqual(len(self.test_instance.reference_sequence), 2676)

    def test_translate_cds(self):
        self.test_instance.translate_cds()
        print(self.test_instance.cds_aminoacids)

    def test_read_to_codon(self):
        cds = "Env"
        read_sequence = "ATGAAAGTGAAGGGGATCAGGAGGAATTAT"
        read_positions = list(range(53, 53 + len(read_sequence)))
        quality = len(read_sequence) * [40]
        indels = {"insertion": [(3, 3)], "deletion": [(6, 3)]}
        codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
                                                                                quality, indels)
        print(codon_positions)
        print(codons)
        print(quality_sum)

    def test_read_to_codon_frame_1(self):
        cds = "Env"
        read_sequence = "TGAAAGTGAAGGGGATCAGGAGGAATTATC"
        read_positions = list(range(54, 54 + len(read_sequence)))
        quality = len(read_sequence) * [40]
        codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
                                                                                quality)
        print(codon_positions)
        print(codons)
        print(quality_sum)

    def test_read_to_codon_frame_2(self):
        cds = "Env"
        read_sequence = "GAAAGTGAAGGGGATCAGGAGGAATTATCA"
        read_positions = list(range(55, 55 + len(read_sequence)))
        quality = len(read_sequence) * [40]
        codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
                                                                                quality)
        print(codon_positions)
        print(codons)
        print(quality_sum)

    def test_read_to_codon_real_read(self):
        cds = "Env"
        read = self.test_instance.reads[5]
        read_sequence, read_positions, quality, indels = parse_read(read)
        codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
                                                                                quality, indels)
        print(codon_positions)
        print(codons)
        # print(quality_sum)

    def test_fix_insert(self):
        read = self.test_instance.reads[6]
        cds = "Env"
        read_sequence, read_positions, quality, indels = parse_read(read)
        insertion = indels["insertion"][0]
        out = fix_insert(read_positions, insertion)
        print(out)

    def test_fix_deletion(self):
        read = self.test_instance.reads[6]
        cds = "Env"
        read_sequence, read_positions, quality, indels = parse_read(read)
        deletion = indels["deletion"][0]
        seq, pos = fix_deletion(read_sequence, read_positions, deletion)
        print(deletion)
        print(read_positions)
        print(read_sequence)
        print(seq)
        print(pos)

    def test_fix_indels(self):
        read = self.test_instance.reads[6]
        cds = "Env"
        read_sequence, read_positions, quality, indels = parse_read(read)
        print(indels)
        print(read_positions)
        print(quality)
        print(read_sequence)

        read_sequence, read_positions, quality = fix_indels(read_sequence, read_positions, quality, indels)

        print(read_positions)
        print(read_sequence)
        print(quality)

    def test_process_read(self):
        read = self.test_instance.reads[6]
        cds = "Env"
        read_sequence, read_positions, quality, indels = parse_read(read)
        print(indels)
        print(read_positions)
        print(read_sequence)
        codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
                                                                                quality, indels)
        print(codon_positions)
        print(codons)
        print(quality_sum)

    def test_process_read2(self):
        read = self.test_instance.reads[6]
        cds = "Env"

        read_sequence, read_positions, quality, indels = parse_read(read)

        codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
                                                                                quality, indels)
        # print(read_sequence)
        # print(read_positions)
        # print(codon_positions)
        # print(codons)
        # print(quality_sum)

        dict_out = self.test_instance.process_read(cds, read)
        print(dict_out)

    def test_process_sample(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        print(self.test_instance.mut_call)

    def test_process_sample(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        print(self.test_instance.aa_counts)


    def test_compiled_reads(self):
        cds = "Env"
        out_dict = self.test_instance.compile_reads(cds)
        print(len(out_dict['77i']) + len(out_dict['77']))
        clean_dict = self.test_instance.fix_insertions(out_dict)
        print(len(clean_dict['77']))
        # print(out_dict.values())


    def test_write_data(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        self.test_instance.aa_count_to_freq(cds)
        self.test_instance.write_to_csv("../test_data/calls/test.csv")

    def test_coverage(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        self.test_instance.aa_count_to_freq(cds)
        print(self.test_instance.aa_coverage)
        self.test_instance.write_coverage_to_csv("../test_data/calls/CD20-m999-8-REJOc-coverage.csv")

    def test_INDEL(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        self.test_instance.aa_count_to_freq(cds)
        print(self.test_instance.aa_coverage)
        # self.test_instance.aa_count_to_freq(cds)

    def test_fix_insert(self):
        read = self.test_instance.reads[5]
        read_sequence, read_positions, quality, indels = parse_read(read)
        sorted_indel = sort_indel_dict(indels)
        print(sorted_indel)
        corrected_positions = read_positions
        for position in sorted_indel.keys():
            if sorted_indel[position][0] == "insertion":
                print(position)
                start = position
                length = sorted_indel[position][1]
                corrected_positions = fix_insert(corrected_positions, start, length)
        print(read_positions)
        count = Counter(corrected_positions)
        print(count)
        assert (count[58] == 2) & (count[142] == 2)

    def test_fix_indels(self):
        read = self.test_instance.reads[6]
        read_sequence, read_positions, quality, indels = parse_read(read)
        print(indels)
        print(read_sequence)
        print(read_positions)
        # print(quality)
        read_sequence, read_positions, quality = fix_indels(read_sequence, read_positions, quality, indels)
        print(read_sequence)
        print(read_positions)
        # print(quality)


    def test_fix_deletion(self):

        read = self.test_instance.reads[5]
        read_sequence, read_positions, quality, indels = parse_read(read)
        sorted_indel = sort_indel_dict(indels)
        print(sorted_indel)
        corrected_positions = read_positions
        total_insert_shift = 3
        for position in sorted_indel.keys():
            if sorted_indel[position][0] == "deletion":
                # print('del:',position)
                start = position
                length = sorted_indel[position][1]
                shift = total_insert_shift
                read_sequence, read_positions, quality = fix_deletion(read_sequence, read_positions, quality, start, length, shift)
        print(read_positions)

#     def test_translate_csf(self):
#         sequence = self.analyzer.read_fasta().seq
#         cds_regions = self.analyzer.load_cds_regions()
#         cds_start, cds_end = cds_regions[0][1], cds_regions[0][2]
#
#         out = self.analyzer.translate_cds(sequence, cds_start, cds_end)
#         self.assertIsInstance(out, Seq())
#
#     def test_is_within_cds(self):
#         position_within_cds = 150
#         position_outside_cds = 30
#         self.assertTrue(self.analyzer.is_within_cds(position_within_cds))
#         self.assertFalse(self.analyzer.is_within_cds(position_outside_cds))
#
#     def test_read_frame_0(self):
#         start_position = self.analyzer.cds_regions[0][1]
#
#         # Simulate Aligned Segments
#         sequence = "ATGAAAGTGAAGGGGATCAGGAGGAATTAT"
#         positions = list(range(53, 53 + len(sequence)))
#         quality = len(sequence) * [25]
#
#         codon_positions, codons, quality_sum = self.analyzer.read_to_codon(start_position, positions, sequence, quality, 0)
#         self.assertEqual(codon_positions, list(range(1, 10)))
#
#     def test_read_frame_1(self):
#         start_position = self.analyzer.cds_regions[0][1]
#
#         # Simulate Aligned Segments
#         sequence = "TGAAAGTGAAGGGGATCAGGAGGAATTATC"
#         positions = list(range(54, 54 + len(sequence)))
#         quality = len(sequence) * [25]
#
#         codon_positions, codons, quality_sum = self.analyzer.read_to_codon(start_position, positions, sequence, quality, 0)
#         self.assertEqual(codon_positions, list(range(2, 10)))
#
#     def test_read_frame_2(self):
#         start_position = self.analyzer.cds_regions[0][1]
#
#         # Simulate Aligned Segments
#         sequence = "GAAAGTGAAGGGGATCAGGAGGAATTATCA"
#         positions = list(range(55, 55 + len(sequence)))
#         quality = len(sequence) * [25]
#
#         codon_positions, codons, quality_sum = self.analyzer.read_to_codon(start_position, positions, sequence, quality, 0)
#         self.assertEqual(codon_positions, list(range(3, 10)))
#
#     def test_read_non_coding(self):
#         start_position = self.analyzer.cds_regions[0][1]
#
#         # Simulate Aligned Segments
#         sequence = "AAAATAGACAGGTTAATTGATAGAATAAGAGATAGAGCAGAAGACAGTGGCAATGAAAGTGAAGGGGATC"
#         positions = list(range(1, 1 + len(sequence)))
#         quality = len(sequence) * [35]
#
#         codon_positions, codons, quality_sum = self.analyzer.read_to_codon(start_position, positions, sequence, quality, 0)
#         self.assertEqual(codon_positions, list(range(1, 10)))
#
#     def test_process_read(self):
#         # Create test data or mock objects as needed
#         read_object = your_mock_read_object
#         quality_threshold = 20
#         result = self.analyzer.process_read(read_object, quality_threshold)
#         self.assertIsInstance(result, dict)
#
#     def test_process_sample(self):
#         result = self.analyzer.process_sample()
#         self.assertIsInstance(result, dict)
#
#     def test_process_sample_call(self):
#         result = self.analyzer.process_sample()
#         out = self.analyzer.call_amino_acid_variants(result)
#         self.assertIsInstance(out, dict)
#
#     def test_write_to_csv(self):
#         # Create test data or mock objects as needed
#         result = self.analyzer.process_sample()
#         out = self.analyzer.call_amino_acid_variants(result)
#         csv_file_path = '../test_data/aa_frequencies.csv'
#         self.analyzer.write_to_csv(frequency_dictionary, csv_file_path)
#         # Add assertions based on the expected output
#
# if __name__ == '__main__':
#     unittest.main()
