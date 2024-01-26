import unittest
from CodonCaller.variant_caller import VariantCaller
import pysam


class TestRead:
    def __init__(self, sequence, position_start, quality):
        self.read = pysam.AlignedSegment()
        self.positions = list(range(position_start + len(sequence)))
        #self.read.reference_start = position_start
        # self.read.reference_end = self.positions[-1]
        self.read.query_sequence = str(sequence)
        self.read.query_qualities = len(sequence) * [quality]



class TestYourCode(unittest.TestCase):

    def setUp(self):
        bed_file = "../test_data/bed_file/REJOc_multiple.bed"
        bam_file = "../test_data/mapped_reads/CD5-m246-20-REJOc.bam"
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
        self.assertEqual(len(self.test_instance.reference_sequence),2676)

    def test_translate_cds(self):
        self.test_instance.translate_cds()
        print(self.test_instance.cds_aminoacids)

    def test_read_to_codon(self):
        cds = "Env"
        read_sequence = "ATGAAAGTGAAGGGGATCAGGAGGAATTAT"
        read_positions = list(range(53, 53 + len(read_sequence)))
        quality = len(read_sequence) * [40]
        codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions, quality)
        print(codon_positions)
        print(codons)
        print(quality_sum)

    def test_read_to_codon_frame_1(self):
        cds = "Env"
        read_sequence = "TGAAAGTGAAGGGGATCAGGAGGAATTATC"
        read_positions = list(range(54, 54 + len(read_sequence)))
        quality = len(read_sequence) * [40]
        codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions, quality)
        print(codon_positions)
        print(codons)
        print(quality_sum)

    def test_read_to_codon_frame_2(self):
        cds = "Env"
        read_sequence = "GAAAGTGAAGGGGATCAGGAGGAATTATCA"
        read_positions = list(range(55, 55 + len(read_sequence)))
        quality = len(read_sequence) * [40]
        codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions, quality)
        print(codon_positions)
        print(codons)
        print(quality_sum)


    def test_process_read(self):
        read = self.test_instance.reads[300]
        cds = "Env"
        read_sequence, read_positions, quality = self.test_instance.parse_read(read)
        codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions, quality)
        print(codon_positions)
        print(codons)
        print(quality_sum)


    def test_process_sample(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        print(self.test_instance.aa_counts)

    def test_process_sample(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        self.test_instance.aa_count_to_freq(cds)
        print(self.test_instance.aa_freq)


    def test_write_data(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        self.test_instance.aa_count_to_freq(cds)
        self.test_instance.write_to_csv("../test_data/calls/CD5-m246-20-REJOc.csv")

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


# import unittest
# from codon_call import *
# import pysam
# from Bio import SeqIO
# import csv
#
#
# class TestYourCode(unittest.TestCase):
#
#
#
#     def test_load_cds_regions(self):
#         bed_file = "../test_data/bed_file/REJOc.bed"
#         cds_regions = load_cds_regions(bed_file)
#         print(cds_regions)
#         self.assertEqual(len(cds_regions), 1)  # Adjust based on the content of your test BED file
#
#     def test_translate_csf(self):
#         reference_sequence = "../test_data/reference_sequence/REJOc-reference.fa"
#         bed_file = "../test_data/bed_file/REJOc.bed"
#         cds_regions = load_cds_regions(bed_file)
#
#         sequence = read_fasta(reference_sequence)
#         cds_regions = load_cds_regions(bed_file)
#         cds_start = cds_regions[0][1]
#         cds_end = cds_regions[0][2]
#
#         out = translate_cds(sequence, cds_start, cds_end)
#         print(str(out))
#         # self.assertEqual()
#
#     def test_is_within_cds(self):
#         position_within_cds = 150
#         position_outside_cds = 30
#         bed_file = "../test_data/bed_file/REJOc.bed"
#         cds_regions = load_cds_regions(bed_file)
#         self.assertTrue(is_within_cds(position_within_cds, cds_regions))
#         self.assertFalse(is_within_cds(position_outside_cds, cds_regions))
#
#     def test_read_frame_0(self):
#         bed_file = "../test_data/bed_file/REJOc.bed"
#         cds_regions = load_cds_regions(bed_file)
#         start_position = cds_regions[0][1]
#         print(start_position)
#
#         # Simulate Aligned Segments
#         mapping_quality = 25
#
#         ## in frame
#         sequence = "ATGAAAGTGAAGGGGATCAGGAGGAATTAT"
#         positions = list(range(53, 53 + len(sequence)))
#         qualty = len(sequence) * [25]
#
#         print(sequence)
#         print(positions)
#
#         codon_positions, codons, qualty_sum = read_to_codon(start_position,positions,sequence,qualty)
#         print(codon_positions)
#         print(codons)
#         print(qualty_sum)
#
#         self.assertTrue(codon_positions, list(range(1, 10)))
#
#     def test_read_frame_1(self):
#         bed_file = "../test_data/bed_file/REJOc.bed"
#         cds_regions = load_cds_regions(bed_file)
#         start_position = cds_regions[0][1]
#         print(start_position)
#
#         # Simulate Aligned Segments
#         mapping_quality = 25
#
#         ## in frame
#         sequence = "TGAAAGTGAAGGGGATCAGGAGGAATTATC"
#         positions = list(range(54, 54 + len(sequence)))
#         qualty = len(sequence) * [25]
#
#         print(sequence)
#         print(positions)
#
#         codon_positions, codons, qualty_sum = read_to_codon(start_position,positions,sequence,qualty)
#         print(codon_positions)
#         print(codons)
#         print(qualty_sum)
#
#         self.assertTrue(codon_positions, list(range(2, 10)))
#
#     def test_read_frame_2(self):
#         bed_file = "../test_data/bed_file/REJOc.bed"
#         cds_regions = load_cds_regions(bed_file)
#         start_position = cds_regions[0][1]
#         print(start_position)
#
#         # Simulate Aligned Segments
#         mapping_quality = 25
#
#         ## in frame 2
#         sequence = "GAAAGTGAAGGGGATCAGGAGGAATTATCA"
#         positions = list(range(55, 55 + len(sequence)))
#         qualty = len(sequence) * [25]
#
#         print(sequence)
#         print(positions)
#
#         codon_positions, codons, qualty_sum = read_to_codon(start_position,positions,sequence,qualty)
#         print(codon_positions)
#         print(codons)
#         print(qualty_sum)
#
#         self.assertTrue(codon_positions, list(range(2, 10)))
#
#     def test_read_non_coding(self):
#         bed_file = "../test_data/bed_file/REJOc.bed"
#         cds_regions = load_cds_regions(bed_file)
#         start_position = cds_regions[0][1]
#         print(start_position)
#
#         # Simulate Aligned Segments
#         mapping_quality = 25
#
#         ## in frame
#         sequence = "AAAATAGACAGGTTAATTGATAGAATAAGAGATAGAGCAGAAGACAGTGGCAATGAAAGTGAAGGGGATC"
#         positions = list(range(1, 1 + len(sequence)))
#         qualty = len(sequence) * [35]
#
#         print(sequence)
#         print(positions)
#
#         codon_positions, codons, quality_sum = read_to_codon(start_position,positions,sequence,qualty)
#         d = dict(zip(codon_positions, zip(codons, quality_sum)))
#         print(d)
#         # print(codons)
#         # print(qualty_sum)
#
#         self.assertTrue(codon_positions, list(range(1, 10)))
#
#     def test_process_read(self):
#         # Create test data or mock objects as needed
#         bam_file = "../test_data/mapped_reads/CD5-m246-20-REJOc.bam"
#         with pysam.AlignmentFile(bam_file, "rb") as bam:
#             reads = list(bam.fetch())
#         read_object = reads[123]
#
#         # Load reference sequence
#         reference_seq = SeqIO.read("../test_data/reference_index/REJOc-reference.fa", "fasta").seq
#
#         quality_threshold = 20
#
#         # Load bed regions
#         bed_file = "../test_data/bed_file/REJOc.bed"
#         cds_regions = load_cds_regions(bed_file)
#
#         result = process_read(read_object, quality_threshold, cds_regions)
#         self.assertIsInstance(result, dict)
#
#     def test_process_sample(self):
#         # Create test data or mock objects as needed
#         bam_file = "../test_data/mapped_reads/CD5-m246-20-REJOc.bam"
#
#         quality_threshold = 40
#
#         # Load bed regions
#         bed_file = "../test_data/bed_file/REJOc.bed"
#         cds_regions = load_cds_regions(bed_file)
#
#         result = process_sample(bam_file, bed_file, quality_threshold)
#         self.assertIsInstance(result, dict)
#
#     def test_process_sample_JRCSF(self):
#         # Create test data or mock objects as needed
#         bam_file = "../test_data/mapped_reads/CD21-m009-12-JRCSF.bam"
#
#         quality_threshold = 40
#
#         # Load bed regions
#         bed_file = "../test_data/bed_file/JRCSF.bed"
#         cds_regions = load_cds_regions(bed_file)
#
#         result = process_sample(bam_file, bed_file, quality_threshold)
#         self.assertIsInstance(result, dict)
#
#
#     def test_process_sample_call(self):
#
#         # Create test data or mock objects as needed
#         bam_file = "../test_data/mapped_reads/CD5-m246-20-REJOc.bam"
#         reference_sequence = "../test_data/reference_sequence/REJOc-reference.fa"
#         sequence = read_fasta(reference_sequence)
#         bed_file = "../test_data/bed_file/REJOc.bed"
#         quality_threshold = 40
#
#         cds_regions = load_cds_regions(bed_file)
#
#         result = process_sample(bam_file, bed_file, quality_threshold)
#
#         out = call_amino_acid_variants(result, sequence, cds_regions)
#         print(out)
#         self.assertIsInstance(result, dict)
#
#     def test_process_sample_call2(self):
#
#         # Create test data or mock objects as needed
#         bam_file = "../test_data/mapped_reads/CD21-m009-12-JRCSF.bam"
#         reference_sequence = "../test_data/reference_sequence/JRCSF-reference.fa"
#         sequence = read_fasta(reference_sequence)
#         bed_file = "../test_data/bed_file/JRCSF.bed"
#         quality_threshold = 40
#
#         cds_regions = load_cds_regions(bed_file)
#
#         result = process_sample(bam_file, bed_file, quality_threshold)
#
#         out = call_amino_acid_variants(result, sequence, cds_regions)
#
#         write_to_csv(out)
#         print(out)
#         self.assertIsInstance(result, dict)
#
# if __name__ == '__main__':
#     unittest.main()
