import unittest
from codon_call import *
import pysam
from Bio import SeqIO
import csv


class TestYourCode(unittest.TestCase):

    def test_load_cds_regions(self):
        bed_file = "../test_data/bed_file/REJOc.bed"
        cds_regions = load_cds_regions(bed_file)
        print(cds_regions)
        self.assertEqual(len(cds_regions), 1)  # Adjust based on the content of your test BED file

    def test_translate_csf(self):
        reference_sequence = "../test_data/reference_sequence/REJOc-reference.fa"
        bed_file = "../test_data/bed_file/REJOc.bed"
        cds_regions = load_cds_regions(bed_file)

        sequence = read_fasta(reference_sequence)
        cds_regions = load_cds_regions(bed_file)
        cds_start = cds_regions[0][1]
        cds_end = cds_regions[0][2]

        out = translate_cds(sequence, cds_start, cds_end)
        print(str(out))
        # self.assertEqual()

    def test_is_within_cds(self):
        position_within_cds = 150
        position_outside_cds = 30
        bed_file = "../test_data/bed_file/REJOc.bed"
        cds_regions = load_cds_regions(bed_file)
        self.assertTrue(is_within_cds(position_within_cds, cds_regions))
        self.assertFalse(is_within_cds(position_outside_cds, cds_regions))

    def test_read_frame_0(self):
        bed_file = "../test_data/bed_file/REJOc.bed"
        cds_regions = load_cds_regions(bed_file)
        start_position = cds_regions[0][1]
        print(start_position)

        # Simulate Aligned Segments
        mapping_quality = 25

        ## in frame
        sequence = "ATGAAAGTGAAGGGGATCAGGAGGAATTAT"
        positions = list(range(53, 53 + len(sequence)))
        qualty = len(sequence) * [25]

        print(sequence)
        print(positions)

        codon_positions, codons, qualty_sum = read_to_codon(start_position,positions,sequence,qualty)
        print(codon_positions)
        print(codons)
        print(qualty_sum)

        self.assertTrue(codon_positions, list(range(1, 10)))

    def test_read_frame_1(self):
        bed_file = "../test_data/bed_file/REJOc.bed"
        cds_regions = load_cds_regions(bed_file)
        start_position = cds_regions[0][1]
        print(start_position)

        # Simulate Aligned Segments
        mapping_quality = 25

        ## in frame
        sequence = "TGAAAGTGAAGGGGATCAGGAGGAATTATC"
        positions = list(range(54, 54 + len(sequence)))
        qualty = len(sequence) * [25]

        print(sequence)
        print(positions)

        codon_positions, codons, qualty_sum = read_to_codon(start_position,positions,sequence,qualty)
        print(codon_positions)
        print(codons)
        print(qualty_sum)

        self.assertTrue(codon_positions, list(range(2, 10)))

    def test_read_frame_2(self):
        bed_file = "../test_data/bed_file/REJOc.bed"
        cds_regions = load_cds_regions(bed_file)
        start_position = cds_regions[0][1]
        print(start_position)

        # Simulate Aligned Segments
        mapping_quality = 25

        ## in frame 2
        sequence = "GAAAGTGAAGGGGATCAGGAGGAATTATCA"
        positions = list(range(55, 55 + len(sequence)))
        qualty = len(sequence) * [25]

        print(sequence)
        print(positions)

        codon_positions, codons, qualty_sum = read_to_codon(start_position,positions,sequence,qualty)
        print(codon_positions)
        print(codons)
        print(qualty_sum)

        self.assertTrue(codon_positions, list(range(2, 10)))

    def test_read_non_coding(self):
        bed_file = "../test_data/bed_file/REJOc.bed"
        cds_regions = load_cds_regions(bed_file)
        start_position = cds_regions[0][1]
        print(start_position)

        # Simulate Aligned Segments
        mapping_quality = 25

        ## in frame
        sequence = "AAAATAGACAGGTTAATTGATAGAATAAGAGATAGAGCAGAAGACAGTGGCAATGAAAGTGAAGGGGATC"
        positions = list(range(1, 1 + len(sequence)))
        qualty = len(sequence) * [35]

        print(sequence)
        print(positions)

        codon_positions, codons, quality_sum = read_to_codon(start_position,positions,sequence,qualty)
        d = dict(zip(codon_positions, zip(codons, quality_sum)))
        print(d)
        # print(codons)
        # print(qualty_sum)

        self.assertTrue(codon_positions, list(range(1, 10)))

    def test_process_read(self):
        # Create test data or mock objects as needed
        bam_file = "../test_data/mapped_reads/CD5-m246-20-REJOc.bam"
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            reads = list(bam.fetch())
        read_object = reads[123]

        # Load reference sequence
        reference_seq = SeqIO.read("../test_data/reference_index/REJOc-reference.fa", "fasta").seq

        quality_threshold = 20

        # Load bed regions
        bed_file = "../test_data/bed_file/REJOc.bed"
        cds_regions = load_cds_regions(bed_file)

        result = process_read(read_object, quality_threshold, cds_regions)
        self.assertIsInstance(result, dict)

    def test_process_sample(self):
        # Create test data or mock objects as needed
        bam_file = "../test_data/mapped_reads/CD5-m246-20-REJOc.bam"

        quality_threshold = 40

        # Load bed regions
        bed_file = "../test_data/bed_file/REJOc.bed"
        cds_regions = load_cds_regions(bed_file)

        result = process_sample(bam_file, bed_file, quality_threshold)
        self.assertIsInstance(result, dict)

    def test_process_sample_JRCSF(self):
        # Create test data or mock objects as needed
        bam_file = "../test_data/mapped_reads/CD21-m009-12-JRCSF.bam"

        quality_threshold = 40

        # Load bed regions
        bed_file = "../test_data/bed_file/JRCSF.bed"
        cds_regions = load_cds_regions(bed_file)

        result = process_sample(bam_file, bed_file, quality_threshold)
        self.assertIsInstance(result, dict)


    def test_process_sample_call(self):
        # Create test data or mock objects as needed
        bam_file = "../test_data/mapped_reads/CD5-m246-20-REJOc.bam"
        reference_sequence = "../test_data/reference_sequence/REJOc-reference.fa"
        sequence = read_fasta(reference_sequence)
        bed_file = "../test_data/bed_file/REJOc.bed"
        quality_threshold = 40

        cds_regions = load_cds_regions(bed_file)

        result = process_sample(bam_file, bed_file, quality_threshold)

        out = call_amino_acid_variants(result, sequence, cds_regions)
        print(out)
        self.assertIsInstance(result, dict)


if __name__ == '__main__':
    unittest.main()
