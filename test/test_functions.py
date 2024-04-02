import unittest
from CodonCaller.variant_caller import *
import pysam
from collections import Counter


### MisMap Testing

class TestMisMap(unittest.TestCase):


    def setUp(self):
        bed_file = "../test_data/bed_file/JRCSF.bed"
        bam_file = "../test_data/sorted_reads/CD27-m787-20-JRCSF.bam"
        reference_fasta = "../test_data/reference_sequence/JRCSF-reference.fa"
        read_quality_threshold = 40
        base_quality_threshold = 30
        self.test_instance = VariantCaller(bed_file, bam_file, reference_fasta, read_quality_threshold,
                                           base_quality_threshold)

    def test_codon_mistmap(self):
        read = self.test_instance.reads[2500]
        cds = "JRCSF"

        codon_starts = [i for i in range(52, 2598, 3)]
        read_sequence, read_positions, quality, indels = parse_read(read)
        print(read_sequence)
        print(read_positions)
        print(indels)

        read_sequence, read_positions, quality, indels = self.test_instance.map_correct(read)
        print(read_sequence)
        print(read_positions)
        print(indels)

        codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions, quality, indels)
        print(codon_positions)
        print(codons)

    def test_pipeline(self):
        cds = "JRCSF"
        self.test_instance.process_sample(cds)
        self.test_instance.aa_count_to_freq(cds)
        print(self.test_instance.aa_counts)
        self.test_instance.write_to_csv("../test_data/calls/map_test.csv")

### WT Testing ###
class TestYourCodeWT(unittest.TestCase):

    def setUp(self):
        bed_file = "../test_data/bed_file/REJOc_mini.bed"
        bam_file = "../test_data/mapped_reads/CD00-m000-WT-REJOc.bam"
        reference_fasta = "../test_data/reference_sequence/REJOc-mini.fa"
        read_quality_threshold = 40
        base_quality_threshold = 30
        self.test_instance = VariantCaller(bed_file, bam_file, reference_fasta, read_quality_threshold,
                                           base_quality_threshold)

    def test_load_cds_regions(self):
        self.test_instance.load_cds_regions()
        self.assertEqual(self.test_instance.cds_regions['Env'][0], 0)
        self.assertEqual(self.test_instance.cds_regions['Env'][1], 585)

    def test_read_fasta(self):
        self.test_instance.read_fasta()
        self.assertEqual(len(self.test_instance.reference_sequence), 597)

    def test_translate_cds(self):
        self.test_instance.translate_cds()
        print(self.test_instance.cds_aminoacids['Env'])
        print(len(self.test_instance.cds_aminoacids['Env']))
        self.assertEqual(len(self.test_instance.cds_aminoacids['Env']), 195)

    def test_parse_read(self):
        read = self.test_instance.reads[1]
        cds = "Env"
        print(read.seq)
        read_sequence, read_positions, quality, indels = parse_read(read)
        print(indels)
        print(read_positions)
        print(quality)
        print(read_sequence)

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

    def test_compiled_reads(self):
        cds = "Env"
        out_dict = self.test_instance.compile_reads(cds)
        print(out_dict)
        # print(len(out_dict['77i']) + len(out_dict['77']))
        clean_dict = fix_insertions(out_dict)
        print(clean_dict)
        # print(len(clean_dict['77']))
        print(clean_dict.values())

    def test_process_samples(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        print(self.test_instance.aa_counts)

    def test_process_aa_count_to_freq_WT(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        self.test_instance.aa_count_to_freq(cds)
        print(self.test_instance.mut_call)

    def test_pipeline(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        self.test_instance.aa_count_to_freq(cds)
        self.test_instance.write_to_csv("../test_data/calls/WT_test.csv")


### SNP Testing ###
class TestYourCodeSNP(unittest.TestCase):

    def setUp(self):
        bed_file = "../test_data/bed_file/REJOc_mini.bed"
        bam_file = "../test_data/mapped_reads/CD00-m000-SNP-REJOc.bam"
        reference_fasta = "../test_data/reference_sequence/REJOc-mini.fa"
        read_quality_threshold = 40
        base_quality_threshold = 30
        self.test_instance = VariantCaller(bed_file, bam_file, reference_fasta, read_quality_threshold,
                                           base_quality_threshold)

    def test_parse_read(self):
        read = self.test_instance.reads[1]
        cds = "Env"
        print(read.seq)
        read_sequence, read_positions, quality, indels = parse_read(read)
        print(indels)
        print(read_positions)
        print(quality)
        print(read_sequence)

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

    def test_compiled_reads(self):
        cds = "Env"
        out_dict = self.test_instance.compile_reads(cds)
        print(out_dict)
        # print(len(out_dict['77i']) + len(out_dict['77']))
        clean_dict = fix_insertions(out_dict)
        print(clean_dict)
        # print(len(clean_dict['77']))
        print(clean_dict.values())

    def test_process_samples(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        print(self.test_instance.aa_counts)

    def test_process_aa_count_to_freq(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        self.test_instance.aa_count_to_freq(cds)
        print(self.test_instance.mut_call)

    def test_pipeline(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        self.test_instance.aa_count_to_freq(cds)
        self.test_instance.write_to_csv("../test_data/calls/SNP_test.csv")


### Insertion Testing ###

class TestYourCodeIN(unittest.TestCase):

    def setUp(self):
        bed_file = "../test_data/bed_file/REJOc_mini.bed"
        bam_file = "../test_data/mapped_reads/CD00-m000-IN-REJOc.bam"
        reference_fasta = "../test_data/reference_sequence/REJOc-mini.fa"
        read_quality_threshold = 30
        base_quality_threshold = 30
        self.test_instance = VariantCaller(bed_file, bam_file, reference_fasta, read_quality_threshold,
                                           base_quality_threshold)

    def test_parse_read(self):
        read = self.test_instance.reads[1]
        cds = "Env"
        print(read.seq)
        read_sequence, read_positions, quality, indels = parse_read(read)
        print(indels)
        print(read_positions)
        print(quality)
        print(read_sequence)

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

    def test_compiled_reads(self):
        cds = "Env"
        out_dict = self.test_instance.compile_reads(cds)
        print(out_dict)
        # print(len(out_dict['77i']) + len(out_dict['77']))
        clean_dict = fix_insertions(out_dict)
        print(clean_dict)
        # print(len(clean_dict['77']))
        print(clean_dict.values())

    def test_process_samples(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        print(self.test_instance.aa_counts)

    def test_process_aa_count_to_freq(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        self.test_instance.aa_count_to_freq(cds)
        print(self.test_instance.mut_call)

    def test_pipeline(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        self.test_instance.aa_count_to_freq(cds)
        self.test_instance.write_to_csv("../test_data/calls/IN_test.csv")

    #
    # def test_read_to_codon(self):
    #     cds = "Env"
    #     read_sequence = "ATGAAAGTGAAGGGGATCAGGAGGAATTAT"
    #     read_positions = list(range(53, 53 + len(read_sequence)))
    #     quality = len(read_sequence) * [40]
    #     indels = {"insertion": [(3, 3)], "deletion": [(6, 3)]}
    #     codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
    #                                                                             quality, indels)
    #     print(codon_positions)
    #     print(codons)
    #     print(quality_sum)
    #
    # def test_read_to_codon_frame_1(self):
    #     cds = "Env"
    #     read_sequence = "TGAAAGTGAAGGGGATCAGGAGGAATTATC"
    #     read_positions = list(range(54, 54 + len(read_sequence)))
    #     quality = len(read_sequence) * [40]
    #     codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
    #                                                                             quality)
    #     print(codon_positions)
    #     print(codons)
    #     print(quality_sum)
    #
    # def test_read_to_codon_frame_2(self):
    #     cds = "Env"
    #     read_sequence = "GAAAGTGAAGGGGATCAGGAGGAATTATCA"
    #     read_positions = list(range(55, 55 + len(read_sequence)))
    #     quality = len(read_sequence) * [40]
    #     codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
    #                                                                             quality)
    #     print(codon_positions)
    #     print(codons)
    #     print(quality_sum)
    #
    # def test_read_to_codon_real_read(self):
    #     cds = "Env"
    #     read = self.test_instance.reads[5]
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #     codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
    #                                                                             quality, indels)
    #     print(codon_positions)
    #     print(codons)
    #     # print(quality_sum)
    #
    # def test_fix_insert(self):
    #     read = self.test_instance.reads[6]
    #     cds = "Env"
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #     insertion = indels["insertion"][0]
    #     out = fix_insert(read_positions, insertion)
    #     print(out)
    #
    # def test_fix_deletion(self):
    #     read = self.test_instance.reads[6]
    #     cds = "Env"
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #     deletion = indels["deletion"][0]
    #     seq, pos = fix_deletion(read_sequence, read_positions, deletion)
    #     print(deletion)
    #     print(read_positions)
    #     print(read_sequence)
    #     print(seq)
    #     print(pos)
    #
    # def test_fix_indels(self):
    #     read = self.test_instance.reads[6]
    #     cds = "Env"
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #     print(indels)
    #     print(read_positions)
    #     print(quality)
    #     print(read_sequence)
    #
    #     read_sequence, read_positions, quality = fix_indels(read_sequence, read_positions, quality, indels)
    #
    #     print(read_positions)
    #     print(read_sequence)
    #     print(quality)
    #
    # def test_process_read(self):
    #     read = self.test_instance.reads[6]
    #     cds = "Env"
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #     print(indels)
    #     print(read_positions)
    #     print(read_sequence)
    #     codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
    #                                                                             quality, indels)
    #     print(codon_positions)
    #     print(codons)
    #     print(quality_sum)
    #
    # def test_process_read2(self):
    #     read = self.test_instance.reads[6]
    #     cds = "Env"
    #
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #
    #     codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
    #                                                                             quality, indels)
    #     # print(read_sequence)
    #     # print(read_positions)
    #     # print(codon_positions)
    #     # print(codons)
    #     # print(quality_sum)
    #
    #     dict_out = self.test_instance.process_read(cds, read)
    #     print(dict_out)
    #
    # def test_process_sample(self):
    #     cds = "Env"
    #     self.test_instance.process_sample(cds)
    #     print(self.test_instance.mut_call)
    #
    # def test_process_sample(self):
    #     cds = "Env"
    #     self.test_instance.process_sample(cds)
    #     print(self.test_instance.aa_counts)

    # def test_compiled_reads(self):
    #     cds = "Env"
    #     out_dict = self.test_instance.compile_reads(cds)
    #     print(len(out_dict['77i']) + len(out_dict['77']))
    #     clean_dict = fix_insertions(out_dict)
    #     print(len(clean_dict['77']))
    #     # print(out_dict.values())

    #
    # def test_write_data(self):
    #     cds = "Env"
    #     self.test_instance.process_sample(cds)
    #     self.test_instance.aa_count_to_freq(cds)
    #     self.test_instance.write_to_csv("../test_data/calls/test.csv")
    #
    # def test_coverage(self):
    #     cds = "Env"
    #     self.test_instance.process_sample(cds)
    #     self.test_instance.aa_count_to_freq(cds)
    #     print(self.test_instance.aa_coverage)
    #     self.test_instance.write_coverage_to_csv("../test_data/calls/CD20-m999-8-REJOc-coverage.csv")
    #
    # def test_INDEL(self):
    #     cds = "Env"
    #     self.test_instance.process_sample(cds)
    #     self.test_instance.aa_count_to_freq(cds)
    #     print(self.test_instance.aa_coverage)
    #     # self.test_instance.aa_count_to_freq(cds)
    #
    # def test_fix_insert(self):
    #     read = self.test_instance.reads[5]
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #     sorted_indel = sort_indel_dict(indels)
    #     print(sorted_indel)
    #     corrected_positions = read_positions
    #     for position in sorted_indel.keys():
    #         if sorted_indel[position][0] == "insertion":
    #             print(position)
    #             start = position
    #             length = sorted_indel[position][1]
    #             corrected_positions = fix_insert(corrected_positions, start, length, del_shift=0)
    #     print(read_positions)
    #     count = Counter(corrected_positions)
    #     print(count)
    #     assert (count[58] == 2) & (count[142] == 2)
    #
    # def test_fix_indels(self):
    #     read = self.test_instance.reads[6]
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #     print(indels)
    #     print(read_sequence)
    #     print(read_positions)
    #     # print(quality)
    #     read_sequence, read_positions, quality = fix_indels(read_sequence, read_positions, quality, indels)
    #     print(read_sequence)
    #     print(read_positions)
    #     # print(quality)
    #
    #
    # def test_fix_deletion(self):
    #
    #     read = self.test_instance.reads[5]
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #     sorted_indel = sort_indel_dict(indels)
    #     print(sorted_indel)
    #     corrected_positions = read_positions
    #     total_insert_shift = 3
    #     for position in sorted_indel.keys():
    #         if sorted_indel[position][0] == "deletion":
    #             # print('del:',position)
    #             start = position
    #             length = sorted_indel[position][1]
    #             shift = total_insert_shift
    #             read_sequence, read_positions, quality = fix_deletion(read_sequence, read_positions, quality, start, length, shift)
    #     print(read_positions)


### Deletion Testing ###

class TestYourCodeDEL(unittest.TestCase):

    def setUp(self):
        bed_file = "../test_data/bed_file/REJOc_mini.bed"
        bam_file = "../test_data/mapped_reads/CD00-m000-DEL-REJOc.bam"
        reference_fasta = "../test_data/reference_sequence/REJOc-mini.fa"
        read_quality_threshold = 30
        base_quality_threshold = 30
        self.test_instance = VariantCaller(bed_file, bam_file, reference_fasta, read_quality_threshold,
                                           base_quality_threshold)

    def test_parse_read(self):
        read = self.test_instance.reads[1]
        cds = "Env"
        print(read.seq)
        read_sequence, read_positions, quality, indels = parse_read(read)
        print(indels)
        print(read_positions)
        print(quality)
        print(read_sequence)

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

    def test_compiled_reads(self):
        cds = "Env"
        out_dict = self.test_instance.compile_reads(cds)
        print(out_dict)
        # print(len(out_dict['77i']) + len(out_dict['77']))
        clean_dict = fix_insertions(out_dict)
        print(clean_dict)
        # print(len(clean_dict['77']))
        print(clean_dict.values())

    def test_process_samples(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        print(self.test_instance.aa_counts)

    def test_process_aa_count_to_freq(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        self.test_instance.aa_count_to_freq(cds)
        print(self.test_instance.mut_call)

    def test_pipeline(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        self.test_instance.aa_count_to_freq(cds)
        self.test_instance.write_to_csv("../test_data/calls/DEL_test.csv")

    #
    # def test_read_to_codon(self):
    #     cds = "Env"
    #     read_sequence = "ATGAAAGTGAAGGGGATCAGGAGGAATTAT"
    #     read_positions = list(range(53, 53 + len(read_sequence)))
    #     quality = len(read_sequence) * [40]
    #     indels = {"insertion": [(3, 3)], "deletion": [(6, 3)]}
    #     codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
    #                                                                             quality, indels)
    #     print(codon_positions)
    #     print(codons)
    #     print(quality_sum)
    #
    # def test_read_to_codon_frame_1(self):
    #     cds = "Env"
    #     read_sequence = "TGAAAGTGAAGGGGATCAGGAGGAATTATC"
    #     read_positions = list(range(54, 54 + len(read_sequence)))
    #     quality = len(read_sequence) * [40]
    #     codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
    #                                                                             quality)
    #     print(codon_positions)
    #     print(codons)
    #     print(quality_sum)
    #
    # def test_read_to_codon_frame_2(self):
    #     cds = "Env"
    #     read_sequence = "GAAAGTGAAGGGGATCAGGAGGAATTATCA"
    #     read_positions = list(range(55, 55 + len(read_sequence)))
    #     quality = len(read_sequence) * [40]
    #     codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
    #                                                                             quality)
    #     print(codon_positions)
    #     print(codons)
    #     print(quality_sum)
    #
    # def test_read_to_codon_real_read(self):
    #     cds = "Env"
    #     read = self.test_instance.reads[5]
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #     codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
    #                                                                             quality, indels)
    #     print(codon_positions)
    #     print(codons)
    #     # print(quality_sum)
    #
    # def test_fix_insert(self):
    #     read = self.test_instance.reads[6]
    #     cds = "Env"
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #     insertion = indels["insertion"][0]
    #     out = fix_insert(read_positions, insertion)
    #     print(out)
    #
    # def test_fix_deletion(self):
    #     read = self.test_instance.reads[6]
    #     cds = "Env"
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #     deletion = indels["deletion"][0]
    #     seq, pos = fix_deletion(read_sequence, read_positions, deletion)
    #     print(deletion)
    #     print(read_positions)
    #     print(read_sequence)
    #     print(seq)
    #     print(pos)
    #
    # def test_fix_indels(self):
    #     read = self.test_instance.reads[6]
    #     cds = "Env"
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #     print(indels)
    #     print(read_positions)
    #     print(quality)
    #     print(read_sequence)
    #
    #     read_sequence, read_positions, quality = fix_indels(read_sequence, read_positions, quality, indels)
    #
    #     print(read_positions)
    #     print(read_sequence)
    #     print(quality)
    #
    # def test_process_read(self):
    #     read = self.test_instance.reads[6]
    #     cds = "Env"
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #     print(indels)
    #     print(read_positions)
    #     print(read_sequence)
    #     codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
    #                                                                             quality, indels)
    #     print(codon_positions)
    #     print(codons)
    #     print(quality_sum)
    #
    # def test_process_read2(self):
    #     read = self.test_instance.reads[6]
    #     cds = "Env"
    #
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #
    #     codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
    #                                                                             quality, indels)
    #     # print(read_sequence)
    #     # print(read_positions)
    #     # print(codon_positions)
    #     # print(codons)
    #     # print(quality_sum)
    #
    #     dict_out = self.test_instance.process_read(cds, read)
    #     print(dict_out)
    #
    # def test_process_sample(self):
    #     cds = "Env"
    #     self.test_instance.process_sample(cds)
    #     print(self.test_instance.mut_call)
    #
    # def test_process_sample(self):
    #     cds = "Env"
    #     self.test_instance.process_sample(cds)
    #     print(self.test_instance.aa_counts)

    # def test_compiled_reads(self):
    #     cds = "Env"
    #     out_dict = self.test_instance.compile_reads(cds)
    #     print(len(out_dict['77i']) + len(out_dict['77']))
    #     clean_dict = fix_insertions(out_dict)
    #     print(len(clean_dict['77']))
    #     # print(out_dict.values())

    #
    # def test_write_data(self):
    #     cds = "Env"
    #     self.test_instance.process_sample(cds)
    #     self.test_instance.aa_count_to_freq(cds)
    #     self.test_instance.write_to_csv("../test_data/calls/test.csv")
    #
    # def test_coverage(self):
    #     cds = "Env"
    #     self.test_instance.process_sample(cds)
    #     self.test_instance.aa_count_to_freq(cds)
    #     print(self.test_instance.aa_coverage)
    #     self.test_instance.write_coverage_to_csv("../test_data/calls/CD20-m999-8-REJOc-coverage.csv")
    #
    # def test_INDEL(self):
    #     cds = "Env"
    #     self.test_instance.process_sample(cds)
    #     self.test_instance.aa_count_to_freq(cds)
    #     print(self.test_instance.aa_coverage)
    #     # self.test_instance.aa_count_to_freq(cds)
    #
    # def test_fix_insert(self):
    #     read = self.test_instance.reads[5]
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #     sorted_indel = sort_indel_dict(indels)
    #     print(sorted_indel)
    #     corrected_positions = read_positions
    #     for position in sorted_indel.keys():
    #         if sorted_indel[position][0] == "insertion":
    #             print(position)
    #             start = position
    #             length = sorted_indel[position][1]
    #             corrected_positions = fix_insert(corrected_positions, start, length, del_shift=0)
    #     print(read_positions)
    #     count = Counter(corrected_positions)
    #     print(count)
    #     assert (count[58] == 2) & (count[142] == 2)
    #
    # def test_fix_indels(self):
    #     read = self.test_instance.reads[6]
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #     print(indels)
    #     print(read_sequence)
    #     print(read_positions)
    #     # print(quality)
    #     read_sequence, read_positions, quality = fix_indels(read_sequence, read_positions, quality, indels)
    #     print(read_sequence)
    #     print(read_positions)
    #     # print(quality)
    #
    #
    # def test_fix_deletion(self):
    #
    #     read = self.test_instance.reads[5]
    #     read_sequence, read_positions, quality, indels = parse_read(read)
    #     sorted_indel = sort_indel_dict(indels)
    #     print(sorted_indel)
    #     corrected_positions = read_positions
    #     total_insert_shift = 3
    #     for position in sorted_indel.keys():
    #         if sorted_indel[position][0] == "deletion":
    #             # print('del:',position)
    #             start = position
    #             length = sorted_indel[position][1]
    #             shift = total_insert_shift
    #             read_sequence, read_positions, quality = fix_deletion(read_sequence, read_positions, quality, start, length, shift)
    #     print(read_positions)


### INDEL Testing ###
class TestYourCodeINDEL(unittest.TestCase):

    def setUp(self):
        bed_file = "../test_data/bed_file/REJOc_mini.bed"
        bam_file = "../test_data/mapped_reads/CD00-m000-INDEL-REJOc.bam"
        reference_fasta = "../test_data/reference_sequence/REJOc-mini.fa"
        read_quality_threshold = 30
        base_quality_threshold = 30
        self.test_instance = VariantCaller(bed_file, bam_file, reference_fasta, read_quality_threshold,
                                           base_quality_threshold)

    def test_parse_read(self):
        read = self.test_instance.reads[1]
        cds = "Env"
        print(read.seq)
        read_sequence, read_positions, quality, indels = parse_read(read)
        print(indels)
        print(read_positions)
        print(quality)
        print(read_sequence)

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

    def test_read_to_codon(self):
        read = self.test_instance.reads[6]
        cds = "Env"
        read_sequence, read_positions, quality, indels = parse_read(read)

        print(indels)
        # fix indels
        codon_positions, codons, quality_sum = self.test_instance.read_to_codon(cds, read_sequence, read_positions,
                                                                                quality, indels)
        print(codon_positions)
        print(codons)

    def test_compiled_reads(self):
        cds = "Env"
        out_dict = self.test_instance.compile_reads(cds)
        print(out_dict)
        # print(len(out_dict['77i']) + len(out_dict['77']))
        clean_dict = fix_insertions(out_dict)
        print(clean_dict)
        # print(len(clean_dict['77']))
        print(clean_dict.values())

    def test_process_samples(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        print(self.test_instance.aa_counts)

    def test_process_aa_count_to_freq(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        self.test_instance.aa_count_to_freq(cds)
        print(self.test_instance.mut_call)

    def test_pipeline(self):
        cds = "Env"
        self.test_instance.process_sample(cds)
        self.test_instance.aa_count_to_freq(cds)
        self.test_instance.write_to_csv("../test_data/calls/INDEL_test.csv")
