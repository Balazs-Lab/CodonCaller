import unittest
from CodonCaller.haplotype import *


class TestYourCodeWT(unittest.TestCase):

    def setUp(self):
        bed_file = "../test_data/bed_file/JRCSF_haplotypes.bed"
        bam_file = "../test_data/sorted_reads/CD27-m787-20-JRCSF.bam"
        reference_fasta = "../test_data/reference_sequence/JRCSF-reference.fa"
        self.test_instance = Haplotype(bed_file, bam_file, reference_fasta)

    def test_load_haplotype_regions(self):
        print(self.test_instance.haplotype_regions)
        self.assertEqual(self.test_instance.haplotype_regions['V1V2'][0], 442)
        self.assertEqual(self.test_instance.haplotype_regions['D'][0], 865)

    def test_fetch_reads(self):
        aa_count = self.test_instance.fetch_sequencing_data(contig='JRCSF',region_start=442, region_end= 631)
        df.to_csv("../test_data/CD27-m787-20-JRCSF-haplotype.csv")

        # print(reads[1850])
        # read = reads[1850]
        # print(read.reference_start)
        # print(read.reference_end)


    def test_all_haplotyoes(self):

        output = self.test_instance.all_haplotyoes()

        output.to_csv("../test_data/annotation.csv")

        print(output)