import pysam
import pandas as pd
from Bio.Seq import Seq
from collections import Counter


class Haplotype:

    def __init__(self, bed_file, bam_file, reference_fasta):
        self.bed_file = bed_file
        self.bam_file = bam_file
        self.fasta_file = reference_fasta
        # self.read_quality_threshold = read_quality_threshold
        # self.base_quality_threshold = base_quality_threshold * 3

        # load in regions for analysis
        self.haplotype_regions = self.load_haplotype_regions()

        # read in mapped data

        # translate mapped data

        # sort mapped data into haplotypes

        # amino acid alignment

    def load_haplotype_regions(self):
        """
        :return: dict of haplotype regions - name, start, stop
        Analysis will be performed for each region
        Haplotype regions should not be longer than 300 bp to ensure sufficient coverage from short reads
        """
        haplotype_regions = dict()
        with open(self.bed_file, 'r') as bed:
            for line in bed:
                chrom, start, end = line.strip().split('\t')
                haplotype_regions[chrom] = (int(start) - 1, int(end) - 1)
        return haplotype_regions

    def fetch_sequencing_data(self, contig, region_start, region_end):
        """
        Fetches reads that overlap with the specified region
        :return:
        """
        nt_list = list()
        aa_list = list()
        # print(region_start, region_end)
        bamfile = pysam.AlignmentFile(self.bam_file, "rb")
        reads = list()
        iter = bamfile.fetch(contig=contig, start=region_start, stop=region_end)
        for x in iter:
            if x.reference_start <= region_start:
                if x.reference_end >= region_end:
                    # print(x.query_sequence())
                    reads.append(x)

        # trim reads
        for read in reads:
            read_sequence = read.query_sequence
            cigar_dict = dict(read.cigartuples)
            deletions = 0
            if 2 in cigar_dict.keys():
                deletions = cigar_dict[2]
            trim_seq = read_sequence[region_start - read.reference_start: region_end - read.reference_start - deletions]
            nt_list.append(trim_seq)
            # print(trim_seq)

            # translate reads
            # print(Seq(trim_seq).translate())
            aa = Seq(trim_seq).translate()
            aa_list.append(str(aa))

        aa_count = Counter(aa_list)
        total = sum(aa_count.values(),0.0)
        for key in aa_count:
            aa_count[key] /= total
        # # print(aa_count)
        # df = pd.DataFrame.from_dict(aa_count,orient="index").reset_index()
        # # print(df)
        return aa_count

    def all_haplotyoes(self):
        haplotypes_out = dict()
        haplotype_df = pd.DataFrame()
        # haplotype_df.columns = ['haplotype', 'count']
        for key, value in self.haplotype_regions.items():
            # print(key, value[0], value [1])
            aa_count = self.fetch_sequencing_data(contig='JRCSF', region_start=value[0], region_end=value[1])
            df = pd.DataFrame.from_dict(aa_count, orient="index").reset_index()
            df['site'] = key
            # print(df.columns)
            # df.columns = ['haplotype', 'count']
            # print(df)
            haplotype_df = haplotype_df.append(df, ignore_index=True)

        haplotype_df.columns = ['haplotype', 'count', 'site']
        return haplotype_df
