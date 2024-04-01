import pysam
import pandas as pd
from Bio.Seq import Seq
from collections import Counter

class Haplotype:

    def __init__(self, bed_file, bam_file, reference_fasta):
        self.bed_file = bed_file
        self.bam_file = bam_file
        self.fasta_file = reference_fasta
        self.contig = self.fasta_file.split("/")[-1].split("-")[0]

        # self.read_quality_threshold = read_quality_threshold
        # self.base_quality_threshold = base_quality_threshold * 3

        # load in regions for analysis
        self.haplotype_regions = self.load_haplotype_regions()

        self.haplotype_df = None

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
            if isinstance(x.reference_start, int) and isinstance(x.reference_end, int):
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
        total = sum(aa_count.values(), 0.0)
        for key in aa_count:
            aa_count[key] /= total

        return aa_count

    @property
    def all_haplotypes(self):
        haplotypes_out = list()
        counts_out = list()
        sites_out = list()

        for key, value in self.haplotype_regions.items():
            aa_count = self.fetch_sequencing_data(contig=self.contig, region_start=value[0], region_end=value[1])
            haplotypes_out.extend(list(aa_count.keys()))
            counts_out.extend(list(aa_count.values()))
            sites_out.extend([key] * len(list(aa_count.keys())))

        haplotype_df = pd.DataFrame(list(zip(haplotypes_out, counts_out, sites_out)),
                                    columns=['haplotype', 'count', 'site'])
        self.haplotype_df = haplotype_df
        return haplotype_df

    def write_to_csv(self,out_file):

        self.haplotype_df.to_csv(out_file)

