#!/usr/bin/env python

import argparse
from argparse import Namespace
from CodonCaller import VaricantCaller

"""
Runs Codon Caller Variant Call Functions
"""

def main():
    parser = argparse.ArgumentParser(description='Description of your tool')
    parser.add_argument('bed_file', help='bed_file')
    parser.add_argument('bam_file', help='bam_file')
    parser.add_argument('cds', help='You MUST specify which CDS you want to call from')
    parser.add_argument('reference_fasta', help='reference_fasta')
    parser.add_argument('read_quality_threshold', help='read_quality_threshold')
    parser.add_argument('base_quality_threshold', help='base_quality_threshold')
    arser.add_argument('output', help='base_quality_threshold')

    parser.add_argument('--optional_arg', help='Description of optional_arg')

    args: Namespace = parser.parse_args()

    File = VariantCaller( args.bed_file, args.bam_file, args.reference_fasta, args.read_quality_threshold, args.base_quality_threshold=35)

    File.process_sample(args.cds)
    File.aa_count_to_freq(args.cds)
    File.write_to_csv(args.output)


if __name__ == '__main__':
    main() # run the script