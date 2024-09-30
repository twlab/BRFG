#!/usr/bin/env python3

import pysam
import argparse


def fix_chrom_names(input_bam, output_bam, prefix):
    # Open input BAM for reading
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        # Modify the header to remove the prefix from reference names
        new_header = infile.header.to_dict()
        reference_name_map = {}

        for seq in new_header['SQ']:
            if seq['SN'].startswith(prefix):
                new_name = seq['SN'].replace(prefix, "", 1)
                reference_name_map[seq['SN']] = new_name
                seq['SN'] = new_name

        # Open output BAM for writing, using the modified header
        with pysam.AlignmentFile(output_bam, "wb", header=new_header) as outfile:
            for read in infile:
                # Update the reference_id based on the modified header
                if read.reference_name in reference_name_map:
                    read.reference_id = outfile.get_tid(reference_name_map[read.reference_name])

                if read.next_reference_name in reference_name_map:
                    read.next_reference_id = outfile.get_tid(reference_name_map[read.next_reference_name])

                # Write the read with the updated reference IDs
                outfile.write(read)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fix chromosome names in BAM files by removing a specified prefix.")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bam", help="Output BAM file with fixed chromosome names")
    parser.add_argument("--prefix", type=str, default="GRCh38.",
                        help="Prefix string to remove from chromosome names (default: 'GRCh38.')")

    args = parser.parse_args()

    fix_chrom_names(args.input_bam, args.output_bam, args.prefix)
