#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import gzip
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Extract 5pUMIs from SmartSeq-FLA-UMI reads (Cogent-style)")
    parser.add_argument("-r1", "--read1", required=True, help="Path to input Read 1 FASTQ (can be .gz)")
    parser.add_argument("-r2", "--read2", required=True, help="Path to input Read 2 FASTQ (can be .gz)")
    parser.add_argument("-t", "--tag_map", required=True, help="Path to recognition_tag_map_5pUMI_with_N.csv")
    parser.add_argument("-o", "--out_prefix", required=True, help="Output directory for processed FASTQs")
    return parser.parse_args()

def open_fastq(filename):
    return gzip.open(filename, "rt") if filename.endswith(".gz") else open(filename, "r")

def main():
    args = parse_args()

    UMI_TAG_LEN = 11
    UMI_LEN = 8
    EXTRA_TRIM = 4

    # Load valid UMI tags
    with open(args.tag_map) as f:
        tag_set = set(line.strip() for line in f if line.strip())

    # Get output filenames
    input_prefix = os.path.basename(args.read1).split("_R1")[0]
    out_r1 = os.path.join(args.out_prefix, f"{input_prefix}_R1.fastq")
    out_r2 = os.path.join(args.out_prefix, f"{input_prefix}_R2.fastq")

    with open_fastq(args.read1) as f1, open_fastq(args.read2) as f2, \
         open(out_r1, "w") as out1, open(out_r2, "w") as out2:

        r1_iter = SeqIO.parse(f1, "fastq")
        r2_iter = SeqIO.parse(f2, "fastq")

        for i, (r1, r2) in enumerate(zip(r1_iter, r2_iter), 1):
            seq = str(r1.seq)
            tag = seq[:UMI_TAG_LEN]

            if tag in tag_set:
                umi = seq[UMI_TAG_LEN:UMI_TAG_LEN + UMI_LEN]
                if "N" not in umi:
                    # Trim R1
                    trim_len = UMI_TAG_LEN + UMI_LEN + EXTRA_TRIM
                    trimmed_seq = seq[trim_len:]
                    trimmed_qual = r1.letter_annotations["phred_quality"][trim_len:]

                    # New R1 with UMI tag
                    r1 = SeqRecord(Seq(trimmed_seq),
                                   id=f"{r1.id}_UM:{umi}",
                                   description="",
                                   letter_annotations={"phred_quality": trimmed_qual})

                    # New R2 with UMI tag
                    r2.id += f"_UM:{umi}"
                    r2.description = ""
                else:
                    # invalid UMI (has N), skip tagging
                    pass
            else:
                # tag not recognized — do nothing
                pass

            SeqIO.write(r1, out1, "fastq")
            SeqIO.write(r2, out2, "fastq")

            if i % 1000000 == 0:
                print(f"Processed {i:,} reads...", flush=True)

    print("Done.")

if __name__ == "__main__":
    main()
