from __future__ import print_function
import os
import os.path
import sys
import gzip
import itertools
import json
import multiprocessing
import glob
import errno
import csv
from datetime import datetime

## Parameters
r1_length = 16
main_path = "/mnt/d/PostdocUnsync/Help/Morgane/dgeSeq/data202306"
sampleSheet = "/mnt/d/PostdocUnsync/Help/Morgane/dgeSeq/data202306/metadata/samplesheet.tsv"
fastq_folder = "/mnt/d/PostdocUnsync/Help/Morgane/dgeSeq/data202306/FASTQ"


fastq_pairs= [
    {
        "read1": "/mnt/d/PostdocUnsync/Help/Morgane/dgeSeq/data202306/multiFASTQ/A001_Frapin_3-SRP_1_H710_CGAGGCTG_run20230601R_S1_R1_001.fastq.gz",
        "read2": "/mnt/d/PostdocUnsync/Help/Morgane/dgeSeq/data202306/multiFASTQ/A001_Frapin_3-SRP_1_H710_CGAGGCTG_run20230601R_S1_R2_001.fastq.gz"
    }
]

def read_tsv_as_dict(filename):
    data = []
    with open(filename, 'r') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            data.append(dict(row))
    return data

# Open .fastq or .fastq.gz files for reading
def open_fastq_or_gz(filename):
    if filename.endswith(".fastq") and os.access(filename, os.F_OK):
        return open(filename, "rU")
    elif filename.endswith(".fastq.gz") and os.access(filename, os.F_OK):
        return gzip.open(filename, "rt")
    elif filename.endswith(".fastq") and os.access(filename + ".gz", os.F_OK):
        return gzip.open(filename + ".gz", "rt")
    elif filename.endswith(".fastq.gz") and os.access(filename[:-3], os.F_OK):
        return open(filename[:-3], "rU")
    raise IOError("Unknown file: " + filename)

def createEmptyFastq(index2sample):
    for index, sample in index2sample.items():
        with open(sample["fastq_path"], "w") as out:
            pass
        out.close

# Write fastq file for each sample
def write_fastq(alignment_dir, sample_well, sample_name, sample_index, sample_name_seq_qual_list):
    fastq_path = sample["fastq_path"]
    if not os.path.exists(os.path.dirname(fastq_path)):
        try:
            os.makedirs(os.path.dirname(fastq_path))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    with open(fastq_path, "a") as out:
        for name, seq, qual in sample_name_seq_qual_list:
            print("\n".join(["@" + name, seq, "+", qual]), file=out)
    out.close

# Mask sequence by quality score
def mask(seq, qual, min_qual=10):
    return "".join((b if (ord(q) - 33) >= min_qual else "N") for b, q in itertools.zip_longest(seq, qual))

# Read fastq and split reads according to barcode


def dischargeBuff(index2sample):
    print(datetime.now())
    for index, sample in index2sample.items():
        with open(sample["fastq_path"], "a") as out:
            for name, seq, qual in sample["buf"]:
                print("\n".join(["@" + name, seq, "+", qual]), file=out)
        out.close

        sample["buf"] = list()


## Main
samples = read_tsv_as_dict(sampleSheet)

fastq_pair = fastq_pairs[0]

index2sample = dict()
for sample in samples:
    sample["buf"] = list()
    sample["index"] = sample["index"][0:5]
    sample["fastq_path"] = os.path.join(fastq_folder, ".".join([sample["name"],"fastq"]))
    index2sample[sample["index"]] = sample

possibleIndex = index2sample.keys()


index2sample["undefined"] = {"name":"sample_unknown","index":"index_unknown","well":"well_unknown","fastq_path":os.path.join(fastq_folder, ".".join(["unknown","fastq"]))}
index2sample["undefined"]["buf"] = list()

createEmptyFastq(index2sample)

readF = fastq_pair["read1"]
readR = fastq_pair["read2"]

with open_fastq_or_gz(readF) as r1_file, open_fastq_or_gz(readR) as r2_file:
    r1_r2 = zip(r1_file, r2_file)
    dischargeBuffindent = 0
    for header1, header2 in r1_r2:
        
        seq1, seq2 = next(r1_r2)
        plus1, plus2 = next(r1_r2)
        qual1, qual2 = next(r1_r2)

        read_name1, read_name2 = header1.split()[0][1:], header2.split()[0][1:]
        assert read_name1 == read_name2
        seq2, qual2 = seq2.rstrip(), qual2.rstrip()
        index, umi, seq, qual = mask(seq1[0:5], qual1[0:5], min_qual=10), mask(seq1[6:r1_length], qual1[6:r1_length]), seq2, qual2   
        barcoded_name = "_".join([read_name2, umi])

        if index not in possibleIndex:
            index = "undefined"

        index2sample[index]["buf"].append((barcoded_name, seq, qual))
        dischargeBuffindent += 1
        if dischargeBuffindent == 500000:
            print(datetime.now())
            dischargeBuff(index2sample)
            dischargeBuffindent = 0

dischargeBuff(index2sample)

