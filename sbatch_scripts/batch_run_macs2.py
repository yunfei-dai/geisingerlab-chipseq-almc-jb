## Yunfei Dai
## 20Sep2023

import os, re, sys

INFILE_DIR = sys.argv[1]
OUTPUT_DIR = sys.argv[2]

for file in os.listdir(INFILE_DIR):
    if re.match(r"^.+BEADS.+.bam$", file):
        exp_name = re.search(r"(^.+_BEADS)_(\d)_.+.bam$", file)
        experiment = exp_name.group(1)
        replicate = exp_name.group(2)
        # control is JBA71_INPUT
        infile = os.path.join(INFILE_DIR, file)
        control = os.path.join(INFILE_DIR, "JBA71_INPUT_"+replicate+"_S"+replicate+"_L001_R1_001.sorted.bam")
        outname = experiment + "_" + replicate + ".ext_size200"
        outname = os.path.join(OUTPUT_DIR, outname)
        cmd = "macs2 callpeak -t " + \
                infile + " -c " + control + \
                " -g 3.8e6 -n " + outname + \
                " --nomodel --extsize 200"
        os.system(cmd)
        print("Processing: " + infile)
