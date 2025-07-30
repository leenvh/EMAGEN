import subprocess as sp
import argparse
import os
from uuid import uuid4
import shutil

argparser = argparse.ArgumentParser(description="Script to generate a bed file with the depth of coverage in bins")
argparser.add_argument("-i","--input", type=str, default="-",help="Comma separated list of values")
argparser.add_argument("-b","--bed", type=str, default="-",help="Comma separated list of values")
argparser.add_argument("-o","--output", type=str, default="-",help="Comma separated list of values")

args = argparser.parse_args()

random_id = str(uuid4())
os.mkdir(random_id)

for i,l in enumerate(sp.Popen(f'bedtools coverage -a {args.bed} -b {args.input}',shell=True,stdout=sp.PIPE).stdout):
    row = l.decode().strip().split()
    total_reads = int(row[5])*2
    if int(row[4]) == 0:
        factor = 1
    else:
        factor = total_reads/int(row[4])
    if factor>1:
        factor = 1
    print(f"Downsampling {row[0]}:{row[1]}-{row[2]} to {factor}")
    sp.run(f"samtools view {args.input} {row[0]}:{row[1]}-{row[2]} -b | samtools view --subsample {factor} -b > {random_id}/region{i}.bam",shell=True)

sp.run(f'samtools merge -f -o {args.output} {random_id}/*.bam',shell=True)
sp.run(f"samtools index {args.output}",shell=True)

shutil.rmtree(random_id)

    
