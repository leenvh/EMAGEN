import pysam
from tqdm import tqdm
from pathogenprofiler.gff import load_gff, Exon
from pathogenprofiler.models import GenomePosition
from typing import List
from collections import defaultdict, Counter
import argparse

#Example:
#python scripts/realign-crt-haplotype.py --vcf genotyped.ann.vcf.gz --ref Pfalciparum.genome.fasta --gff Pfalciparum.genome.modified.new.gff3 --out genotyped_v2.vcf 

parser = argparse.ArgumentParser()
parser.add_argument('--vcf', required=True, help='VCF file')
parser.add_argument('--ref', required=True, help='Fasta file')
parser.add_argument('--gff', required=True, help='GFF file')
parser.add_argument('--out', required=True, help='New vcf file')
parser.add_argument('--region', default='Pf3D7_07_v3:403615-403635', help='Region')
args = parser.parse_args()

def get_overlapping_exons(chrom: str,pos: int,exons: List[Exon]):
    exons = [e for e in exons if e.start<=pos and e.end>=pos and chrom==e.chrom]
    if len(exons)==0:
        return None
    else:
        return exons[0]
    
def get_codon_pos(chrom: str,pos: int,exons: List[Exon]):
    e = get_overlapping_exons(chrom,pos,exons)
    if e==None:
        return (None,None)

    if e.strand=="+":
        codon_pos = (pos-e.start-e.phase)//3 + 1
    else:
        codon_pos = (e.end + e.phase - pos )//3 + 1
    return (e.id,codon_pos)
        


gff = load_gff(args.gff)
exons = []
for gene in gff:
    for transcript in gene.transcripts:
        for i,exon in enumerate(transcript.exons):
            exon.id = f"{gene.gene_id}_exon{i+1}"
            exons.append(exon)


vcf = pysam.VariantFile(args.vcf)



data = []
for sample in vcf.header.samples:
    variants = []
    for var in vcf.fetch('Pf3D7_07_v3',403610,403630):
        gt = var.samples[sample]['GT']
        alleles = [var.ref] + list(var.alts)
        if gt[0]==None or gt[1]==None:
            allele = None
        elif gt[0]!=gt[1]:
            allele = None
        else:
            allele = alleles[gt[0]]
        variants.append((var.pos,var.ref,allele))

    data.append((sample,tuple(variants)))

print('Variants extracted:',Counter([(x[1]) for x in data]))

def realign_variants(input_variants, ref_seq, seq_start):
    if any([x[2] is None for x in input_variants]):
        return None
    alt_seq = "".join(list(ref_seq))
    input_variants = sorted(input_variants, key=lambda x: x[0], reverse=True)
    for start, ref, alt in input_variants:
        alt_seq = alt_seq[:start-seq_start] + alt + alt_seq[start-seq_start+len(ref):]

    new_variants = []

    for i,(ref,alt) in enumerate(zip(ref_seq,alt_seq)):
        if ref != alt:
            new_variants.append((i+seq_start,ref,alt))
    
    coding_variants = defaultdict(list)
    for pos,ref,alt in new_variants:
        gene,cpos = get_codon_pos(chrom,pos,exons)
        coding_variants[(gene,cpos)].append((pos,ref,alt))
    
    final_variants = []
    for variants in coding_variants.values():
        if len(variants)==1:
            final_variants.append(variants[0])
        else:
            all_positions = [GenomePosition(chrom=chrom,pos=v[0]) for v in variants]
            positions = sorted(list(set(all_positions)))

            if len(positions)==2 and positions[0].pos == positions[1].pos-2:
                positions.insert(1,GenomePosition(chrom=positions[0].chrom,pos=positions[0].pos+1))

            ref_hap = {p.pos:ref_fasta.fetch(p.chrom,p.pos-1,p.pos).upper() for p in positions}
            alt_hap = {p.pos:ref_fasta.fetch(p.chrom,p.pos-1,p.pos).upper() for p in positions}
            for v in variants:
                alt_hap[v[0]] = v[2]
            ref_hap = ''.join(ref_hap.values())
            alt_hap = ''.join(alt_hap.values())

            final_variants.append((positions[0].pos,ref_hap,alt_hap))

    return final_variants


chrom,start_end = args.region.split(':')
start,end = start_end.split('-')
start = int(start)
end = int(end)

ref_fasta = pysam.FastaFile(args.ref)
ref_seq = ref_fasta.fetch(chrom, start-1, end).upper()
print('Reference sequence:',ref_seq)



realigned_variants = []
for d in data:
    realigned_variants.append((d[0],realign_variants(d[1], ref_seq, start)))

variants_found = set()
for s,variants in realigned_variants:
    if variants:
        for v in variants:
            variants_found.add(v)

print('New variants:',variants_found)
new_vcf = pysam.VariantFile(args.out,'w',header=vcf.header)


vcf = pysam.VariantFile(args.vcf)
print("Writing new VCF")
print("Writing variants before region")
for var in tqdm(vcf):
    if var.chrom=='Pf3D7_07_v3' and var.pos>=start:
        break
    new_vcf.write(var)

print("Writing new variants")
vcf_lines = []
for var in sorted(list(variants_found), key=lambda x: x[0]):
    new_vcf_variant = new_vcf.new_record()
    new_vcf_variant.chrom = 'Pf3D7_07_v3'
    new_vcf_variant.pos = var[0]
    new_vcf_variant.ref = var[1]
    new_vcf_variant.alts = [var[2]]
    new_vcf_variant.id = '.'
    new_vcf_variant.qual = 0
    new_vcf_variant.filter.add('PASS')
    for s,variants in realigned_variants:
        if variants is None:
            new_vcf_variant.samples[s]['GT'] = (None,None)
        elif var in variants:
            new_vcf_variant.samples[s]['GT'] = (1,1)
        else:
            new_vcf_variant.samples[s]['GT'] = (0,0)
    new_vcf.write(new_vcf_variant)

print("Writing variants after region")
for var in tqdm(vcf):
    if var.chrom<='Pf3D7_07_v3' and var.pos<=end:
        continue
    new_vcf.write(var)
new_vcf.close()