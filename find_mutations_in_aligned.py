#!/usr/bin/env python3
"""
find_mutations_in_aligned_columns.py

Usage:
  python3 find_mutations_in_aligned_columns.py --alignment dhfr_codon_aln.fasta --mutations hotspots.txt --out report.tsv

Input:
 - alignment: codon alignment FASTA (with '-' gaps)
 - mutations: list like C466T (where 466 = alignment column number, 1-based)
Output:
 - TSV file showing which sequences carry each mutation, codon, AA change, synonymous?
"""
import argparse, math
from Bio import SeqIO
from Bio.Seq import Seq
import re

p = argparse.ArgumentParser()
p.add_argument("--alignment", required=True)
p.add_argument("--mutations", required=True)
p.add_argument("--out", default="mutation_report.tsv")
args = p.parse_args()

# Load alignment
seqs = list(SeqIO.parse(args.alignment, "fasta"))
if not seqs:
    raise SystemExit("No sequences loaded from alignment.")

aligned_length = len(seqs[0].seq)
print(f"Loaded {len(seqs)} sequences, alignment length = {aligned_length} columns.")

# Load mutation list
mutlist = []
pat = re.compile(r'([ACGT])(\d+)([ACGT])', re.I)
with open(args.mutations) as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        m = pat.match(line)
        if not m:
            print(f"Skipping invalid mutation format: {line}")
            continue
        ref = m.group(1).upper()
        col = int(m.group(2))  # alignment column number
        alt = m.group(3).upper()
        mutlist.append((ref, col, alt))

# Write output
with open(args.out, "w") as outf:
    outf.write("mutation\tseq_id\tbase_at_column\tmatch_mutant\tcodon_index_in_alignment\taligned_codon\tmutated_codon\tAA_before\tAA_after\tsynonymous\n")
    for ref, col, alt in mutlist:
        # Convert to 0-based for Python indexing
        idx = col - 1
        for rec in seqs:
            seq = str(rec.seq).upper()
            if idx >= len(seq):
                base = ""
                match = "NA"
                codon_idx = aligned_codon = mutated_codon = AA_before = AA_after = syn = ""
            else:
                base = seq[idx]
                match = "YES" if base == alt else "NO"
                # Figure codon index and codon
                codon_idx = math.ceil((idx+1)/3)  # 1-based codon number in alignment
                codon_start = (codon_idx-1)*3
                if codon_start+3 <= len(seq):
                    aligned_codon = seq[codon_start:codon_start+3]
                    # mutate codon if base not a gap
                    codon_list = list(aligned_codon)
                    pos_in_codon = (idx) % 3
                    codon_list[pos_in_codon] = alt
                    mutated_codon = "".join(codon_list)
                    # Remove gaps for translation
                    if "-" in aligned_codon or "-" in mutated_codon:
                        AA_before = AA_after = syn = ""
                    else:
                        AA_before = str(Seq(aligned_codon).translate())
                        AA_after = str(Seq(mutated_codon).translate())
                        syn = "YES" if AA_before == AA_after else "NO"
                else:
                    aligned_codon = mutated_codon = AA_before = AA_after = syn = ""
            outf.write(f"{ref}{col}{alt}\t{rec.id}\t{base}\t{match}\t{codon_idx}\t{aligned_codon}\t{mutated_codon}\t{AA_before}\t{AA_after}\t{syn}\n")

print(f"Report written to {args.out}")
