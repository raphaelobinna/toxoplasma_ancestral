import pandas as pd
import numpy as np
from Bio import Phylo, SeqIO

def read_fasta(file_path):
    record = next(SeqIO.parse(file_path, "fasta"))  # Get first sequence only
    seq_str = str(record.seq)
    print(seq_str[:15])  # Print first 15 characters
    return seq_str

def map_mutation_seq(pdb_seq, mutation_file):
    mut_file = pd.read_csv(mutation_file, sep="\t")

    for _, row in mut_file.iterrows():
        if row['match_mutant'] == "NO":
            continue
        
        codon_idx = int(row['codon_index_in_alignment'])
        pdb_pos = codon_idx

        if pdb_pos is None:
            print(f"Codon {codon_idx} is a gap in PDB sequence")
            continue

        ref_aa = pdb_seq[pdb_pos - 1]  # Convert to 0-based
        if ref_aa == row['AA_before']:
            print(f"✓ Match: Mutation {row['mutation']} at PDB AA {pdb_pos}: Ref {ref_aa} → {row['AA_after']}")
        else:
            print(f"⚠ Mismatch: PDB AA {pdb_pos} = {ref_aa}, table says {row['AA_before']}")




def main():
    try:
        protein_fasta ="rcsb_pdb_4ECK (1).fasta"
        mutation_file = "mutation_report.tsv"

        fasta_seq = read_fasta(protein_fasta)

        final = map_mutation_seq(fasta_seq, mutation_file)
    except FileNotFoundError as e:
        print(f"File not found: {e}")
        print("Please check your file paths and make sure all files exist")
    except Exception as e:
        print(f"Error during analysis: {e}")
        return None

if __name__ == "__main__":
    result = main()



# ⚠ Mismatch: PDB AA 591 = E, table says D
# ⚠ Mismatch: PDB AA 591 = E, table says D
# ⚠ Mismatch: PDB AA 591 = E, table says D

