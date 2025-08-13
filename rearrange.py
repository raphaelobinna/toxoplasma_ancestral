#!/usr/bin/env python3
"""
Align protein and mRNA FASTA files by gene ID
Reads protein sequences, extracts gene IDs, finds matching mRNA sequences,
and creates aligned output files where indexes correspond.
"""

import re
import argparse
import sys
from pathlib import Path

def accumulate_header(filename):
    """
    Parse a FASTA file and return a dictionary of sequences
    
    Args:
        filename (str): Path to FASTA file
        
    Returns:
        dict: Dictionary with headers as keys and sequences as values
    """
    sequences = {}
    current_header = None
    all_headers = []
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    
                    
                    current_header = line[1:]  # Remove '>' character
                    all_headers.append(current_header)
                
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found")
        return None
    except Exception as e:
        print(f"Error reading headers '{filename}': {e}")
        return None
    
    return all_headers

def parse_fasta(filename):
    """
    Parse a FASTA file and return a dictionary of sequences
    
    Args:
        filename (str): Path to FASTA file
        
    Returns:
        dict: Dictionary with headers as keys and sequences as values
    """
    sequences = {}
    current_header = None
    current_seq = []
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_header:
                        sequences[current_header] = ''.join(current_seq)
                    
                    current_header = line[1:]  # Remove '>' character
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Save last sequence
            if current_header:
                sequences[current_header] = ''.join(current_seq)
                
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found")
        return None
    except Exception as e:
        print(f"Error reading file '{filename}': {e}")
        return None
    
    return sequences

def extract_gene_id_from_protein_header(header):
    """
    Extract gene ID from protein FASTA header
    
    Args:
        header (str): Protein sequence header
        
    Returns:
        str: Gene ID or None if not found
    """
    # Look for gene=GENE_ID pattern
    gene_match = re.search(r'gene=([^\s|]+)', header)
    if gene_match:
        gene_id = gene_match.group(1)
        # Convert TGP89_249180 format to TGARI_249180 format if needed
        # This handles the case where protein uses TGP89 and mRNA uses TGARI

        return gene_id
    return None

def find_matching_mrna_gene(gene_id, mrna_headers):
    """
    Find matching mRNA header for a given gene ID
    
    Args:
        gene_id (str): Gene ID to search for
        mrna_headers (list): List of mRNA headers
        
    Returns:
        str: Matching mRNA header or None
    """
    # Direct match
    for header in mrna_headers:
        if gene_id in header:
            return header
    
    return None

def align_sequences(protein_file, mrna_file, output_mrna=None):
    """
    Align protein and mRNA sequences by gene ID
    
    Args:
        protein_file (str): Path to protein FASTA file
        mrna_file (str): Path to mRNA FASTA file
        output_protein (str): Output aligned protein FASTA file
        output_mrna (str): Output aligned mRNA FASTA file
    """
    print("Reading protein sequences...")
    protein_sequences = accumulate_header(protein_file)
    if protein_sequences is None:
        return False
    
    print("Reading mRNA sequences...")
    mrna_sequences = parse_fasta(mrna_file)
    if mrna_sequences is None:
        return False
    
    print(f"Found {len(protein_sequences)} protein sequences")
    print(f"Found {len(mrna_sequences)} mRNA sequences")
    
    # Extract gene IDs from protein sequences

    gene_order = []
    
    for header in protein_sequences:
        gene_id = extract_gene_id_from_protein_header(header)
        if gene_id:
            
            gene_order.append(gene_id)
            print(f"Extracted gene ID: {gene_id}")
        else:
            print(f"Warning: Could not extract gene ID from header: {header}")
    
    # Find matching mRNA sequences
    matched_pairs = []
    unmatched_genes = []
    mrna_headers_list = list(mrna_sequences.keys())
    
    for gene_id in gene_order:
       
        matching_mrna_header = find_matching_mrna_gene(gene_id, mrna_headers_list)
        
        if matching_mrna_header:
            matched_pairs.append(( matching_mrna_header, gene_id))
            print(f"Matched: {gene_id} -> {matching_mrna_header}")
        else:
            unmatched_genes.append(gene_id)
            print(f"Warning: No mRNA sequence found for gene: {gene_id}")
    
    print(f"\nMatched {len(matched_pairs)} sequences")
    if unmatched_genes:
        print(f"Unmatched genes: {len(unmatched_genes)}")
        for gene in unmatched_genes[:5]:  # Show first 5 unmatched
            print(f"  - {gene}")
        if len(unmatched_genes) > 5:
            print(f"  ... and {len(unmatched_genes) - 5} more")
    

    
    if output_mrna:
        print(f"Writing aligned mRNA sequences to {output_mrna}")
        with open(output_mrna, 'w') as f:
            for  mrna_header, gene_id in matched_pairs:
                f.write(f">{mrna_header}\n")
                sequence = mrna_sequences[mrna_header]
                # Write sequence in 80-character lines
                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + '\n')
    
    # Create summary report
    print(f"\n=== ALIGNMENT SUMMARY ===")
    print(f"Total protein sequences: {len(protein_sequences)}")
    print(f"Total mRNA sequences: {len(mrna_sequences)}")
    print(f"Successfully matched: {len(matched_pairs)}")
    print(f"Unmatched proteins: {len(unmatched_genes)}")
    
    return True

def main():
    parser = argparse.ArgumentParser(description='Align protein and mRNA FASTA files by gene ID')
    parser.add_argument('protein_fasta', help='Input protein FASTA file')
    parser.add_argument('mrna_fasta', help='Input mRNA FASTA file')

    parser.add_argument('--output-mrna', help='Output aligned mRNA FASTA file (optional)')
    
    args = parser.parse_args()
    
    # Validate input files
    if not Path(args.protein_fasta).exists():
        print(f"Error: Protein file '{args.protein_fasta}' does not exist")
        sys.exit(1)
    
    if not Path(args.mrna_fasta).exists():
        print(f"Error: mRNA file '{args.mrna_fasta}' does not exist")
        sys.exit(1)
    
    # Set default output filenames if not provided

    output_mrna = args.output_mrna or 'aligned_mrna_sequences.fasta'
    
    success = align_sequences(args.protein_fasta, args.mrna_fasta, output_mrna)
    
    if not success:
        sys.exit(1)
    
    print(f"\nAlignment complete! Check output files:")

    print(f"  - {output_mrna}")

if __name__ == "__main__":
    main()

# Example usage:
# python align_sequences.py protein_sequences.fasta mrna_sequences.fasta
# python align_sequences.py protein_sequences.fasta mrna_sequences.fasta --output-protein aligned_proteins.fasta --output-mrna aligned_mrna.fasta