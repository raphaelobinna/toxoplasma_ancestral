#!/usr/bin/env python3
"""
Convert BLAST XML output to FASTA format
Extracts hit sequences from BLAST XML results and saves them as FASTA
"""

import xml.etree.ElementTree as ET
import argparse
import sys
from pathlib import Path

def parse_blast_xml_to_fasta(xml_file, output_file, max_hits=None, evalue_threshold=None):
    """
    Parse BLAST XML file and convert to FASTA format
    
    Args:
        xml_file (str): Path to BLAST XML file
        output_file (str): Path to output FASTA file
        max_hits (int, optional): Maximum number of hits per query
        evalue_threshold (float, optional): E-value threshold filter
    """
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
    except ET.ParseError as e:
        print(f"Error parsing XML file: {e}")
        return False
    except FileNotFoundError:
        print(f"File not found: {xml_file}")
        return False

    sequences_written = 0
    
    with open(output_file, 'w') as fasta_out:
        # Iterate through iterations (queries)
        for iteration in root.findall('.//Iteration'):
            query_def = iteration.find('Iteration_query-def')
            query_id = query_def.text if query_def is not None else "Unknown_Query"
            
            hits = iteration.findall('.//Hit')
            hit_count = 0
            
            # Process each hit
            for hit in hits:
                if max_hits and hit_count >= max_hits:
                    break
                    
                hit_id = hit.find('Hit_id')
                hit_def = hit.find('Hit_def')
                hit_accession = hit.find('Hit_accession')
                
                # Get the best HSP for this hit
                hsps = hit.findall('.//Hsp')
                if not hsps:
                    continue
                    
                best_hsp = hsps[0]  # Usually sorted by score
                
                # Check E-value threshold if specified
                if evalue_threshold:
                    evalue_elem = best_hsp.find('Hsp_evalue')
                    if evalue_elem is not None:
                        try:
                            evalue = float(evalue_elem.text)
                            if evalue > evalue_threshold:
                                continue
                        except ValueError:
                            continue
                
                # Extract sequence information
                hsp_hseq = best_hsp.find('Hsp_hseq')
                if hsp_hseq is None:
                    continue
                    
                sequence = hsp_hseq.text.replace('-', '')  # Remove gaps
                
                # Create FASTA header
                header_parts = []
                if hit_id is not None:
                    header_parts.append(hit_id.text)
                if hit_accession is not None and hit_accession.text != hit_id.text:
                    header_parts.append(f"acc={hit_accession.text}")
                if hit_def is not None:
                    header_parts.append(hit_def.text)
                
                header = ' '.join(header_parts) if header_parts else f"hit_{sequences_written + 1}"
                
                # Write FASTA entry
                fasta_out.write(f">{header}\n")
                
                # Write sequence in 80-character lines
                for i in range(0, len(sequence), 80):
                    fasta_out.write(sequence[i:i+80] + '\n')
                
                sequences_written += 1
                hit_count += 1

    print(f"Successfully converted {sequences_written} sequences to {output_file}")
    return True

def parse_blast_xml_query_sequences(xml_file, output_file):
    """
    Extract query sequences from BLAST XML file
    
    Args:
        xml_file (str): Path to BLAST XML file
        output_file (str): Path to output FASTA file
    """
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
    except ET.ParseError as e:
        print(f"Error parsing XML file: {e}")
        return False
    except FileNotFoundError:
        print(f"File not found: {xml_file}")
        return False

    sequences_written = 0
    
    with open(output_file, 'w') as fasta_out:
        # Extract query sequences
        for iteration in root.findall('.//Iteration'):
            query_def = iteration.find('Iteration_query-def')
            query_seq = iteration.find('Iteration_query-seq')
            
            if query_def is not None and query_seq is not None:
                header = query_def.text
                sequence = query_seq.text
                
                # Write FASTA entry
                fasta_out.write(f">{header}\n")
                
                # Write sequence in 80-character lines
                for i in range(0, len(sequence), 80):
                    fasta_out.write(sequence[i:i+80] + '\n')
                
                sequences_written += 1

    print(f"Successfully extracted {sequences_written} query sequences to {output_file}")
    return True

def main():
    parser = argparse.ArgumentParser(description='Convert BLAST XML to FASTA format')
    parser.add_argument('input_xml', help='Input BLAST XML file')
    parser.add_argument('output_fasta', help='Output FASTA file')
    parser.add_argument('--mode', choices=['hits', 'queries'], default='hits',
                      help='Extract hit sequences (default) or query sequences')
    parser.add_argument('--max-hits', type=int, 
                      help='Maximum number of hits per query (hits mode only)')
    parser.add_argument('--evalue-threshold', type=float,
                      help='E-value threshold for filtering hits (hits mode only)')
    
    args = parser.parse_args()
    
    # Validate input file
    if not Path(args.input_xml).exists():
        print(f"Error: Input file '{args.input_xml}' does not exist")
        sys.exit(1)
    
    # Convert based on mode
    if args.mode == 'hits':
        success = parse_blast_xml_to_fasta(
            args.input_xml, 
            args.output_fasta, 
            args.max_hits, 
            args.evalue_threshold
        )
    else:  # queries
        success = parse_blast_xml_query_sequences(args.input_xml, args.output_fasta)
    
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()

# Example usage:
# python blast_xml_to_fasta.py blast_results.xml output_sequences.fasta
# python blast_xml_to_fasta.py blast_results.xml output_sequences.fasta --mode hits --max-hits 10 --evalue-threshold 1e-5
# python blast_xml_to_fasta.py blast_results.xml query_sequences.fasta --mode queries