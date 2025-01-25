#!/usr/bin/env python3

"""
AIOalign5.6.py - Automated DNA sequence alignment with ORF optimization

This program performs a DNA sequence alignment by:
1. Finding all possible open reading frames (ORFs) in input sequences
2. Trimming untranslated regions (UTRs) for each ORF
3. Translating sequences to amino acids
4. Performing initial protein alignment
5. Selecting optimal ORFs based on alignment quality
6. Running final alignment on selected sequences
7. Back-translating to DNA while maintaining alignment

The process preserves codon degeneracy and silent mutation information while
leveraging the advantages of amino acid alignment for proper reading frame maintenance.
"""

import argparse
import logging
import multiprocessing as mp
import os
import subprocess
import sys
import re
from typing import List, Tuple, Dict
from pathlib import Path

from Bio import SeqIO, AlignIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class ORFFinder:
    """
    Identifies potential open reading frames (ORFs) in DNA sequences.
    Handles both forward and reverse complement strands.
    """
    
    @staticmethod
    def find_orfs(seq: Seq, min_protein_length: int = 50) -> List[Tuple[int, int, int]]:
        """
        Finds all possible ORFs in a DNA sequence that meet minimum length requirements.
        
        Args:
            seq (Seq): DNA sequence to analyze
            min_protein_length (int): Minimum protein length (in amino acids) to consider
                                    Default: 50 amino acids (150 nucleotides)
        
        Returns:
            List[Tuple[int, int, int]]: List of ORFs as (start, end, strand) tuples
                                       strand is 1 for forward, -1 for reverse complement
        """
        # Remove alignment gaps before processing
        clean_seq = str(seq).replace('-', '')
        seq = Seq(clean_seq)
        seq_len = len(seq)
        orfs = []
        
        # Check both forward and reverse complement strands
        for frame in range(3):  # Each possible reading frame (0, 1, 2)
            for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
                for start in range(frame, seq_len - 2, 3):  # Step by codons
                    # Look for start codon
                    if nuc[start:start+3] == 'ATG':
                        # Search for stop codon
                        for end in range(start + 3, seq_len - 2, 3):
                            codon = nuc[end:end+3]
                            if codon in ['TAA', 'TAG', 'TGA']:  # Stop codons
                                if end - start >= min_protein_length * 3:
                                    # Convert coordinates for reverse complement
                                    if strand == 1:
                                        orf = (start, end + 3, strand)
                                    else:
                                        orf = (seq_len - end - 3, seq_len - start, strand)
                                    orfs.append(orf)
                                break
        return orfs

class AlignmentProcessor:
    """
    Handles sequence alignment operations and scoring.
    Includes BLOSUM62 matrix for amino acid similarity scoring.
    """
    
    # BLOSUM62 scoring matrix for amino acid similarity calculations
    # Each tuple key represents a pair of amino acids and their corresponding score
    BLOSUM62 = {
        ('A', 'A'):  4, ('R', 'A'): -1, ('N', 'A'): -2, ('D', 'A'): -2, ('C', 'A'):  0,
        # ... [rest of BLOSUM62 matrix remains unchanged]
    }

    @staticmethod
    def evaluate_start_region(ref_seq: str, test_seq: str, region_length: int = 30) -> float:
        """
        Evaluates the quality of sequence start regions using progressive weighting.
        
        Args:
            ref_seq (str): Reference sequence
            test_seq (str): Test sequence to compare
            region_length (int): Length of region to evaluate (default: 30)
            
        Returns:
            float: Score multiplier between 0.1 and 4.0 based on start region quality
        """
        if len(test_seq) < region_length or len(ref_seq) < region_length:
            return 0.1
            
        ref_start = ref_seq[:region_length]
        test_start = test_seq[:region_length]
        
        # Define sections for granular analysis with different weights
        sections = [
            (0, 20),    # First 20 AA - most critical
            (20, 30),   # Next 10 AA - very important
            (30, 40)    # Last 10 AA - important
        ]
        
        section_weights = [4.0, 2.0, 1.5]  # Progressive weighting
        total_score = 0
        
        # Analyze each section separately
        for (start, end), weight in zip(sections, section_weights):
            ref_section = ref_start[start:end]
            test_section = test_start[start:end]
            
            # Calculate sequence identity in this section
            matches = sum(1 for a, b in zip(ref_section, test_section) if a == b)
            identity = matches / (end - start)
            
            # Penalize gaps
            gaps = test_section.count('-')
            gap_penalty = max(0.1, 1.0 - (gaps * 0.3))
            
            # Calculate weighted section score
            section_score = identity * gap_penalty * weight
            total_score += section_score
        
        # Normalize and scale final score
        final_multiplier = (total_score / sum(section_weights)) * 4.0
        return max(0.1, min(4.0, final_multiplier))

    @staticmethod
    def score_sequence_similarity(seq1: str, seq2: str, ref_length: int) -> float:
        """
        Scores the similarity between two sequences with emphasis on N-terminal region.
        Uses regional weighting and gap penalties for more accurate ORF selection.
        
        Args:
            seq1 (str): Reference sequence
            seq2 (str): Sequence to compare
            ref_length (int): Length of reference sequence
            
        Returns:
            float: Similarity score weighted by region importance
        """
        # Clean sequences of gaps
        clean_seq1 = seq1.replace('-', '')
        clean_seq2 = seq2.replace('-', '')
        
        # Get start region quality multiplier
        start_quality = AlignmentProcessor.evaluate_start_region(clean_seq1, clean_seq2)
        
        # Configure local alignment with strict gap penalties
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = 2
        aligner.mismatch_score = -3
        aligner.open_gap_score = -65
        aligner.extend_gap_score = -40
        
        alignment = aligner.align(clean_seq1, clean_seq2)[0]
        
        # Define regional weights for scoring
        n_term_weight = 8.0     # Extreme weight for first region
        early_weight = 4.0      # Very high weight for early region
        mid_weight = 1.0        # Normal weight for middle
        c_term_weight = 2.0     # Moderate weight for C-terminal
        
        # Define region boundaries
        aligned_length = len(alignment[0])
        n_term_region = 30      # First 30 amino acids
        early_region = 50       # First 50 amino acids
        c_term_region = int(aligned_length * 0.8)
        
        # Define position-specific gap penalties
        gap_penalties = {
            'n_term': -10,      # Severe penalty for N-terminal gaps
            'early': -5,        # High penalty for early region gaps
            'mid': -2,          # Normal penalty for middle region
            'c_term': -3        # Moderate penalty for C-terminal gaps
        }
    
        # Calculate alignment score with regional weighting
        raw_score = 0
        for i, (a, b) in enumerate(zip(alignment[0], alignment[1])):
            if a == '-' or b == '-':
                # Apply position-specific gap penalties
                if i < n_term_region:
                    raw_score += gap_penalties['n_term']
                elif i < early_region:
                    raw_score += gap_penalties['early']
                elif i > c_term_region:
                    raw_score += gap_penalties['c_term']
                else:
                    raw_score += gap_penalties['mid']
                continue
            
            # Get BLOSUM62 score for amino acid pair
            pair = (a, b)
            score = (AlignmentProcessor.BLOSUM62.get(pair) or 
                    AlignmentProcessor.BLOSUM62.get((b, a), 0))
            
            # Apply regional weights with emphasis on N-terminal
            if i < n_term_region:
                weight = n_term_weight
                if a == b:  # Bonus for exact matches in N-terminal
                    score *= 4.0
                elif score < 0:  # Extra penalty for N-terminal mismatches
                    score *= 2
            elif i < early_region:
                weight = early_weight
                if a == b:  # Bonus for early region matches
                    score *= 2
            elif i > c_term_region:
                weight = c_term_weight
            else:
                weight = mid_weight
                
            raw_score += score * weight
        
        # Apply length penalty to discourage significantly different lengths
        length_ratio = len(clean_seq2) / ref_length
        length_penalty = 1 / (1 + abs(1 - length_ratio) ** 3.5)
        
        # Apply perfect start bonus if first 10 amino acids match exactly
        perfect_start_bonus = 2
        if clean_seq1[:10] == clean_seq2[:10]:
            perfect_start_bonus = 4.0
        
        # Calculate final score combining all factors
        final_score = raw_score * length_penalty * start_quality * perfect_start_bonus
        
        return final_score
        
    @staticmethod
    def run_clustalo(input_file: str, output_file: str):
        """
        Runs Clustal Omega multiple sequence alignment with optimized parameters.
        
        Args:
            input_file (str): Path to input sequences
            output_file (str): Path for aligned output
            
        Raises:
            subprocess.CalledProcessError: If Clustal Omega fails
            FileNotFoundError: If Clustal Omega is not installed
        """
        try:
            # Configure Clustal Omega command with optimized parameters
            cmd = [
                'clustalo',
                '-i', input_file,
                '-o', output_file,
                '--force',  # Overwrite existing files
                '--threads', str(mp.cpu_count() - 1),  # Use all but one CPU
                '--guidetree-out=guidetree.txt',  # Save guide tree for debugging
                '--full',  # Use full distance matrix for better accuracy
                '--iter', '2',  # Number of iterations
                '--clustering-out=clusters.txt'  # Save clustering information
            ]
            
            logging.info(f"Running Clustal Omega command: {' '.join(cmd)}")
            
            # Execute Clustal Omega
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            # Log output for debugging
            if result.stdout:
                logging.debug(f"Clustal Omega stdout: {result.stdout}")
            if result.stderr:
                logging.debug(f"Clustal Omega stderr: {result.stderr}")
                
        except subprocess.CalledProcessError as e:
            logging.error(f"Clustal Omega failed with error code {e.returncode}")
            logging.error(f"stdout: {e.stdout}")
            logging.error(f"stderr: {e.stderr}")
            raise
        except FileNotFoundError:
            logging.error("Clustal Omega not found. Please ensure it is installed and in PATH.")
            raise

class SequenceProcessor:
    """
    Handles sequence processing operations including ORF extraction,
    translation, and back-translation.
    """
    
    @staticmethod
    def extract_orfs(dna_file: str, output_file: str) -> Tuple[str, SeqRecord]:
        """
        Extracts all possible ORFs from input sequences.
        
        Args:
            dna_file (str): Path to input DNA sequences
            output_file (str): Path for extracted ORFs output
            
        Returns:
            Tuple[str, SeqRecord]: Path to output file and reference sequence
        """
        orf_records = []
        
        # Read all sequences from input file
        all_records = list(SeqIO.parse(dna_file, "fasta"))
        if not all_records:
            raise ValueError("No sequences found in input file")
        
        # Process reference sequence (first sequence in file)
        reference_record = all_records[0]
        reference_seq = str(reference_record.seq).replace('-', '')
        reference_record.seq = Seq(reference_seq)
        orf_records.append(reference_record)
        
        # Extract ORFs from remaining sequences
        for record in all_records[1:]:
            clean_seq = str(record.seq).replace('-', '')
            record.seq = Seq(clean_seq)
            
            # Find all possible ORFs
            orfs = ORFFinder.find_orfs(record.seq)
            for i, (start, end, strand) in enumerate(orfs):
                # Extract and store ORF sequence
                orf_seq = record.seq[start:end] if strand == 1 else record.seq[start:end].reverse_complement()
                new_id = f"{record.id}_ORF{i+1}"
                new_record = SeqRecord(
                    orf_seq,
                    id=new_id,
                    description=f"ORF {start}-{end} Strand {strand}"
                )
                orf_records.append(new_record)
        
        SeqIO.write(orf_records, output_file, "fasta")
        return output_file, reference_record

    @staticmethod
    def translate_sequences(sequences: List[SeqRecord]) -> List[SeqRecord]:
        """
        Translates DNA sequences to protein sequences.
        
        Args:
            sequences (List[SeqRecord]): DNA sequences to translate
            
        Returns:
            List[SeqRecord]: Translated protein sequences
        """
        translated = []
        for seq in sequences:
            # Remove gaps and ensure sequence length is divisible by 3
            clean_seq = str(seq.seq).replace('-', '')
            remainder = len(clean_seq) % 3
            if remainder:
                clean_seq = clean_seq[:-remainder]
            
            # Translate to protein sequence
            protein_seq = Seq(clean_seq).translate()
            translated.append(SeqRecord(
                protein_seq,
                id=seq.id,
                description="translated"
            ))
        return translated

    @staticmethod
    def back_translate(aligned_protein_file: str, original_dna: Dict[str, str], output_file: str):
        """
        Back-translates aligned protein sequences to DNA while maintaining alignment.
        
        Args:
            aligned_protein_file (str): Path to aligned protein sequences
            original_dna (Dict[str, str]): Original DNA sequences keyed by sequence ID
            output_file (str): Path for output aligned DNA sequences
        """
        aligned_dna = []
        
        for protein_record in SeqIO.parse(aligned_protein_file, "fasta"):
            if protein_record.id not in original_dna:
                continue
                
            dna_seq = original_dna[protein_record.id]
            aligned_dna_seq = ""
            dna_pos = 0
            
            # Convert each amino acid position back to corresponding codon
            for aa in protein_record.seq:
                if aa == '-':
                    aligned_dna_seq += '---'  # Maintain alignment gaps
                else:
                    if dna_pos + 3 <= len(dna_seq):
                        aligned_dna_seq += dna_seq[dna_pos:dna_pos + 3]
                        dna_pos += 3
            
            aligned_dna.append(SeqRecord(
                Seq(aligned_dna_seq),
                id=protein_record.id,
                description="back_translated"
            ))
        
        SeqIO.write(aligned_dna, output_file, "fasta")

def process_sequence_chunk(args: Tuple[SeqRecord, List[SeqRecord], int]) -> Tuple[str, SeqRecord]:
    """
    Processes a chunk of sequences to find the best ORF match to reference.
    
    Args:
        args (Tuple): Contains reference sequence, list of ORFs, and reference length
        
    Returns:
        Tuple[str, SeqRecord]: Original sequence ID and best matching ORF
    """
    reference_seq, orfs, ref_length = args
    best_score = float('-inf')
    best_orf = None
    ref_protein = str(reference_seq.seq.translate())
    
    # Score all ORFs and maintain original order for stable sorting
    score_data = []
    for orf in orfs:
        try:
            orf_protein = str(orf.seq.translate())
            score = AlignmentProcessor.score_sequence_similarity(
                ref_protein,
                orf_protein,
                ref_length
            )
            score_data.append((score, len(score_data), orf))
        except Exception as e:
            logging.warning(f"Error processing ORF {orf.id}: {str(e)}")
            continue
    
    if score_data:
        # Sort by score while maintaining stable order
        score_data.sort(key=lambda x: x[0], reverse=True)
        best_score, _, best_orf = score_data[0]
        
        # Log score difference between best and second best
        if len(score_data) > 1:
            score_diff = score_data[0][0] - score_data[1][0]
            logging.debug(f"Score difference between best and second best ORF for {orfs[0].id}: {score_diff}")
    
    # Use first ORF as fallback if scoring fails
    if best_orf is None and orfs:
        logging.warning(f"Using first ORF for {orfs[0].id} as fallback")
        best_orf = orfs[0]
    
    return (orfs[0].id.split('_ORF')[0], best_orf) if best_orf else None

class AIOalign:
    """
    Main class that coordinates the entire alignment pipeline.
    Handles file management and orchestrates the alignment process.
    """
    
    def __init__(self, input_dna_file: str, output_dir: str):
        """
        Initializes AIOalign with input file and output directory.
        
        Args:
            input_dna_file (str): Path to input DNA sequences
            output_dir (str): Directory for output files
        """
        self.input_dna_file = input_dna_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Define output file paths
        self.orf_dna_file = self.output_dir / "extracted_orfs.fasta"
        self.selected_orf_file = self.output_dir / "selected_orfs.fasta"
        self.protein_file = self.output_dir / "proteins.fasta"
        self.aligned_file = self.output_dir / "aligned_proteins.fasta"
        self.final_dna_file = self.output_dir / "final_aligned_dna.fasta"

    def select_best_orfs(self, reference_seq: SeqRecord, orf_records: List[SeqRecord]) -> List[SeqRecord]:
        """
        Selects the best matching ORF for each input sequence.
        
        Args:
            reference_seq (SeqRecord): Reference sequence for comparison
            orf_records (List[SeqRecord]): All possible ORFs
            
        Returns:
            List[SeqRecord]: Best matching ORFs for each sequence
        """
        # Prepare reference sequence
        clean_ref_seq = str(reference_seq.seq).replace('-', '')
        ref_protein_length = len(clean_ref_seq) // 3
        
        # Group ORFs by original sequence ID
        orf_groups = {}
        for record in orf_records:
            if '_ORF' in record.id:
                seq_id = record.id.split('_ORF')[0]
                if seq_id not in orf_groups:
                    orf_groups[seq_id] = []
                orf_groups[seq_id].append(record)
        
        # Prepare data for parallel processing
        chunks = [
            (reference_seq, orfs, ref_protein_length) 
            for orfs in orf_groups.values()
        ]
        
        # Process ORFs in parallel
        with mp.Pool(processes=mp.cpu_count() - 1) as pool:
            results = pool.map(process_sequence_chunk, chunks)
        
        # Collect results and include reference sequence
        selected_orfs = [reference_seq]
        for result in results:
            if result and result[1]:
                selected_orfs.append(result[1])
            else:
                logging.warning(f"No ORF selected for sequence {result[0] if result else 'unknown'}")
        
        # Verify all sequences are accounted for
        expected_count = len(orf_groups) + 1  # +1 for reference sequence
        if len(selected_orfs) != expected_count:
            logging.warning(f"Expected {expected_count} sequences but got {len(selected_orfs)}")
        
        return selected_orfs

    def run(self):
        """
        Executes the complete AIOalign pipeline:
        1. Extracts ORFs
        2. Selects best ORFs
        3. Translates to protein
        4. Performs alignment
        5. Back-translates to DNA
        """
        try:
            logging.info("Starting ORF extraction")
            _, reference_seq = SequenceProcessor.extract_orfs(
                self.input_dna_file, self.orf_dna_file)
            
            logging.info("Reading ORF records")
            orf_records = list(SeqIO.parse(self.orf_dna_file, "fasta"))
            
            logging.info("Selecting best ORFs")
            selected_orfs = self.select_best_orfs(reference_seq, orf_records)
            logging.info(f"Selected {len(selected_orfs)} ORFs")
            
            # Save selected ORFs
            SeqIO.write(selected_orfs, self.selected_orf_file, "fasta")
            
            # Store original DNA for back-translation
            dna_dict = {record.id: str(record.seq) for record in selected_orfs}
            
            # Translate to protein
            logging.info("Translating sequences")
            protein_sequences = SequenceProcessor.translate_sequences(selected_orfs)
            SeqIO.write(protein_sequences, self.protein_file, "fasta")
            
            # Perform final alignment
            logging.info("Running Clustal Omega alignment")
            AlignmentProcessor.run_clustalo(str(self.protein_file), str(self.aligned_file))
            
            # Back-translate aligned proteins to DNA
            logging.info("Back-translating to DNA")
            SequenceProcessor.back_translate(self.aligned_file, dna_dict, self.final_dna_file)
            
            logging.info("Program completed successfully")

        except Exception as e:
            logging.error(f"An error occurred during execution: {e}", exc_info=True)
            raise

def setup_logging(log_file: str):
    """
    Configures logging system with both file and console output.
    
    Args:
        log_file (str): Path to log file
        
    Configures:
        - INFO level logging
        - Timestamps and log levels in output
        - Simultaneous output to file and console
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

def parse_arguments():
    """
    Parses command line arguments using argparse.
    
    Returns:
        argparse.Namespace: Parsed command line arguments
        
    Arguments:
        input_file: Input FASTA file path
        -o/--output_dir: Output directory (default: 'output')
        -l/--log_file: Log file path (default: 'aioalign.log')
        --min-protein-length: Minimum protein length for ORF detection (default: 50)
    """
    parser = argparse.ArgumentParser(
        description='AIOalign v5.4: Align and select best ORFs from DNA sequences.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        'input_file',
        help='Input DNA sequences in FASTA format. First sequence must be the reference.'
    )
    parser.add_argument(
        '-o', '--output_dir',
        default='output',
        help='Output directory for results'
    )
    parser.add_argument(
        '-l', '--log_file',
        default='aioalign.log',
        help='Log file name'
    )
    parser.add_argument(
        '--min-protein-length',
        type=int,
        default=50,
        help='Minimum protein length for ORF detection'
    )
    return parser.parse_args()

if __name__ == "__main__":
    """
    Main execution block:
    1. Parse command line arguments
    2. Set up logging
    3. Initialize and run AIOalign
    4. Handle any errors that occur during execution
    """
    args = parse_arguments()
    setup_logging(args.log_file)
    
    try:
        aligner = AIOalign(args.input_file, args.output_dir)
        aligner.run()
    except Exception as e:
        logging.error(f"Program failed: {str(e)}")
        sys.exit(1)