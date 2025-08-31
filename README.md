# AIOalign

DNA sequence alignment tool optimized for phylogenetic analysis. AIOalign prevents frame-shift misalignments by performing intelligent ORF selection and protein-level alignment while preserving codon information.

## ⚠️ Important Notice
This code represents a partially reconstructed version of a previously working tool. While the core algorithms and approach are sound, this specific implementation may contain bugs and may not function as intended. The original version successfully handled DNA sequence alignment with ORF optimization for datasets of 300+ sequences, but this reconstruction has not been thoroughly tested.
This repository is provided for:
- Reference and educational purposes
- Documentation of the approach
- Potential future development by others
- Historical record of the methodology

While pull requests to fix issues are welcome, please note that this repository is not actively maintained.

## Features
- Smart ORF detection and selection using reference-based scoring
- Strong N-terminal region weighting for biological relevance
- Preserves codon degeneracy and silent mutation information
- Parallel processing for efficient handling of large datasets (>100 sequences)
- Comprehensive logging system for debugging and quality control
- Gap handling optimized for biological sequences

## Installation

### Prerequisites
- Python 3.6 or higher
- Biopython
- Clustal Omega

### Windows
1. Install Python from [python.org](https://python.org)
2. Install required Python packages:
```bash
pip install biopython
```
3. Install Clustal Omega:
   - Download the Windows binary from [Clustal Omega](http://www.clustal.org/omega/)
   - Add Clustal Omega to your PATH or place it in the same directory as AIOalign

### Linux
```bash
# Install Python and pip if not already installed
sudo apt-get update
sudo apt-get install python3 python3-pip

# Install Clustal Omega
sudo apt-get install clustalo

# Install Biopython
pip3 install biopython
```

### macOS
```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install Python3 and Clustal Omega
brew install python3 clustal-omega

# Install Biopython
pip3 install biopython
```

## Usage

### Basic Usage
```bash
python AIOalign.py input.fasta -o output_dir
```

### Full Options
```bash
python AIOalign.py input.fasta \
    -o output_directory \     # Output directory (default: output)
    -l logfile.log \         # Log file (default: aioalign.log)
    --min-protein-length 50  # Minimum protein length for ORF detection (default: 50)
```

### Input Requirements
- FASTA format DNA sequences
- First sequence must be the reference sequence

### Output Files
- `extracted_orfs.fasta`: All potential ORFs found
- `selected_orfs.fasta`: Best ORF selected for each sequence
- `proteins.fasta`: Translated protein sequences
- `aligned_proteins.fasta`: Aligned protein sequences
- `final_aligned_dna.fasta`: Final DNA alignment
- `guidetree.txt`: Phylogenetic guide tree from Clustal Omega
- `clusters.txt`: Clustering information from alignment

## How It Works
1. Finds all possible open reading frames (ORFs) in each sequence
2. Translates ORFs to amino acid sequences
3. Selects best ORF based on:
   - Similarity to reference sequence
   - N-terminal region quality
   - Sequence length compatibility
   - Gap patterns
4. Performs protein-level alignment using Clustal Omega
5. Back-translates to DNA while maintaining alignment

## Acknowledgments

- Built using BioPython and Clustal Omega
- Developed with assistance from Anthropic's Claude 3.5 Sonnet AI

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
