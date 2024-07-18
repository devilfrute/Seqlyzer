import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC, molecular_weight
import re

# Function to calculate GC content and ATGC count
def calculate_gc_atgc(sequence):
    gc_content = GC(sequence)
    atgc_counts = {
        'A': sequence.count('A'),
        'T': sequence.count('T'),
        'G': sequence.count('G'),
        'C': sequence.count('C')
    }
    return gc_content, atgc_counts

# Function to find ORFs
def find_orfs(sequence):
    orfs = []
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    
    for match in re.finditer(r'(?=(ATG(?:...)*?)(?=TAA|TAG|TGA))', str(sequence)):
        orf = match.group()
        orfs.append(orf)
    
    return orfs

# Function to detect SNPs
def detect_snps(sequence, reference):
    snps = []
    for i, (nuc_seq, nuc_ref) in enumerate(zip(sequence, reference)):
        if nuc_seq != nuc_ref:
            snps.append(f"Position {i+1}: {nuc_ref} -> {nuc_seq}")
    return snps

# Function to detect restriction sites
def detect_restriction_sites(sequence):
    restriction_sites = {
        'EcoRI': 'GAATTC',
        'BamHI': 'GGATCC',
        'HindIII': 'AAGCTT'
    }
    found_sites = []
    for name, site in restriction_sites.items():
        matches = [m.start() for m in re.finditer(site, str(sequence))]
        found_sites.extend([(name, match) for match in matches])
    return found_sites

# Main Streamlit app
def main():
    st.title("Seqlyzer: DNA Sequence Analyzer")
    
    # Upload a FASTA file
    uploaded_file = st.file_uploader("Upload a FASTA file", type="fasta")
    
    if uploaded_file is not None:
        # Parse the uploaded file
        seq_record = SeqIO.read(uploaded_file, "fasta")
        sequence = seq_record.seq
        
        # Calculate GC content and ATGC count
        gc_content, atgc_counts = calculate_gc_atgc(sequence)
        
        # Display sequence information
        st.header("Sequence Information")
        st.write(f"Sequence ID: {seq_record.id}")
        st.write(f"Sequence Length: {len(sequence)}")
        
        # Display GC content and ATGC counts
        st.header("Nucleotide Analysis")
        st.write(f"GC Content: {gc_content:.2f}%")
        st.write("ATGC Counts:")
        st.write(f"A: {atgc_counts['A']}")
        st.write(f"T: {atgc_counts['T']}")
        st.write(f"G: {atgc_counts['G']}")
        st.write(f"C: {atgc_counts['C']}")
        
        # Find ORFs
        orfs = find_orfs(sequence)
        st.header("Open Reading Frames (ORFs)")
        if orfs:
            st.write(f"Number of ORFs found: {len(orfs)}")
            for i, orf in enumerate(orfs):
                st.write(f"ORF {i+1}: {orf}")
        else:
            st.write("No ORFs found.")
        
        # Translate DNA sequence to protein sequence
        protein_seq = sequence.translate(to_stop=True)
        
        # Display translated protein sequence
        st.header("Translated Protein Sequence")
        st.text_area("Protein Sequence", str(protein_seq), height=200)
        
        # Detect SNPs (compare with a reference sequence)
        reference_sequence = Seq("ATGCGTACGCGT")
        snps = detect_snps(sequence, reference_sequence)
        
        # Display SNPs
        st.header("Single Nucleotide Polymorphisms (SNPs)")
        if snps:
            for snp in snps:
                st.write(snp)
        else:
            st.write("No SNPs detected.")
        
        # Detect restriction sites
        restriction_sites = detect_restriction_sites(sequence)
        
        # Display restriction sites
        st.header("Restriction Sites")
        if restriction_sites:
            for site in restriction_sites:
                st.write(f"{site[0]} site are found at position {site[1]}")
        else:
            st.write("No restriction sites found.")
    
    else:
        st.write("Upload a FASTA file to get started.")

if __name__ == "__main__":
    main()
