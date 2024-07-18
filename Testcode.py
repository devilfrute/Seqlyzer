import streamlit as st
from Bio import SeqIO
from io import BytesIO, StringIO

st.title("Seqlyzer: DNA Sequence Analyzer")

# Upload a FASTA file
uploaded_file = st.file_uploader("Upload a FASTA file", type="fasta")

if uploaded_file is not None:
    # Decode the uploaded file to a string
    file_contents = uploaded_file.getvalue().decode("utf-8")
    
    # Convert to StringIO object
    stringio = StringIO(file_contents)
    
    try:
        # Read the sequence record from the StringIO object
        seq_record = SeqIO.read(stringio, "fasta")
        
        # Display sequence information
        st.header("Sequence Information")
        st.write(f"Sequence ID: {seq_record.id}")
        st.write(f"Sequence Length: {len(seq_record.seq)}")
        
        # Calculate GC content
        g_count = seq_record.seq.count("G")
        c_count = seq_record.seq.count("C")
        gc_content = (g_count + c_count) / len(seq_record.seq) * 100
        st.write(f"GC Content: {gc_content:.2f}%")
        
        # Detect ORFs (example: finding start codons)
        orfs = []
        for i in range(len(seq_record.seq) - 2):
            if seq_record.seq[i:i+3] == "ATG":
                orf = seq_record.seq[i:]
                orfs.append(orf)
        
        # Display ORFs
        st.header("Open Reading Frames (ORFs)")
        if orfs:
            st.write(f"Number of ORFs found: {len(orfs)}")
            for i, orf in enumerate(orfs):
                st.write(f"ORF {i+1}: {orf}")
        else:
            st.write("No ORFs found.")
        
        # Translate DNA sequence to protein sequence
        protein_seq = seq_record.seq.translate(to_stop=True)
        
        # Display translated protein sequence
        st.header("Translated Protein Sequence")
        st.text_area("Protein Sequence", str(protein_seq), height=200)
    
    except Exception as e:
        st.error(f"Error reading file: {str(e)}")
