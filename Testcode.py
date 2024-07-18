import streamlit as st
from Bio import SeqIO
from io import StringIO

st.set_page_config(layout="wide")  # Set wide layout for better aesthetics

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
        
        # Identify ORFs (simple example considering start codons only)
        orfs = []
        for i in range(len(seq_record.seq) - 2):
            if seq_record.seq[i:i+3] == "ATG":
                orf = seq_record.seq[i:]
                orfs.append(orf)
        
        # Display all ORFs
        st.header("Open Reading Frames (ORFs)")
        
        for idx, orf in enumerate(orfs, start=1):
            st.subheader(f"ORF {idx}")
            st.write(f"Sequence: {orf}")
            st.download_button(
                label=f"Download ORF {idx} Sequence",
                data=str(orf),
                file_name=f"orf_{idx}.txt",
                mime="text/plain"
            )
    
    except Exception as e:
        st.error(f"Error reading file: {str(e)}")
