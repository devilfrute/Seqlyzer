import streamlit as st
from Bio import SeqIO
from io import StringIO

# Set page configuration
st.set_page_config(page_title="Seqlyzer", page_icon="🔬")

# Title and introduction
st.title("🔬 Seqlyzer")
st.markdown("""
Welcome to Seqlyzer! This tool allows you to analyze DNA sequences from FASTA files. 
Upload a FASTA file to get started and see detailed information about your DNA sequence.
""")

# Upload file section
st.sidebar.header("Upload a FASTA file")
uploaded_file = st.sidebar.file_uploader("Choose a FASTA file", type="fasta")

if uploaded_file is not None:
    # Decode the uploaded file to a string
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    seq_record = SeqIO.read(stringio, "fasta")
    
    # Display sequence information
    st.subheader("Sequence Information")
    st.markdown(f"""
    - **Sequence ID:** {seq_record.id}
    - **Sequence Length:** {len(seq_record.seq)}
    """)

    # Calculate and display nucleotide counts
    a_count = seq_record.seq.count("A")
    t_count = seq_record.seq.count("T")
    g_count = seq_record.seq.count("G")
    c_count = seq_record.seq.count("C")

    st.subheader("Nucleotide Counts")
    st.markdown(f"""
    - **Adenine (A) count:** {a_count}
    - **Thymine (T) count:** {t_count}
    - **Guanine (G) count:** {g_count}
    - **Cytosine (C) count:** {c_count}
    """)

    # Calculate and display total nucleotide count
    total_nucleotides = a_count + t_count + g_count + c_count
    st.markdown(f"**Total Nucleotide Count:** {total_nucleotides}")

    # Calculate and display GC content
    gc_content = (g_count + c_count) / len(seq_record.seq) * 100
    st.markdown(f"**GC Content:** {gc_content:.2f}%")

    # Identify ORFs (simple example considering start codons only)
    orfs = [str(seq_record.seq[i:i+3]) for i in range(0, len(seq_record.seq)-2, 3) if seq_record.seq[i:i+3] == "ATG"]
    st.markdown(f"**Number of ORFs:** {len(orfs)}")

    # Create a text output for download
    output = (
        f"Sequence ID: {seq_record.id}\n"
        f"Sequence Length: {len(seq_record.seq)}\n"
        f"Nucleotide Counts:\n"
        f"  Adenine (A): {a_count}\n"
        f"  Thymine (T): {t_count}\n"
        f"  Guanine (G): {g_count}\n"
        f"  Cytosine (C): {c_count}\n"
        f"Total Nucleotide Count: {total_nucleotides}\n"
        f"GC Content: {gc_content:.2f}%\n"
        f"Number of ORFs: {len(orfs)}\n"
    )
    
    # Provide a download button for the complete analysis
    st.sidebar.download_button(
        label="Download Complete Analysis as TXT",
        data=output,
        file_name="dna_sequence_analysis.txt",
        mime="text/plain"
    )

    st.sidebar.success("Analysis complete! Check the main page for results and download options.")

else:
    st.warning("Please upload a FASTA file to analyze.")

st.sidebar.markdown("""
---
Developed by VAMSI
""")
