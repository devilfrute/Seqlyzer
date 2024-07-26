import streamlit as st
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq

# Set page configuration
st.set_page_config(page_title="Seqlyzer", page_icon="🧬")

# Title and introduction
st.title("🧬 Seqlyzer")
st.markdown("""
Welcome to Seqlyzer! This tool allows you to analyze DNA sequences from FASTA files.
Upload a FASTA file to get started and see detailed information about your DNA sequence.
""")

# Upload file section
st.sidebar.header("🧬 Upload a FASTA file")
uploaded_file = st.sidebar.file_uploader("Choose a FASTA file", type="fasta")

# Initialize the session state for toggling visibility
if "show_details" not in st.session_state:
    st.session_state.show_details = False

if "show_translation" not in st.session_state:
    st.session_state.show_translation = False

if uploaded_file is not None:
    # Decode the uploaded file to a string
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    seq_record = SeqIO.read(stringio, "fasta")

    # Display sequence information
    st.subheader("🧬 Sequence Information")
    st.markdown(f"""
    - **Sequence ID:** {seq_record.id}
    - **Sequence Length:** {len(seq_record.seq)}
    """)

    # Toggle visibility of details
    if st.button("🧬 Toggle Details"):
        st.session_state.show_details = not st.session_state.show_details

    # Show or hide details based on the session state
    if st.session_state.show_details:
        # Calculate and display nucleotide counts
        a_count = seq_record.seq.count("A")
        t_count = seq_record.seq.count("T")
        g_count = seq_record.seq.count("G")
        c_count = seq_record.seq.count("C")

        st.subheader("🧬 Nucleotide Counts")
        st.markdown(f"""
        - **Adenine (A) count:** {a_count}
        - **Thymine (T) count:** {t_count}
        - **Guanine (G) count:** {g_count}
        - **Cytosine (C) count:** {c_count}
        """)

        # Calculate and display total nucleotide count
        total_nucleotides = a_count + t_count + g_count + c_count
        st.markdown(f"🧬 Total Nucleotide Count: {total_nucleotides}")

        # Calculate and display GC content
        gc_content = (g_count + c_count) / len(seq_record.seq) * 100
        st.markdown(f"**🧬 GC Content:** {gc_content:.2f}%")

    # Toggle visibility of translation
    if st.button("🧬 Show Translation"):
        st.session_state.show_translation = not st.session_state.show_translation

    # Show or hide translation based on the session state
    if st.session_state.show_translation:
        # Translate the sequence
        protein_seq = seq_record.seq.translate()
        st.subheader("🧬 Translated Protein Sequence")
        st.text(protein_seq)

# Sidebar footer
st.sidebar.markdown("""
---
Developed by VAMSI S
""")
