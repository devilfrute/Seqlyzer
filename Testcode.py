import streamlit as st
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Function to find ORFs (Open Reading Frames)
def find_orfs(sequence):
    orfs = []
    for i in range(0, len(sequence)-2, 3):
        if sequence[i:i+3] == "ATG":
            orf_start = i
            for j in range(i+3, len(sequence)-2, 3):
                codon = sequence[j:j+3]
                if codon in ["TAA", "TAG", "TGA"]:
                    orf_end = j + 3
                    orfs.append((orf_start, orf_end))
                    break
    return orfs

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

    # Find ORFs
    orfs = find_orfs(seq_record.seq)

    # Display number of ORFs found
    st.subheader("Open Reading Frames (ORFs)")
    st.markdown(f"**Number of ORFs:** {len(orfs)}")

    # Provide option to view individual ORFs
    orf_to_view = st.number_input("Enter the ORF number to view (1 to {len(orfs)})", min_value=1, max_value=len(orfs), step=1, value=1)

    if st.button("View Selected ORF"):
        if orf_to_view > len(orfs) or orf_to_view < 1:
            st.warning(f"Invalid ORF number. Enter a number between 1 and {len(orfs)}")
        else:
            orf_start, orf_end = orfs[orf_to_view-1]
            selected_orf = seq_record.seq[orf_start:orf_end]
            st.write(f"**Selected ORF {orf_to_view}**: Start={orf_start}, End={orf_end}, Length={len(selected_orf)}")
            st.write(selected_orf)

            # Provide a download button for the selected ORF
            orf_record = SeqRecord(selected_orf, id=f"ORF_{orf_to_view}", description=f"ORF {orf_to_view}")
            st.sidebar.download_button(
                label=f"Download ORF {orf_to_view} as FASTA",
                data=orf_record.format("fasta"),
                file_name=f"orf_{orf_to_view}.fasta",
                mime="text/plain"
            )

# Sidebar footer
st.sidebar.markdown("""
---
Developed by VAMSI
""")
