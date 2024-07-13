import streamlit as st
from Bio import SeqIO
from io import StringIO

st.title("Seqlyzer")

uploaded_file = st.file_uploader("Upload a FASTA file", type="fasta")

if uploaded_file is not None:
    # Decode the uploaded file to a string
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    seq_record = SeqIO.read(stringio, "fasta")
    
    st.write(f"**Sequnece ID  :** {seq_record.id}")
    st.write(f"**Sequence Length :** {len(seq_record.seq)}")
    
    # Calculate and display nucleotide counts
    a_count = seq_record.seq.count("A")
    t_count = seq_record.seq.count("T")
    g_count = seq_record.seq.count("G")
    c_count = seq_record.seq.count("C")
    
    st.write(f"**Nucleotide Counts as follows :**")
    st.write(f"Adenine (A) count: {a_count}")
    st.write(f"Thymine (T) count: {t_count}")
    st.write(f"Guanine (G) count: {g_count}")
    st.write(f"Cytosine (C) count: {c_count}")
    
    # Calculate and display total nucleotide count
    total_nucleotides = a_count + t_count + g_count + c_count
    st.write(f"**Total Nucleotide Count:** {total_nucleotides}")
    
    # Calculate and display GC content
    gc_content = (g_count + c_count) / len(seq_record.seq) * 100
    st.write(f"**GC Content:** {gc_content:.2f}%")
    
    # Identify ORFs (simple example considering start codons only)
    orfs = [str(seq_record.seq[i:i+3]) for i in range(0, len(seq_record.seq)-2, 3) if seq_record.seq[i:i+3] == "ATG"]
    st.write(f"**Number of ORFs:** {len(orfs)}")
    
    # Translate DNA sequence to protein sequence
    protein_seq = seq_record.seq.translate(to_stop=True)
    
    # Display the translated protein sequence in a text area
    st.write("**Translated Protein Sequence:**")
    st.text_area("Protein Sequence", str(protein_seq), height=200)
