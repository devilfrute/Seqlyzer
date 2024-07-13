import streamlit as st
from Bio import SeqIO
from io import StringIO

st.title("Seqlyzer")

uploaded_file = st.file_uploader("Upload a FASTA file", type="fasta")

if uploaded_file is not None:
    # Decode the uploaded file to a string
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    seq_record = SeqIO.read(stringio, "fasta")
    
    seq_id = f"**Sequence ID  :** {seq_record.id}"
    seq_length = f"**Sequence Length :** {len(seq_record.seq)}"
    
    # Calculate and display nucleotide counts
    a_count = seq_record.seq.count("A")
    t_count = seq_record.seq.count("T")
    g_count = seq_record.seq.count("G")
    c_count = seq_record.seq.count("C")
    
    nucleotide_counts = (
        f"**Nucleotide Counts as follows :**\n"
        f"Adenine (A) count: {a_count}\n"
        f"Thymine (T) count: {t_count}\n"
        f"Guanine (G) count: {g_count}\n"
        f"Cytosine (C) count: {c_count}\n"
    )
    
    # Calculate and display total nucleotide count
    total_nucleotides = a_count + t_count + g_count + c_count
    total_nucleotides_str = f"**Total Nucleotide Count:** {total_nucleotides}"
    
    # Calculate and display GC content
    gc_content = (g_count + c_count) / len(seq_record.seq) * 100
    gc_content_str = f"**GC Content:** {gc_content:.2f}%"
    
    # Identify ORFs (simple example considering start codons only)
    orfs = [str(seq_record.seq[i:i+3]) for i in range(0, len(seq_record.seq)-2, 3) if seq_record.seq[i:i+3] == "ATG"]
    num_orfs = f"**Number of ORFs:** {len(orfs)}"
    
    # Translate DNA sequence to protein sequence, stopping at the first stop codon
    protein_seq = seq_record.seq.translate(to_stop=True)
    translated_seq_str = f"**Translated Protein Sequence:**\n{protein_seq}"
    
    # Display results in Streamlit
    st.write(seq_id)
    st.write(seq_length)
    st.write(nucleotide_counts)
    st.write(total_nucleotides_str)
    st.write(gc_content_str)
    st.write(num_orfs)
    st.write(translated_seq_str)
    
    # Display the translated protein sequence in a text area
    st.text_area("Translated Protein Sequence", str(protein_seq), height=200)
    
    # Create a text output for download
    output = (
        f"{seq_id}\n"
        f"{seq_length}\n"
        f"{nucleotide_counts}\n"
        f"{total_nucleotides_str}\n"
        f"{gc_content_str}\n"
        f"{num_orfs}\n"
        f"{translated_seq_str}\n"
    )
    
    # Provide a download button for the complete analysis
    st.download_button(
        label="Download Complete Analysis as TXT",
        data=output,
        file_name="dna_sequence_analysis.txt",
        mime="text/plain"
    )
    
    # Provide a download button for the translated protein sequence only
    st.download_button(
        label="Download Translated Protein Sequence as TXT",
        data=str(protein_seq),
        file_name="translated_protein_sequence.txt",
        mime="text/plain"
    )
