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
        
        # Visualize the DNA sequence
        st.header("Sequence Visualization")
        seq_str = str(seq_record.seq)
        st.text_area("DNA Sequence", seq_str, height=300)
        
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
        
        # Display total number of ORFs detected
        st.write(f"Total ORFs found: {len(orfs)}")
        
        # Get user input for ORF numbers to display in detail
        orf_numbers_input = st.text_input("Enter ORF numbers (comma-separated)", "")
        orf_numbers = [int(num.strip()) for num in orf_numbers_input.split(',') if num.strip().isdigit()]
        
        # Validate ORF numbers entered by the user
        invalid_orfs = [num for num in orf_numbers if num <= 0 or num > len(orfs)]
        
        if invalid_orfs:
            st.error(f"Invalid ORF numbers: {', '.join(map(str, invalid_orfs))}. Enter valid ORF numbers.")
        
        else:
            for orf_num in orf_numbers:
                orf_index = orf_num - 1  # Convert to zero-indexed for list access
                if orf_index < len(orfs):
                    st.subheader(f"ORF {orf_num}")
                    st.write(f"Sequence: {orfs[orf_index]}")
                    st.download_button(
                        label=f"Download ORF {orf_num} Sequence",
                        data=str(orfs[orf_index]),
                        file_name=f"orf_{orf_num}.txt",
                        mime="text/plain"
                    )
    
    except Exception as e:
        st.error(f"Error reading file: {str(e)}")
