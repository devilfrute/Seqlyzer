import streamlit as st
from Bio import Entrez, SeqIO
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.feature_extraction import DictVectorizer
from sklearn.metrics import accuracy_score, confusion_matrix

# Function to fetch sequences from NCBI
def fetch_sequences(email, gene="BRCA1", organism="Homo sapiens", retmax=10):
    Entrez.email = email
    handle = Entrez.esearch(db="nucleotide", term=f"{gene}[Gene] AND {organism}[Organism]", retmax=retmax)
    record = Entrez.read(handle)
    ids = record["IdList"]
    
    sequences = []
    for seq_id in ids:
        handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
        seq_record = SeqIO.read(handle, "fasta")
        sequences.append(seq_record.seq)
    return sequences

# Function to extract k-mer features from sequences
def kmer_count(seq, k=3):
    kmer_dict = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer in kmer_dict:
            kmer_dict[kmer] += 1
        else:
            kmer_dict[kmer] = 1
    return kmer_dict

# Streamlit app interface
st.title("SNP Detection and Classification")

email = st.text_input("Enter your email:")
gene = st.text_input("Enter the gene (default: BRCA1):", "BRCA1")
organism = st.text_input("Enter the organism (default: Homo sapiens):", "Homo sapiens")
retmax = st.number_input("Number of sequences to fetch (default: 10):", 1, 100, 10)

if st.button("Fetch and Analyze Sequences"):
    with st.spinner("Fetching sequences from NCBI..."):
        sequences = fetch_sequences(email, gene, organism, retmax)
    
    if sequences:
        st.success("Sequences fetched successfully!")

        # Extracting features
        with st.spinner("Extracting k-mer features..."):
            features = [kmer_count(str(seq)) for seq in sequences]

        # Simulated labels (for simplicity, let's create dummy labels)
        labels = [1 if i % 2 == 0 else 0 for i in range(len(features))]

        # Vectorizing features
        vectorizer = DictVectorizer(sparse=False)
        X = vectorizer.fit_transform(features)
        y = labels

        # Splitting data
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

        # Training model
        model = LogisticRegression()
        model.fit(X_train, y_train)

        # Making predictions
        predictions = model.predict(X_test)

        # Evaluating model
        accuracy = accuracy_score(y_test, predictions)
        confusion = confusion_matrix(y_test, predictions)

        st.subheader("Model Evaluation")
        st.write(f"Accuracy: {accuracy}")
        st.write("Confusion Matrix:")
        st.write(confusion)
    else:
        st.error("No sequences found. Please check your input and try again.")
