import streamlit as st
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from io import StringIO
import matplotlib.pyplot as plt
import numpy as np

# Function to pad sequences to the same length
def pad_sequences(sequences):
    max_len = max(len(record.seq) for record in sequences)
    for record in sequences:
        if len(record.seq) < max_len:
            record.seq = record.seq + "-" * (max_len - len(record.seq))
    return sequences

# Title and description
st.title("Simple Phylogenetic Tree Builder")
st.write("Upload a DNA FASTA file to generate a phylogenetic tree.")
st.write("Created by Ahmed Saif")

# Upload FASTA file
uploaded_file = st.file_uploader("Choose a FASTA file", type="fasta")

if uploaded_file is not None:
    # Decode the uploaded file from binary to text mode
    fasta_text = StringIO(uploaded_file.getvalue().decode("utf-8"))
    
    # Parse the FASTA file
    sequences = list(SeqIO.parse(fasta_text, "fasta"))

    if len(sequences) > 1:
        # Pad sequences to the same length (simulating alignment)
        sequences = pad_sequences(sequences)

        # Create a MultipleSeqAlignment object
        alignment = MultipleSeqAlignment(sequences)

        # Calculate the distance matrix based on sequence identity
        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)

        # Build a phylogenetic tree using the Neighbor-Joining method
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(distance_matrix)

        # Display the tree using matplotlib
        fig = plt.figure(figsize=(10, 6))
        Phylo.draw(tree, do_show=False)
        st.pyplot(fig)

        # Download option for the tree in Newick format
        newick_str = StringIO()
        Phylo.write(tree, newick_str, "newick")
        st.download_button(label="Download Tree in Newick Format", 
                           data=newick_str.getvalue(), 
                           file_name="phylogenetic_tree.newick",
                           mime="text/plain")
    else:
        st.error("Please upload a FASTA file with more than one sequence.")
else:
    st.info("Waiting for a FASTA file upload...")

# Footer
st.write("---")
st.write("Developed with Streamlit")
