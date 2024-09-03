import streamlit as st
from Bio import Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio.Align import MultipleSeqAlignment
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import matplotlib.pyplot as plt

# Function to perform pairwise sequence alignment and convert to SeqRecord
def align_sequences(sequences):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    
    # Reference sequence for alignment
    ref_seq = sequences[0]
    aligned_sequences = []
    
    # Align each sequence to the reference sequence
    for seq in sequences:
        alignment = aligner.align(ref_seq.seq, seq.seq)[0]
        aligned_seq = SeqRecord(Seq(str(alignment.target)), id=seq.id)
        aligned_sequences.append(aligned_seq)
    
    return aligned_sequences

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
        # Align sequences using pairwise aligner
        aligned_sequences = align_sequences(sequences)

        # Create a MultipleSeqAlignment object
        alignment = MultipleSeqAlignment(aligned_sequences)

        # Calculate the distance matrix based on sequence identity
        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)

        # Build a phylogenetic tree using the Neighbor-Joining method
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(distance_matrix)

        # Display the tree using matplotlib
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, ax=ax)
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
