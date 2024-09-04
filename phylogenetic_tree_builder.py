import streamlit as st
from Bio import Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio.Align import MultipleSeqAlignment
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import matplotlib.pyplot as plt
import os

# Apply custom CSS for theming
st.markdown("""
    <style>
        body {
            background-color: #f0f0f0; /* Light grey background */
        }
        .title {
            color: #1f77b4; /* Blue color for the title */
        }
        .header {
            background-color: #ffffff; /* White background for the header */
            padding: 20px;
            border-bottom: 3px solid #1f77b4; /* Blue bottom border for the header */
        }
        .subheader {
            color: #ff7f0e; /* Orange color for subheaders */
        }
        .text-area {
            border: 2px solid #1f77b4; /* Blue border for text areas */
            border-radius: 5px;
        }
        .footer {
            text-align: center;
            color: #888888; /* Grey color for the footer */
        }
        .table {
            background-color: #ffffff; /* White background for tables */
            border-radius: 5px;
        }
    </style>
""", unsafe_allow_html=True)

# Add a full-width logo at the top
logo = "logo.jpg"
st.image(logo, use_column_width=True)  # Logo spans the full width of the column


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
st.write("**Created by Ahmed Saif**")

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

        # Generate tree visualization as an image file
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, do_show=False, axes=ax)
        
        # Save tree image to a file
        image_path = "phylogenetic_tree.png"
        fig.savefig(image_path)

        # Display the image in Streamlit
        st.image(image_path, caption="Phylogenetic Tree", use_column_width=True)

        # Download option for the tree in Newick format
        newick_str = StringIO()
        Phylo.write(tree, newick_str, "newick")
        st.download_button(label="Download Tree in Newick Format", 
                           data=newick_str.getvalue(), 
                           file_name="phylogenetic_tree.newick",
                           mime="text/plain")
        
        # Clean up the saved image file
        os.remove(image_path)

    else:
        st.error("Please upload a FASTA file with more than one sequence.")
else:
    st.info("Waiting for a FASTA file upload...")

# Footer
st.write("---")
st.write("Developed with Streamlit")
