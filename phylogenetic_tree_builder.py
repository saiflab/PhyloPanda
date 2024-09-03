
import streamlit as st
from Bio import Phylo, AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from io import StringIO
import matplotlib.pyplot as plt

# Title and description
st.title("Phylogenetic Tree Builder")
st.write("Upload a DNA FASTA file and generate a phylogenetic tree.")
st.write("Created by Ahmed Saif")

# Upload FASTA file
uploaded_file = st.file_uploader("Choose a FASTA file", type="fasta")

if uploaded_file is not None:
    # Save the uploaded file content to a string
    fasta_content = uploaded_file.read().decode("utf-8")
    
    # Write the content to a temporary file
    with open("sequences.fasta", "w") as f:
        f.write(fasta_content)
    
    # Run Clustal Omega to align sequences (you need to have Clustal Omega installed)
    clustalomega_cline = ClustalOmegaCommandline(infile="sequences.fasta", outfile="aligned.aln", verbose=True, auto=True)
    clustalomega_cline()
    
    # Read the aligned sequences
    aligned = AlignIO.read("aligned.aln", "clustal")
    
    # Create a phylogenetic tree using the aligned sequences
    with open("aligned.aln", "r") as aln_file:
        alignment = AlignIO.read(aln_file, "clustal")
    
    # Generate the tree from alignment
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    
    # Display the tree
    Phylo.draw(tree)
    st.pyplot(plt)
    
    # Optionally download the tree
    st.download_button(label="Download Newick Tree", data=tree.format("newick"), file_name="phylogenetic_tree.newick", mime="text/plain")

# Footer
st.write("---")
st.write("Developed with Streamlit")
