import streamlit as st
from chem_utils import compute_properties, detect_functional_groups, compute_similarity

st.title("Chemical Properties Calculator")

smiles1 = st.text_input("Enter a valid SMILES notation:")
if smiles1:
    props = compute_properties(smiles1)
    st.write("Computed Properties:", props)

if smiles1:
    funct = detect_functional_groups(smiles1)
    st.write("Functional Group(s):", funct)

smiles2 = st.text_input("Enter another SMILES for similarity comparison:")
if smiles1 and smiles2:
    sim = compute_similarity(smiles1, smiles2)
    st.write("Tanimoto Similarity:", sim)
