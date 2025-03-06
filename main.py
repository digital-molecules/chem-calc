import streamlit as st
from chem_calc import compute_properties, render_molecule, compute_lipinski, compute_similarity

st.title("Thally's QSAR Toolkit")
st.write("*Thanks to the brilliant Suliman Sharif for his public resources.*")

smiles1 = st.text_input("**Enter a valid SMILES notation please:**")
if smiles1:
    props = compute_properties(smiles1)
    st.write("Computed Properties:", props)

mol_img = render_molecule(smiles1)
if mol_img:
    st.image(mol_img, caption="Molecular Structure", use_container_width=True)
    
if smiles1:
    lip_info = compute_lipinski(smiles1)
    if qed_info:
        st.success("This molecule passes Lipinski's Rule of 5! 😺")
    else:
        st.warning("This molecule does **NOT** pass Lipinski's Rule of 5. 😿")

smiles2 = st.text_input("**Enter another SMILES for similarity comparison:**")
if smiles1 and smiles2:
    sim = compute_similarity(smiles1, smiles2)
    st.write("Tanimoto Similarity:", sim)






    
#if smiles1:
    #funct = detect_functional_groups(smiles1)
    #st.write("Functional Group(s):", funct)
