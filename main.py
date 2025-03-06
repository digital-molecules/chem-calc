import streamlit as st
from chem_calc import compute_properties, compute_lip, compute_similarity

st.title("Thally's QSAR Toolkit")
st.write("Thanks to the brilliant Suliman Sharif for his public resources.")

smiles1 = st.text_input("**Enter a valid SMILES notation please:**")

if smiles1:
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Molecular Properties")
        props = compute_properties(smiles1)
        st.write("Computed Properties:", props)

        lip_info = compute_lip(smiles1)
        if lip_info:
            st.success("This molecule follows Lipinski's Rule of 5! 😺")
        else:
            st.warning("This molecule does **NOT** follow Lipinski's Rule of 5. 😿")

smiles2 = st.text_input("**Enter another SMILES for similarity comparison:**")
if smiles1 and smiles2:
    sim = compute_similarity(smiles1, smiles2)
    st.write("Tanimoto Similarity:", sim)



    #img = render_molecule_image(smiles1)
    #if img:
        #st.image(img, caption="Molecule Structure")


    
#if smiles1:
    #funct = detect_functional_groups(smiles1)
    #st.write("Functional Group(s):", funct)
