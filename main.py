import streamlit as st
from chem_calc import compute_properties, compute_qed, compute_similarity

st.title("Thally's QSAR Toolkit")
st.write("*Thanks to the brilliant Suliman Sharif for his public resources and my younger brother Chris who was (mostly) available to help with the debugging.*")

smiles1 = st.text_input("**Enter a valid SMILES notation please:**")
if smiles1:
    props = compute_properties(smiles1)
    st.write("Computed Properties:", props)
    
if smiles1:
    qed_info = qed_drug_likeness(smiles1)
    if "Error" in qed_info:
            st.error(qed_info["Error"])
        else:
            if qed_info["Passes QED Drug-Likeness"]:
                    st.success("✅ This molecule meets the QED drug-likeness criteria! 😺")
                else:
                    st.warning("❌ This molecule does **NOT** meet the QED drug-likeness criteria. 😿")

smiles2 = st.text_input("**Enter another SMILES for similarity comparison:**")
if smiles1 and smiles2:
    sim = compute_similarity(smiles1, smiles2)
    st.write("Tanimoto Similarity:", sim)






    
#if smiles1:
    #funct = detect_functional_groups(smiles1)
    #st.write("Functional Group(s):", funct)
