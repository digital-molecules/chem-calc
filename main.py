import streamlit as st
from chem_calc import compute_properties, compute_qed, compute_lip, molecular_formula, compute_similarity


with st.sidebar:
    st.header("QSAR")
    st.subheader("🤔 What is QSAR and why is it relevant to cheminformatics?")
    st.write("💡 This is text")
    st.subheader("🤔 What are the parameters of QSAR?")
    st.write("💡 This is text")
    st.subheader("🤔 Why is QSAR useful?")
    st.write("💡 This is text")
    st.subheader("🤔 What is the SAR paradox?")
    st.write("💡 This is text")

    st.header("LogP")
    st.subheader("🤔 What is LogP and why is it important?")
    st.write("💡 This is text")
    st.subheader("🤔 How do we interpet LogP values?")
    st.write("💡 This is text")

    st.header("QED")
    st.subheader("🤔 What is QED and why is it important?")
    st.write("💡 This is text")
    st.subheader("🤔 How do we interpet QED values?")
    st.write("💡 This is text")
    
    st.header("Tanimoto Index")
    st.subheader("🤔 What is Tanimoto Index and why is it useful?")
    st.write("💡 This is text")
    st.subheader("🤔 How do we interpet the results of the comparison?")
    st.write("💡 This is text")
    
    st.header("Lipinski's Rule of 5")
    st.subheader("🤔 What is 'Lipinski's Rule of 5'?")
    st.write("💡 This is text")
    st.subheader("🤔 How accurate is it?")
    st.write("💡 This is text")

st.title("Metrichemical: A chemist's pocket toolbox")
st.write("Thanks to the brilliant Suliman Sharif and Vaneet Saini for their public resources.")

smiles1 = st.text_input("**Please enter a valid SMILES:**")

if smiles1:
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Molecular Properties")
        props = compute_properties(smiles1)
        st.write(props)

    with col2:
        st.subheader("Molecular Formula")
        formula = molecular_formula(smiles1)
        st.write(f"Molecular Formula: {formula}")
        st.subheader("Tanimoto Similarity")
        smiles2 = st.text_input("**Enter another SMILES for similarity comparison:**")
        if smiles1 and smiles2:
            sim = compute_similarity(smiles1, smiles2)
            st.write(sim)
            
if smiles1:
    qed = compute_qed(smiles1)
    st.markdown(
        f"<p style='font-size: 20px; font-weight: bold; color: #a20f66;'>QED Drug-Likeness: {qed}</p>",
        unsafe_allow_html=True
    )
            

if smiles1:
    lip_info = compute_lip(smiles1)
    if lip_info:
        st.success("This molecule follows Lipinski's Rule of 5! 😺")
    else:
        st.warning("This molecule does **NOT** follow Lipinski's Rule of 5. 😿")

st.markdown("---")
st.markdown("<h1 style='text-align: center; font-size: 20px; font-weight: 100; font-style: italic;'>Make sure to refer to the sidebar for more information</h1>", unsafe_allow_html=True)


    #img = render_molecule_image(smiles1)
    #if img:
        #st.image(img, caption="Molecule Structure")


    
#if smiles1:
    #funct = detect_functional_groups(smiles1)
    #st.write("Functional Group(s):", funct)
