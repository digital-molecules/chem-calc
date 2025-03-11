import streamlit as st
from chem_calc import compute_mw, compute_logp, compute_tpsa, compute_rotbond, compute_hbd, compute_hba, compute_qed, compute_rules, compute_similarity


with st.sidebar:
    st.title("QSAR")
    st.subheader("🤔 What is QSAR and why is it relevant to chemoinformatics?")
    with open("info/1qsar.txt", "r") as file:
        qsar = file.read()
    st.markdown(qsar, unsafe_allow_html=True)
    st.subheader("🤔 What are the parameters of QSAR?")
    with open("info/2descriptors.txt", "r") as file:
        descriptors = file.read()
    st.markdown(descriptors, unsafe_allow_html=True)
    st.subheader("🤔 What is the SAR paradox?")
    with open("info/3sar_paradox.txt", "r") as file:
        sar = file.read()
    st.markdown(sar, unsafe_allow_html=True)
    st.markdown("---")
    st.header("LogP")
    st.subheader("🤔 What is logP and why is it important?")
    with open("info/4logp.txt", "r") as file:
        logp = file.read()
    st.markdown(logp, unsafe_allow_html=True)
    st.subheader("🤔 How do we interpret LogP values?")
    with open("info/5logp_values.txt", "r") as file:
        logp_values = file.read()
    st.markdown(logp_values, unsafe_allow_html=True)
    st.markdown("---")
    st.header("TPSA")
    st.subheader("🤔 What is TPSA and why is it important?")
    with open("info/6tpsa.txt", "r") as file:
        tpsa = file.read()
    st.markdown(tpsa, unsafe_allow_html=True)
    st.subheader("🤔 How do we interpret TPSA values?")
    with open("info/7tpsa_values.txt", "r") as file:
        tpsa_values = file.read()
    st.markdown(tpsa_values, unsafe_allow_html=True)
    st.markdown("---")    
    st.header("QED")
    st.subheader("🤔 What is QED and why is it important?")
    with open("info/8qed.txt", "r") as file:
        qed = file.read()
    st.markdown(qed, unsafe_allow_html=True)
    st.subheader("🤔 How do we interpret QED values?")
    with open("info/9qed_values.txt", "r") as file:
        qed_values = file.read()
    st.markdown(qed_values, unsafe_allow_html=True)
    st.markdown("---")    
    st.header("Tanimoto Index")
    st.subheader("🤔 What is Tanimoto Index and why is it useful?")
    with open("info/10tanimoto.txt", "r") as file:
        tanimoto = file.read()
    st.markdown(tanimoto, unsafe_allow_html=True)
    st.subheader("🤔 What are Morgan Fingerprints?")
    with open("info/11morgan.txt", "r") as file:
        morgan = file.read()
    st.markdown(morgan, unsafe_allow_html=True)
    st.subheader("🤔 How do we interpret the results of the comparison?")
    with open("info/12comparison_values.txt", "r") as file:
        comparison_values = file.read()
    st.markdown(comparison_values, unsafe_allow_html=True)
    st.markdown("---")    
    st.header("Lipinski's Rule of 5")
    st.subheader("🤔 What is 'Lipinski's Rule of 5'?")
    with open("info/13lipinski.txt", "r") as file:
        lipinski = file.read()
    st.markdown(lipinski, unsafe_allow_html=True)
    st.markdown("---")
    st.header("Veber's rule")
    st.subheader("🤔 What is 'Veber's Rule'?")
    with open("info/14veber.txt", "r") as file:
        veber = file.read()
    st.markdown(veber, unsafe_allow_html=True)
    st.markdown("---")
    st.header("Interesting examples to try")
    st.subheader("💊Ibuprofen:")
    st.write("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
    st.subheader("💊Paracetamol:")
    st.write("CC(=O)NC1=CC=C(C=C1)O")
    
    st.caption("Have fun with chemoinformatics!")

st.title("Metricular: A pocket molecular calculator")
st.caption("Made by 'ares to cat'")

smiles1 = st.text_input("**Please enter a valid SMILES:**")

if smiles1:
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Molecular Descriptors")
        mw = compute_mw(smiles1)
        st.write("Molecular Weight:", mw)
        logp = compute_logp(smiles1)
        st.write("logP:", logp)
        tpsa = compute_tpsa(smiles1)
        st.write("TPSA:", tpsa)
        rotbond = compute_rotbond(smiles1)
        st.write("Number of Rotatable Bonds:", rotbond)
        hbd = compute_hbd(smiles1)
        st.write("H-Bond Donors:", hbd)
        hba = compute_hba(smiles1)
        st.write("H-Bond Acceptors:", hba)
        
    with col2:
        st.subheader("Quantitative Druglikeness")
        qed = compute_qed(smiles1)
        st.write("QED:", qed)
        st.subheader("Molecular Similarity")
        smiles2 = st.text_input("**Enter another SMILES for comparison:**")
        if smiles1 and smiles2:
            sim = compute_similarity(smiles1, smiles2)
            st.write("Tanimoto Index:", sim)
            

if smiles1:
    rules = compute_rules(smiles1)
    if not isinstance(rules, dict):
        st.error(rules if isinstance(rules, str) else "An unknown error occurred. 😿")
    elif rules.get("This compound passes Lipinski's Rule of 5 and Veber's Rule"):
        st.success("This molecule follows both Lipinski's and Veber's Rules! 😺")
    else:
        st.warning("This molecule does not follow Lipinski's and Verber's Rules. 😿")

st.markdown("---")
