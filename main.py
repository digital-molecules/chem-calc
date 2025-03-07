import streamlit as st
from chem_calc import compute_mw, compute_logp, compute_tpsa, compute_rotbond, compute_hbd, compute_hba, compute_qed, compute_lip, molecular_formula, compute_similarity


with st.sidebar:
    st.title("QSAR")
    st.subheader("🤔 What is QSAR and why is it relevant to cheminformatics?")
    st.write("💡 Quantitative structure-activity relationship (QSAR) is a computational method which helps us determine the relationship between chemical structures and biological activity. It is based upon the idea that changes in bioactivity are connected to structural and molecular changes in a group of chemicals. [[1]](https://www.sciencedirect.com/science/article/abs/pii/B9780443186387000025)")
    st.subheader("🤔 What are the parameters of QSAR?")
    st.write("💡 There are many physicochemical parameters that we can take into consideration when predicting the behaviour of a substance within a biological system. In this app, we calculate only some of them, like lipophilicity (logP), but we could also calculate the dipole moment, the pKa, the LUMO, the HOMO, and many more. [[2]](https://www.frontiersin.org/journals/drug-discovery/articles/10.3389/fddsv.2024.1424402/full)")
    st.subheader("🤔 What is the SAR paradox?")
    st.write("💡 The reliability of QSAR methods have been challenged for well over a decade [[3]](https://pubs.acs.org/doi/10.1021/jm020155c). The SAR paradox refers to the fact that structurally similar molecules may, in fact, not have similar biological properties. Thus, it is important to develop dynamic and adaptable models that can be trained with experimental data, using methods such as machine learning. [[4]](https://pmc.ncbi.nlm.nih.gov/articles/PMC6270197/)")

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
        mw = compute_mw(smiles1)
        st.write("MW:", mw)
        logp = compute_logp(smiles1)
        st.write("LogP:", logp)
        tpsa = compute_tpsa(smiles1)
        st.write("TPSA:", tpsa)
        rotbond = compute_rotbond(smiles1)
        st.write("Number of Rotatable Bonds:", rotbond)
        hbd = compute_hbd(smiles1)
        st.write("H-Bond Donors:", hbd)
        hba = compute_hba(smiles1)
        st.write("H-Bond Acceptors:", hba)
        qed = compute_qed(smiles1)
        st.write("QED Drug-Likeness:", qed)

    with col2:
        st.subheader("Molecular Formula")
        formula = molecular_formula(smiles1)
        st.write(f"Molecular Formula: {formula}")
        st.subheader("Tanimoto Similarity")
        smiles2 = st.text_input("**Enter another SMILES for similarity comparison:**")
        if smiles1 and smiles2:
            sim = compute_similarity(smiles1, smiles2)
            st.write("Tanimoto Index:", sim)
            

if smiles1:
    lip_info = compute_lip(smiles1)
    if not isinstance(lip_info, dict):  # If it's not a dictionary, assume it's an error message
        st.error(lip_info if isinstance(lip_info, str) else "An unknown error occurred. 😿")
    elif lip_info.get("This compound passes Lipinski's Rule of 5 and Verber's Rule"):
        st.success("This molecule follows Lipinski's Rule of 5 and Verber's Rule! 😺")
    else:
        st.warning("This molecule does **NOT** follow EITHER Lipinski's Rule of 5 OR Verber's Rule. 😿")

st.markdown("---")
st.caption("<h1 style='text-align: center; font-size: 20px; font-weight: 100; font-style: italic;'>Make sure to refer to the sidebar for more information</h1>", unsafe_allow_html=True)


    #img = render_molecule_image(smiles1)
    #if img:
        #st.image(img, caption="Molecule Structure")


    
#if smiles1:
    #funct = detect_functional_groups(smiles1)
    #st.write("Functional Group(s):", funct)
