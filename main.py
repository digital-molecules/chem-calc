import streamlit as st
from chem_calc import compute_mw, compute_logp, compute_tpsa, compute_rotbond, compute_hbd, compute_hba, compute_qed, compute_rules, compute_similarity


with st.sidebar:
    st.title("QSAR")
    st.subheader("🤔 What is QSAR and why is it relevant to cheminformatics?")
    st.write("💡 Quantitative structure-activity relationship (QSAR) is the main pillar of ligand-based drug design (LBDD). It is a computational method which helps us determine the relationship between chemical structures and biological activity. QSAR is based upon the idea that changes in bioactivity are connected to structural and molecular changes in a group of chemicals. [[1]](https://www.sciencedirect.com/science/article/abs/pii/B9780443186387000025)")
    st.subheader("🤔 What are the parameters of QSAR?")
    with open("info/1qsar.txt", "r") as file:
        qsar = file.read()
    st.markdown(qsar, unsafe_allow_html=True)
    st.markdown("💡 There are many parameters, also known as descriptors, that we can take into consideration when predicting the behaviour of a substance within a biological system. In this app, we calculate only some of them, like lipophilicity (logP), but we could also calculate other physicochemical descriptors like the dipole moment and the pKa, or electronic descriptors like the LUMO and the HOMO, and many more. [[2]](https://www.frontiersin.org/journals/drug-discovery/articles/10.3389/fddsv.2024.1424402/full)")
    st.subheader("🤔 What is the SAR paradox?")
    st.write("💡 The reliability of QSAR methods has been challenged for well over a decade [[3]](https://pubs.acs.org/doi/10.1021/jm020155c). The SAR paradox refers to how structurally similar molecules may, in fact, not have similar biological properties. Thus, it is important to develop dynamic and adaptable models that can be trained with experimental data, using methods such as machine learning. [[4]](https://pmc.ncbi.nlm.nih.gov/articles/PMC6270197/)")
    st.markdown("---")
    st.header("LogP")
    st.subheader("🤔 What is logP and why is it important?")
    st.write("💡 _Lipophilicity_ represents the compound's preference for lipid phases over aqueous ones. It's a physicochemical property and it plays a major role in a molecule's absorption, distribution within the body, metabolism and excretion (ADME properties), as well as its ability to penetrate biological membranes and barriers. [[5]](https://www.acdlabs.com/wp-content/uploads/download/app/physchem/making_sense.pdf) It is mathematically defined by the partition coefficient P, the measurement of the differential solubility of a compound in two immiscible solvents. *LogP* describes lipophilicity in neutral compounds only, so we should be careful to measure the logP only at a pH where the molecule is in its neutral form and not ionized. [[6]](https://www.acdlabs.com/wp-content/uploads/download/app/physchem/logp_vs_logd.pdf)")
    st.subheader("🤔 How do we interpret LogP values?")
    st.write("💡 Higher logP means higher lipophilicity. A negative logP value indicates that the compound is hydrophilic. A logP of zero indicates the molecule is equally dispersed among the two phases. [[5]](https://www.acdlabs.com/wp-content/uploads/download/app/physchem/making_sense.pdf)")
    st.markdown("---")
    st.header("TPSA")
    st.subheader("🤔 What is TPSA and why is it important?")
    st.write("💡 Polar Surface Area (PSA) is the sum of the contributions to the molecular van der Waals surface area of polar atoms. _Topological Polar Surface Area (TPSA)_ is the simplified measurement of PSA, which helps us reduce the required computational power by avoiding the necessity to calculate ligand 3D structures or biological conformations. TPSA is an important topological descriptor for large and diverse pharmaceutical data sets. Just like logP, it allows us to make ADME predictions. [[7]](https://pmc.ncbi.nlm.nih.gov/articles/PMC7549127/)")
    st.subheader("🤔 How do we interpret TPSA values?")
    st.write("💡 In general, a drug presumed for intestinal absorption should have a TPSA of less than 140Å² and less than 90Å² if it's intended for blood-brain barrier (BBB) penetration. [[8]](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5146717)")
    st.markdown("---")    
    st.header("QED")
    st.subheader("🤔 What is QED and why is it important?")
    st.write("💡 _Quantitative Estimate of Druglikeness (QED)_ is based on 'desirability functions', which are used for multi-criteria optimizations. QED  shows us a deterministic score depending on how closely a molecule aligns with the properties of known oral pharmaceuticals. It is based on the following parameters: molecular weight (MW), lipophilicity (logP), number of hydrogen bond donors (HBD), number of hydrogen bond acceptors (HBA), molecular polar surface area (PSA), number of rotatable bonds (ROTB), number of aromatic rings (AROM) and number of structural alerts (ALERTS). [[9]](https://pmc.ncbi.nlm.nih.gov/articles/PMC3524573/)")
    st.subheader("🤔 How do we interpret QED values?")
    st.write("💡 QED gives us a score from 1 to 0; the closer to 1, the higher the 'drug-likeness'. [[10]](https://academic.oup.com/bib/article/25/4/bbae321/7709089)")
    st.markdown("---")    
    st.header("Tanimoto Index")
    st.subheader("🤔 What is Tanimoto Index and why is it useful?")
    st.write("💡 A molecular fingerprint is structural descriptor and a powerful and compact representation of a molecule as a binary vector, describing the presence or absence of a specific chemical property or structure within the compound. [[11]](https://onlinelibrary.wiley.com/doi/10.1002/minf.201900130) The *Tanimoto Index*, or Jaccard coefficient, divides the shared features between two molecules to the total number of unique features present in either molecule. [[12]](https://pmc.ncbi.nlm.nih.gov/articles/PMC8479812/)")
    st.subheader("🤔 What are Morgan Fingerprints?")
    st.write("💡 _Morgan Fingerprints_ are based on the Morgan algorithm and were created to solve the molecular isomorphism problem — to identify when two molecules, with different atom numberings, are the same. Individual numerical identifiers are assigned to each of the atoms in a molecule, providing more accurate comparisons compared to Daylight fingerprints. [[13]](https://pubs.acs.org/doi/abs/10.1021/ci100050t)")
    st.subheader("🤔 How do we interpret the results of the comparison?")
    st.write("💡 The Tanimoto Index gives a score from 0 to 1, the higher the index the more similar the molecules are and thus more likely to have similar properties.")
    st.markdown("---")    
    st.header("Lipinski's Rule of 5")
    st.subheader("🤔 What is 'Lipinski's Rule of 5'?")
    st.write("💡 '_Lipinski's Rule of 5'_ says, in regards to a drug's oral bioavailability, that '**poor** absorption and permeation are more likely when there are more than 5 H-bond donors, the molecular weight is over 500Da, the logP is over 5 and there are more than 10 H-bond acceptors' and that 'compound classes that are substrates for biological transporters are exceptions to the rule.' [[14]](https://www.sciencedirect.com/science/article/abs/pii/S0169409X96004231) We should be careful when using it, as it is not intended to exclude molecules from being potential orally administered drug candidates, but rather to serve as a guideline for the design of such pharmaceuticals.")
    st.markdown("---")
    st.header("Veber's rule")
    st.subheader("🤔 What is 'Veber's Rule'?")
    st.write("💡 _'Veber's Rule'_ proposes that a compound is more likely to have sufficient oral bioavailability when it has 10 or fewer rotatable bonds and PSA of less than or equal to 140Å². [[15]](https://pubs.acs.org/doi/10.1021/jm020017n) It's often used together with Lipinski's Rule to develop molecules with more 'drug-like' traits.")
    st.markdown("---")
    st.header("Interesting examples to try")
    st.subheader("💊Ibuprofen:")
    st.write("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
    st.subheader("💊Paracetamol:")
    st.write("CC(=O)NC1=CC=C(C=C1)O")
    
    st.caption("Have fun with cheminformatics!")

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
