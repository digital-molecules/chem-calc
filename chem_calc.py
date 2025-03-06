from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, DataStructs, QED
#from rdkit.Chem.Fragments import fr_Al_OH, fr_ketone, fr_amide, fr_aldehyde, fr_COO, fr_ester, fr_ether, fr_nitrile, fr_halogen, fr_phenol

def smiles_to_mol(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    return mol

def compute_properties(smiles: str):
    mol = smiles_to_mol(smiles)
    if mol is None:
        return ("Invalid compound 😿")

    properties = {
        "Molecular Weight": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "Number of Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        "Ring Count": rdMolDescriptors.CalcNumRings(mol),
        "H-Bond Donors": rdMolDescriptors.CalcNumHBD(mol),
        "H-Bond Acceptors": rdMolDescriptors.CalcNumHBA(mol),
        "Synthetic Accessibility": QED.qed(mol),
    }
    return properties

def qed(smiles: str):
    mol = smiles_to_mol(smiles)
    if mol is None:
        return ("Invalid compound 😿")

    molecular_weight = Descriptors.MolWt(mol)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    h_bond_donors = rdMolDescriptors.CalcNumHBD(mol)
    h_bond_acceptors = rdMolDescriptors.CalcNumHBA(mol)
    log_p = Descriptors.MolLogP(mol)

    successful_parameters = (
        molecular_weight < 400 and
        ring_count > 0 and
        rotatable_bonds < 5 and
        h_bond_donors <= 5 and
        h_bond_acceptors <= 10 and
        log_p < 5
    )

    return {"This compound's QED score is close to 1, it passes the drug-likeness test": successful_parameters}


def compute_similarity(smiles1: str, smiles2: str):
    mol1 = smiles_to_mol(smiles1)
    mol2 = smiles_to_mol(smiles2)

    if mol1 is None or mol2 is None:
        return ("Invalid compound 😿")

    fp1 = Chem.RDKFingerprint(mol1) #smiles fingerprints
    fp2 = Chem.RDKFingerprint(mol2) #takes both 'canonical' and 'isomeric' smiles into account

    similarity = DataStructs.FingerprintSimilarity(fp1, fp2)
    return {"Tanimoto Similarity": similarity}










#removed function from the streamlit app, kept for reference purposes
def detect_functional_groups(smiles: str):
    mol = smiles_to_mol(smiles)
    if mol is None:
        return ("Invalid compound 😿")

    functional_groups = {
        "Alcohol": fr_Al_OH(mol),
        "Ketone": fr_ketone(mol),
        "Amide": fr_amide(mol),
        "Aldehyde": fr_aldehyde(mol),
        "Carboxylic Acid": fr_COO(mol),
        "Ester": fr_ester(mol),
        "Ether": fr_ether(mol),
        "Nitrile": fr_nitrile(mol),
        "Halogen": fr_halogen(mol),
        "Phenol": fr_phenol(mol),
    }
    return functional_groups
