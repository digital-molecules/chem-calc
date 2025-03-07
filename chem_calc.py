from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, DataStructs, QED
#from rdkit.Chem.Fragments import fr_Al_OH, fr_ketone, fr_amide, fr_aldehyde, fr_COO, fr_ester, fr_ether, fr_nitrile, fr_halogen, fr_phenol

def smiles_to_mol(smiles1: str):
    mol = Chem.MolFromSmiles(smiles1)
    return mol

def compute_properties(smiles1: str):
    mol = smiles_to_mol(smiles1)
    if mol is None:
        return ("Invalid compound 😿")

    properties = {
        "Molecular Weight": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "Number of Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        "H-Bond Donors": rdMolDescriptors.CalcNumHBD(mol),
        "H-Bond Acceptors": rdMolDescriptors.CalcNumHBA(mol),
    }
    return properties

def compute_mw(smiles1: str):
    mol = smiles_to_mol(smiles1)
    if mol is None:
        return ("Invalid compound 😿")
    return Descriptors.MolWt(mol)

def compute_logp(smiles1: str):
    mol = smiles_to_mol(smiles1)
    if mol is None:
        return ("Invalid compound 😿")
    return Descriptors.MolLogP(mol)

def compute_tpsa(smiles1: str):
    mol = smiles_to_mol(smiles1)
    if mol is None:
        return ("Invalid compound 😿")
    return rdMolDescriptors.CalcTPSA(mol)

def compute_rotbond(smiles1: str):
    mol = smiles_to_mol(smiles1)
    if mol is None:
        return ("Invalid compound 😿")
    return Descriptors.NumRotatableBonds

def compute_hbd(smiles1: str):
    mol = smiles_to_mol(smiles1)
    if mol is None:
        return ("Invalid compound 😿")
    return rdMolDescriptors.CalcNumHBD(mol)

def compute_hba(smiles1: str):
    mol = smiles_to_mol(smiles1)
    if mol is None:
        return ("Invalid compound 😿")
    return rdMolDescriptors.CalcNumHBA(mol)

def compute_qed(smiles1: str):
    mol = smiles_to_mol(smiles1)
    if mol is None:
        return ("Invalid compound 😿")
    return QED.qed(mol)

def compute_lip(smiles1: str):
    mol = smiles_to_mol(smiles1)
    if mol is None:
        return ("Invalid compound 😿")

    molecular_weight = Descriptors.MolWt(mol)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    h_bond_donors = rdMolDescriptors.CalcNumHBD(mol)
    h_bond_acceptors = rdMolDescriptors.CalcNumHBA(mol)
    log_p = Descriptors.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)

    successful_parameters = (
        molecular_weight < 500 and
        rotatable_bonds <= 10 and
        h_bond_donors <= 5 and
        h_bond_acceptors <= 10 and
        log_p < 5 and
        tpsa <= 140
    )

    return {"This compound passes Lipinski's Rule of 5 and Verber's Rule": successful_parameters}

def molecular_formula(smiles1):
    mol = Chem.MolFromSmiles(smiles1)
    if mol is None:
        return "Invalid compound 😿"
    return Chem.rdMolDescriptors.CalcMolFormula(mol)


def compute_similarity(smiles1: str, smiles2: str):
    mol1 = smiles_to_mol(smiles1)
    mol2 = smiles_to_mol(smiles2)

    if mol1 is None or mol2 is None:
        return ("Invalid compound 😿")

    fp1 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=2048) #morgan fingerprints
    fp2 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=2048)

    similarity = DataStructs.FingerprintSimilarity(fp1, fp2)
    return {"Tanimoto Similarity": similarity}












#doesnt work :(
def render_molecule_image(smiles: str):
    mol = smiles_to_mol(smiles)
    if mol is None:
        return None
    img = Draw.MolToImage(mol)  # Generate the image without needing GUI
    return img

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
