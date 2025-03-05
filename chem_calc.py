from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, DataStructs, QED
from rdkit.Chem.Fragments import fr_Al_OH, fr_ketone, fr_amide, fr_aldehyde, fr_COO, fr_ester, fr_ether, fr_nitrile, fr_halogen, fr_phenol

def smiles_to_mol(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    return mol


def compute_properties(smiles: str):
    mol = smiles_to_mol(smiles)
    if mol is None:
        return ("Make sure you enter a valid SMILES notation")

    properties = {
        "Molecular Weight": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "Number of Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        "Ring Count": rdMolDescriptors.CalcNumRings(mol)
        "H-Bond Donors": rdMolDescriptors.CalcNumHBD(mol),
        "H-Bond Acceptors": rdMolDescriptors.CalcNumHBA(mol),
        "Synthetic Accessibility": QED.qed(mol),
    }
    return properties

def detect_functional_groups(smiles: str):
    mol = smiles_to_mol(smiles)
    if mol is None:
        return ("Make sure you enter a valid SMILES notation")

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

def compute_similarity(smiles1: str, smiles2: str):
    mol1 = smiles_to_mol(smiles1)
    mol2 = smiles_to_mol(smiles2)

    if mol1 is None or mol2 is None:
        return ("Make sure you enter a valid SMILES notation")

    fp1 = Chem.RDKFingerprint(mol1)
    fp2 = Chem.RDKFingerprint(mol2)

    similarity = DataStructs.FingerprintSimilarity(fp1, fp2)
    return {"Tanimoto Similarity": similarity}
