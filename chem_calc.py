from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, DataStructs, QED

def smiles_to_mol(smiles1: str):
    mol = Chem.MolFromSmiles(smiles1)
    return mol

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
    return Descriptors.NumRotatableBonds(mol)

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

def compute_rules(smiles1: str):
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

    return {"This compound passes Lipinski's Rule of 5 and Veber's Rule": successful_parameters}


def compute_similarity(smiles1: str, smiles2: str):
    mol1 = smiles_to_mol(smiles1)
    mol2 = smiles_to_mol(smiles2)

    if mol1 is None or mol2 is None:
        return ("Invalid compound 😿")

    fp1 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=2048) #morgan fingerprints
    fp2 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=2048)

    return DataStructs.FingerprintSimilarity(fp1, fp2)


#not implemented on the app
def molecular_formula(smiles1):
    mol = Chem.MolFromSmiles(smiles1)
    if mol is None:
        return "Invalid compound 😿"
    return Chem.rdMolDescriptors.CalcMolFormula(mol)
