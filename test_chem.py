from chem_utils import compute_properties, detect_functional_groups, compute_similarity

#example for verification
smiles1 = "CCO"  # Ethanol
smiles2 = "CCCO"  # Propanol

#basic properties
props = compute_properties(smiles1)
print("Ethanol Properties:", props)

#funct groups
funct = detect_functional_groups(smiles1)
print("Ethanol Properties:", funct)

#similarity
sim = compute_similarity(smiles1, smiles2)
print("Ethanol vs Propanol Similarity:", sim)
