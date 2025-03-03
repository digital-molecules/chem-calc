from chem_calc import compute_properties, detect_functional_groups, compute_similarity

#the examples used are testing both the case sensitivity and double-character elements
smiles1 = "Cl"  #chlorobenzene
smiles2 = "Clc1ccccc1"  #benzene

props = compute_properties(smiles1)
print("Chlorobenzene Properties:", props)

funct = detect_functional_groups(smiles1)
print("Chlorobenzene's Functional Groups:", funct)

sim = compute_similarity(smiles1, smiles2)
print("Similarity between Chlorobenzene and Benzene:", sim)
