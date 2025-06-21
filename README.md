# smiles_info
This script takes a SMILES string as input, parses it using RDKit to calculate chemical properties,
and queries PubChem via InChIKey to retrieve the compound's IUPAC name and common synonyms.

Outputs include:
- Canonical SMILES
- Molecular formula, weight, LogP, TPSA, H-bond counts, ring count
- IUPAC name
- Up to 5 known synonyms (e.g., trade names or alternate names)

# Example:
    python smiles_info_with_names.py "CC(=O)NC1=CC=C(C=C1)O"
    python smiles_info_with_names.py "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
    python smiles_info_with_names.py "CC(CC1=CC=CC=C1)CC(C(=O)O)NC(=O)CCCN"

# Dependencies:
    - rdkit
    - requests
