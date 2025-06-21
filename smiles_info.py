import sys
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen, inchi

def get_cid_from_inchikey(inchikey):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"
    r = requests.get(url)
    if r.status_code != 200:
        return None
    try:
        data = r.json()
        cids = data.get("IdentifierList", {}).get("CID", [])
        return cids[0] if cids else None
    except:
        return None

def get_synonyms_from_cid(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
    r = requests.get(url)
    if r.status_code != 200:
        return []
    try:
        data = r.json()
        return data.get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])
    except:
        return []

def get_iupac_from_cid(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName/JSON"
    r = requests.get(url)
    if r.status_code != 200:
        return None
    try:
        data = r.json()
        return data["PropertyTable"]["Properties"][0].get("IUPACName")
    except:
        return None

def print_pubchem_info(iupac_name, synonyms):
    print(f"IUPAC Name: {iupac_name if iupac_name else 'N/A'}")
    if synonyms:
        print("Synonyms:")
        for s in synonyms[:5]:
            print(f"  - {s}")
    else:
        print("Synonyms: N/A")

def print_rdkit_info(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print(f"Invalid SMILES: {smiles}")
        return
    print(f"\n=== RDKit Properties ===")
    print(f"SMILES: {smiles}")
    print(f"Molecular Formula: {rdMolDescriptors.CalcMolFormula(mol)}")
    print(f"Molecular Weight: {Descriptors.MolWt(mol):.2f}")
    print(f"LogP (octanol-water): {Crippen.MolLogP(mol):.2f}")
    print(f"H-Bond Donors: {rdMolDescriptors.CalcNumHBD(mol)}")
    print(f"H-Bond Acceptors: {rdMolDescriptors.CalcNumHBA(mol)}")
    print(f"TPSA: {rdMolDescriptors.CalcTPSA(mol):.2f}")
    print(f"Rings: {rdMolDescriptors.CalcNumRings(mol)}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: python3 {sys.argv[0]} '<SMILES_STRING>'")
        sys.exit(1)

    input_smiles = sys.argv[1]
    mol = Chem.MolFromSmiles(input_smiles)

    if mol is None:
        print(f"Invalid SMILES: {input_smiles}")
        sys.exit(1)

    canonical_smiles = Chem.MolToSmiles(mol)
    inchikey = inchi.MolToInchiKey(mol)

    cid = get_cid_from_inchikey(inchikey)

    if cid:
        iupac = get_iupac_from_cid(cid)
        synonyms = get_synonyms_from_cid(cid)
        print("=== PubChem Information ===")
        print_pubchem_info(iupac, synonyms)
    else:
        print("No PubChem data found.")

    print_rdkit_info(canonical_smiles)
