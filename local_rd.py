from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, AllChem, Draw
import rdkit


smile = 'CC1(CCC(=C(C1)C2=CC=C(C=C2)Cl)CN3CCN(CC3)C4=CC(=C(C=C4)C(=O)NS(=O)(=O)C5=CC(=C(C=C5)NCC6CCOCC6)[N+](=O)[O-])OC7=CN=C8C(=C7)C=CN8)C'

def main():
    mol_descriptors = {}
    #Mol
    m = Chem.MolFromSmiles(smile)
    #print(f"RDKIT Obj: {m}")

    #LogP
    mol_descriptors.update({"log_p": Crippen.MolLogP(m)})
    #print(f"LogP: {log_p}")

    #Mol Wt
    mol_descriptors.update({"mol_wt": Descriptors.ExactMolWt(m)})
    #print(f"Molecular Weight: {mol_wt}")

    #NumHAcceptors
    mol_descriptors.update({"ha_accepters": Lipinski.NumHAcceptors(m)})
    #print(f"Lipinkski # H's acceptors: {ha_accepters}")

    #NumHDoners
    mol_descriptors.update({"ha_doners": Lipinski.NumHDonors(m)})
    #print(f"Lipinkski # H's doners: {ha_doners}")

    #MOL

    print(mol_descriptors)

    if Lipinski_five(mol_descriptors):
        print(f"Meets Lipinski Rule of 5")
    print(f"Does not meet Lipinski Rule of 5")

    if Lipinski_three(mol_descriptors):
        print(f"Meets Lipinski Rule of 3")
    print(f"Does not meet Lipinski Rule of 3")

    contour(m)
    two_d(m)

def Lipinski_five(mol_descriptors):
    lip_five = {"lipinski_five": {"Result": None, "Reason": None }}
    lipenski_descriptors = {"mol_wt": 500, "ha_accepters": 10, "mol_descriptors": 5, "log_p": 5}

    # Reason = []
    # for descriptors in lipenski_descriptors:

    if mol_descriptors['mol_wt'] > 500 and mol_descriptors['ha_accepters'] >= 10 and mol_descriptors['ha_doners'] <=5 and mol_descriptors['log_p'] <= 5:
        return True
    return False

def Lipinski_three(mol_descriptors):
    if mol_descriptors['mol_wt'] > 300 and mol_descriptors['ha_accepters'] >= 3 and mol_descriptors['ha_doners'] <=3 and mol_descriptors['log_p'] <= 3:
        return True
    return False


def mol_to_molblock(mol):
    m2 = Chem.MolToMolBlock(mol)
    return m2

def contour(mol):
    contribs = Chem.rdMolDescriptors._CalcCrippenContribs(mol)
    fig = rdkit.Chem.Draw.SimilarityMaps.GetSimilarityMapFromWeights(mol,[x for x,y in contribs], colorMap='jet', contourLines=10)

def two_d(mol):
    this = Draw.MolToImage(mol)
    this = Draw.MolToFile(mol, filename='test.png', imageType='png')


main()