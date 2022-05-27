from textwrap import indent
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski
import json

"""
RDKIT makes an object that is throwing errors when serialized.  For simplicity the object is craeted, used,
then removed before the json creation in the lipinski tests and the final json serialization.  If the order is changed or
a json tries to searialize the rdkit object it will throw an error.
"""

lipenski_descriptors = {"5": {"mol_wt": 500, "ha_accepters": 10, "ha_donors": 5, "log_p": 5}, 
"3": {"mol_wt": 300, "ha_accepters": 3, "ha_donors": 3, "log_p": 3, "rotatable_bonds": 3}}

class Molecule_obj:
    
    def __init__(self, smiles):
        self.smiles = smiles
        self.errors = []
        self.get_descriptors()
        self.get_visuals()
        if not self.errors:
            del self.Mol
            self.get_lipinski_tests()
        else:
            del self.Mol

    def get_descriptors(self):
        self.mol()
        self.log_p()
        self.molecular_wt()
        self.ha_accepters()
        self.ha_donors()
        self.rotatable_bonds()
        
    def get_visuals(self):
        self.to_mol_block()

    def get_lipinski_tests(self):
        self.Lepinski_test("5")
        self.Lepinski_test("3")

    def to_mol_block(self):
        try:
            self.mol_block = Chem.MolToMolBlock(self.Mol)
        except:
            self.add_error("mol_block")

    def add_error(self, descriptor):
        try:
            self.errors
        except AttributeError:
            self.errors = []
        finally:
            self.errors.append(descriptor)

    def mol(self):
        try:
            self.Mol = Chem.MolFromSmiles(self.smiles)
        except:
            self.add_error("mol")

    def log_p(self):
        try:
            self.log_p = Crippen.MolLogP(Chem.MolFromSmiles(self.smiles))

        except:
            self.add_error("log_p")

    def molecular_wt(self):
        try:
            self.mol_wt = Descriptors.ExactMolWt(self.Mol)

        except:
            self.add_error("mol_wt")

    def ha_accepters(self):
        try:
            self.ha_accepters = Lipinski.NumHAcceptors(self.Mol)
        except:
            self.add_error("ha_accepters")

    def ha_donors(self):
        try:
            self.ha_donors = Lipinski.NumHDonors(self.Mol)
        except:
            self.add_error("ha_donors")

    def rotatable_bonds(self):
        try:
            self.rotatable_bonds = Lipinski.NumRotatableBonds(self.Mol)
        except:
            self.add_error("rotatable_bonds")

    def Lepinski_test(self, version):
        json_obj = json.loads(json.dumps(self.__dict__, indent=4))
        local_Lipinski_reasons = []

        for properties in lipenski_descriptors[version]:
            if float(lipenski_descriptors[version][properties]) < float(json_obj[properties]):
                local_Lipinski_reasons.append(f'{properties} of {json_obj[properties]} above Range of {lipenski_descriptors[version][properties]}')

        try:
            self.Lipinski
        except AttributeError:
            self.Lipinski = {}
        finally:
            if local_Lipinski_reasons:
                self.Lipinski.update({version: {"Result": False, "Reasoning": local_Lipinski_reasons}})
            else:
                self.Lipinski.update({version: {"Result": True, "Reasoning": local_Lipinski_reasons}})

    def to_json(self):
        return json.dumps(self.__dict__, indent=4)
# smile = 'CC1(CCC(=C(C1)C2=CC=C(C=C2)Cl)CN3CCN(CC3)C4=CC(=C(C=C4)C(=O)NS(=O)(=O)C5=CC(=C(C=C5)NCC6CCOCC6)[N+](=O)[O-])OC7=CN=C8C(=C7)C=CN8)C'


#obj = molecule_obj(smile)
#obj.log_p()
#print(json.dumps(obj.__dict__, indent=4))