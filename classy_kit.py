from multiprocessing.sharedctypes import Value
from textwrap import indent
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski
import json

"""
RDKIT makes an object that is throwing errors when serialized.  For simplicity the object is created, used,
then removed before the json creation in the lipinski tests and the final json serialization.  If the order is changed or
a json tries to searialize the rdkit object it will throw an error.
"""

lipenski_descriptors = {"5": {"mol_wt": 500, "ha_accepters": 10, "ha_donors": 5, "log_p": 5}, 
"3": {"mol_wt": 300, "ha_accepters": 3, "ha_donors": 3, "log_p": 3, "rotatable_bonds": 3}}


class Molecule_helper:
    def __init__(self, smiles_list):
        self.results = []
        self.process_samples(smiles_list)
        
    def validate_smiles(self, smile) -> bool:
        """ Currently the only validation step"""
        try:
            Chem.MolFromSmiles(smile)
            return True
        except:
            self.results.append({"smiles" : smile, "errors" : [{"SMILES" : "Could not parse SMILES " }]})
            return False

    def process_samples(self, smiles_list):
        for smile in smiles_list:
            if self.validate_smiles(smile):
                try:
                    self.results.append(Molecule_obj(smile).results)
                except:
                    self.results.append({"smiles" : smile, "errors" : [{"Unknown" : "An error occured that crashed RDKIT " }]})

    def to_json(self):
        return json.dumps(self.results)

    def get_results(self):
        return self.results

class Molecule_obj:
    
    def __init__(self, smiles):
        self.smiles = smiles
        self.errors = {}

        self.descriptors = {}
        self.get_descriptors()
        
        self.Lipinski = {}
        self.get_lipinski_tests()

        self.visualizations = {}
        self.get_visuals()
        
        del self.descriptors["mol"]
        self.prepare_json()

    def get_descriptors(self):
        self.mol()
        self.log_p()
        self.molecular_wt()
        self.ha_accepters()
        self.ha_donors()
        self.rotatable_bonds()
        self.heavy_atom_count()
        self.aromatic_atoms()
        self.aromatic_proportion()
        
    def get_visuals(self):
        self.to_mol_block()

    def get_lipinski_tests(self):
        self.Lepinski_test("5")
        self.Lepinski_test("3")

    def try_method(self, methodToRun, args, descriptor_str, update_descriptors=True):
        descriptor_result = None
        try:
            if not args:
                descriptor_result = methodToRun()
            else:
                descriptor_result = methodToRun(args)
        except Exception as e:
            self.add_error(descriptor_str, str(e))
        finally:
            if update_descriptors:
                self.update_descriptors(descriptor_str, descriptor_result)
            else:
                return descriptor_result

    def add_error(self, descriptor, message):
        try:
            self.errors
        except AttributeError:
            self.errors = {}
        finally:
            self.errors.update({descriptor: message})
            
    def update_descriptors(self, descriptor_str, descriptor):
        self.descriptors.update({descriptor_str: descriptor})  

    def mol(self):
        self.try_method(Chem.MolFromSmiles, self.smiles, "mol")

    def log_p(self):
        self.try_method(Descriptors.MolLogP, self.descriptors["mol"], "log_p")

    def molecular_wt(self):
        self.try_method(Descriptors.ExactMolWt, self.descriptors["mol"], "mol_wt")

    def ha_accepters(self):
        self.try_method(Lipinski.NumHAcceptors, self.descriptors["mol"], "ha_accepters")

    def ha_donors(self):
        self.try_method(Lipinski.NumHDonors, self.descriptors["mol"], "ha_donors")

    def rotatable_bonds(self):
        self.try_method(Lipinski.NumRotatableBonds, self.descriptors["mol"], "rotatable_bonds")

    def heavy_atom_count(self):
        self.try_method(Lipinski.HeavyAtomCount, self.descriptors["mol"], "heavy_atom_count")

    def aromatic_atoms_method(self, mol):
        aromatic_atoms = [mol.GetAtomWithIdx(i).GetIsAromatic() for i in range(mol.GetNumAtoms())]
        aa_count = []
        for i in aromatic_atoms:
            if i == True:
                aa_count.append(i)
        sum_aa_count = sum(aa_count)
        return sum_aa_count

    def aromatic_atoms(self):
        self.try_method(self.aromatic_atoms_method, self.descriptors["mol"], "aromatic_atoms")

    def aromatic_proportion(self):
        aromatic_prop = None
        try:
            aromatic_prop = self.descriptors["aromatic_atoms"] / self.descriptors["heavy_atom_count"]
        except Exception as e:
            self.add_error("aromatic_proportion", str(e))
        finally:
            self.update_descriptors("aromatic_proportion", aromatic_prop)

    def Lepinski_test(self, version):
        local_Lipinski_reasons = []
        try:
            self.Lipinski
        except AttributeError:
            self.Lipinski = {}
        finally:

            for properties in lipenski_descriptors[version]:
                try:
                    if float(lipenski_descriptors[version][properties]) < float(self.descriptors[properties]):
                        local_Lipinski_reasons.append(f'{properties} = {self.descriptors[properties]} above Range of {lipenski_descriptors[version][properties]}')
                except TypeError as e:
                    local_Lipinski_reasons.append(f'{properties} = -> {self.descriptors[properties]} <- cannot be compared again Lipinski value')
                except ValueError as e:
                    local_Lipinski_reasons.append(f'{properties} = -> {self.descriptors[properties]} <- cannot be compared again Lipinski value')

            if local_Lipinski_reasons:
                self.Lipinski.update({version: {"Result": False, "Reasoning": local_Lipinski_reasons}})
            else:
                self.Lipinski.update({version: {"Result": True, "Reasoning": local_Lipinski_reasons}})

    def to_mol_block(self):
        self.visualizations.update({"mol_block" : self.try_method(Chem.MolToMolBlock, self.descriptors["mol"], "mol_block", False)})

    def prepare_json(self):
        self.results = {"smiles" : self.smiles, "errors" : self.errors, "descriptors" : self.descriptors, "Lipinski" : self.Lipinski, "visualizations" : self.visualizations, }

    def to_json(self):
        return json.dumps(self.results)

    

# smile = 'CC1(CCC(=C(C1)C2=CC=C(C=C2)Cl)CN3CCN(CC3)C4=CC(=C(C=C4)C(=O)NS(=O)(=O)C5=CC(=C(C=C5)NCC6CCOCC6)[N+](=O)[O-])OC7=CN=C8C(=C7)C=CN8)C'
# # smile = 'CC1(CCC(=C(C1)C2=CC=C(C=C2)Cl)CN3CCN(CC3)C4=CC(=C(C=C4)C(=O)NS(=O)(=O)C5=CC(=C(C=C5)NCC6CCOCC6)[N+](=O)[O-])OC7=CN=C8C(=C7)C=CN8)Ca'

# # obj = Molecule_obj(smile)
# obj = Molecule_helper([smile])
# print(obj.get_processed_data())
# # #obj.log_p()
# # #print(json.dumps(obj.__dict__, indent=4))
# # print(obj.to_json())