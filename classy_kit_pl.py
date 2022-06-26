from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
import concurrent


#### - Concurrent version ####
"""
Call the multiprocessor method via get_descriptors(submited_smiles_from_post_request)

Smiles size (n)  -- Results returned in (s)  -> Chunk size 128
50                      0.04
5000                    0.90
50000                   9.01
500000                  94.34

RDKIT makes an object that is throwing errors when serialized.  For simplicity the object is created, used,
then removed before the json creation in the lipinski tests and the final json serialization.  If the order is changed or
a json tries to searialize the rdkit object it will throw an error.
"""

lipenski_descriptors = {"5": {"mol_wt": 500, "ha_accepters": 10, "ha_donors": 5, "log_p": 5}, 
"3": {"mol_wt": 300, "ha_accepters": 3, "ha_donors": 3, "log_p": 3, "rotatable_bonds": 3}}

def main(smiles_list):
    return Molecule_obj(smiles_list)

class Molecule_obj:
    """In the parallel version visuals and lipenski methods are not created"""
    
    def __init__(self, smiles):
        self.smiles = smiles
        self.descriptors = {}
        self.errors = {}
        self.get_descriptors()
        
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

    def add_error(self, descriptor, message):
        try:
            self.errors
        except AttributeError:
            self.errors = {}
        finally:
            self.errors.update({descriptor: message})

    def try_method(self, methodToRun, args, descriptor_str, update_descriptors=True):
        descriptor_result = None
        try:
            if not args:
                descriptor_result = methodToRun()
            else:
                descriptor_result = methodToRun(args)
            self.update_descriptors(descriptor_str, descriptor_result)

        except Exception as e:
            self.add_error(descriptor_str, str(e))
            self.update_descriptors(descriptor_str, None)

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

    def aromatic_atoms(self):
        self.try_method(self.aromatic_atoms_method, self.descriptors["mol"], "aromatic_atoms")

    def aromatic_atoms_method(self, mol):
        aromatic_atoms = [mol.GetAtomWithIdx(i).GetIsAromatic() for i in range(mol.GetNumAtoms())]
        aa_count = []
        for i in aromatic_atoms:
            if i == True:
                aa_count.append(i)
        sum_aa_count = sum(aa_count)
        return sum_aa_count

    def aromatic_proportion(self):
        aromatic_prop = None
        try:
            aromatic_prop = self.descriptors["aromatic_atoms"] / self.descriptors["heavy_atom_count"]
        except Exception as e:
            self.add_error("aromatic_proportion", str(e))
            self.update_descriptors("aromatic_proportion", None)
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
        """Instead of jsonifying the entire object - this function organizes a dictionary which is later dumped into a json"""
        self.results = {"smiles" : self.smiles, "errors" : self.errors, "descriptors" : self.descriptors }

    def to_json(self):
        """No longer used but could be adapted for the future"""
        return json.dumps(self.results)

def get_descriptors(smiles):

    """Multiprocessing handler for the json - runs 1 less than the max CPUs available"""
    workers = multiprocessing.cpu_count() - 1
    if workers == 0:
        workers = 1
    with concurrent.futures.ProcessPoolExecutor(workers) as executor:
        # Tested with CPU size of 3 and 31, chunk size of 128 was top performer
        chunksize = round(128)
        if chunksize == 0:
            chunksize = (multiprocessing.cpu_count() - 1)
        result_iter = executor.map(main, smiles["smiles"], chunksize = chunksize)

    return [rd_mol.results for rd_mol in result_iter]