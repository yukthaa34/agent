
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

mol = Chem.MolFromSmiles("C(C1C(C(C(C(O1)O)O)O)O)O")
mol_h = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
view = py3Dmol.view()
view.addModel(Chem.MolToMolBlock(mol_h),'mol')
view.setStyle({'stick':{}}).setBackgroundColor('0xeeeeee').zoomTo().show()
