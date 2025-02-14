
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import os

# Create 'chem' folder if it doesn't exist
if not os.path.exists('chem'):
    os.makedirs('chem')

mol = Chem.MolFromSmiles("C(C1C(C(C(C(O1)O)O)O)O)O")
mol_h = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
view = py3Dmol.view()
view.addModel(Chem.MolToMolBlock(mol_h),'mol')
view.setStyle({'stick':{}}).setBackgroundColor('0xeeeeee').zoomTo()

# Save the visualization to a file in 'chem' folder
image_path = os.path.join('chem', 'glucose_3d.png')
view.renderToFile(image_path, width=400, height=300)

view.show()
