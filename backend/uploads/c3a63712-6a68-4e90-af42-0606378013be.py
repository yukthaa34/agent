
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import os

mol = Chem.MolFromSmiles("C(C1C(C(C(C(O1)O)O)O)O)O")
mol_h = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())

view = py3Dmol.view()
view.addModel(Chem.MolToMolBlock(mol_h),'mol')
view.setStyle({'stick':{}}).setBackgroundColor('0xeeeeee').zoomTo()
view.show()

# Create 'chem' folder if it doesn't exist
if not os.path.exists('chem'):
    os.makedirs('chem')

# Save the image to 'chem' folder
image_path = os.path.join('chem', 'glucose_3d.png')
view.render_image(width=400, height=400, filename=image_path)
print(f"3D molecule image saved to: {image_path}")
