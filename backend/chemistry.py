import os
from openai import OpenAI
import re
import asyncio
import uuid
from jupyter_tt import execute_in_jupyter
from dotenv import load_dotenv

load_dotenv()

client = OpenAI(
    base_url="https://openrouter.ai/api/v1",
    api_key="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",
)

chem_cont ="""Here’s a compressed version of the 10 unique RDKit use cases, retaining all essential information while being concise and suitable as context for other models:
    ---

    ### **1. 2D Molecular Visualization**
    - **Problem**: Visualize 2D molecular structures with atom indices.
    - **Solution**: Use `Draw.MolToImage`.
    - **Code**:
    python
    from rdkit import Chem
    from rdkit.Chem import Draw
    mol = Chem.MolFromSmiles('CCO')  # Ethanol
    Draw.MolToImage(mol, highlightAtoms=[0, 1], legend='Ethanol').save('ethanol_2d.png')
    

    ---

    ### **2. 3D Conformer Generation**
    - **Problem**: Generate and visualize 3D conformers.
    - **Solution**: Use `EmbedMultipleConfs` and `py3Dmol`.
    - **Code**:
    python
    from rdkit.Chem import AllChem
    import py3Dmol
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))  # Ethanol with hydrogens
    AllChem.EmbedMultipleConfs(mol, numConfs=5)
    view = py3Dmol.view(mol.ToXYZ(confId=0))
    view.setStyle({'stick': {}}).show()
    

    ---

    ### **3. Molecular Descriptor Calculation**
    - **Problem**: Compute logP and TPSA.
    - **Solution**: Use `Descriptors` module.
    - **Code**:
    python
    from rdkit.Chem import Descriptors
    mol = Chem.MolFromSmiles('C1=CC=CN=C1')  # Pyridine
    logP, tpsa = Descriptors.MolLogP(mol), Descriptors.TPSA(mol)
    print(f"logP: {logP}, TPSA: {tpsa}")
    

    ---

    ### **4. QSAR Model for Bioactivity**
    - **Problem**: Predict IC50 using molecular descriptors.
    - **Solution**: Train regression models with `scikit-learn`.
    - **Code**:
    python
    from rdkit.ML.Descriptors import MoleculeDescriptors
    from sklearn.linear_model import LinearRegression
    descriptor_calc = MoleculeDescriptors.MolecularDescriptorCalculator(['MolWt', 'NumRotatableBonds'])
    X = [descriptor_calc.CalcDescriptors(mol) for mol in molecules_list]
    model = LinearRegression().fit(X, y_ic50)  # y_ic50: experimental values
    

    ---

    ### **5. Chemical Reaction Simulation**
    - **Problem**: Simulate ester hydrolysis.
    - **Solution**: Use `ReactionFromSmarts`.
    - **Code**:
    python
    rxn = AllChem.ReactionFromSmarts('[C:1](=[O:2])-[O:3]-[C:4]>>[C:1](=[O:2])-[O:3].[C:4][OH]')
    aspirin = Chem.MolFromSmiles('CC(=O)OC1=CC=CC=C1C(=O)O')
    products = rxn.RunReactants((aspirin,))  # Salicylic acid + acetic acid
    

    ---

    ### **6. Molecular Similarity Search**
    - **Problem**: Find similar molecules using fingerprints.
    - **Solution**: Compute Tanimoto similarity.
    - **Code**:
    python
    from rdkit.DataStructs import TanimotoSimilarity
    fp1, fp2 = AllChem.GetMorganFingerprint(mol1, 2), AllChem.GetMorganFingerprint(mol2, 2)
    similarity = TanimotoSimilarity(fp1, fp2)  # Range: 0 (dissimilar) to 1 (identical)
    

    ---

    ### **7. Tautomer Enumeration**
    - **Problem**: Generate all tautomers.
    - **Solution**: Use `EnumerateTautomers`.
    - **Code**:
    python
    from rdkit.Chem.MolStandardize import tautomer
    mol = Chem.MolFromSmiles('C1=CNC=C1')  # Imidazole
    tautomers = tautomer.TautomerEnumerator().Enumerate(mol)
    Draw.MolsToGridImage(tautomers, molsPerRow=4).save('tautomers.png')
    

    ---

    ### **8. 3D Molecular Alignment**
    - **Problem**: Align molecules to a reference scaffold.
    - **Solution**: Use `AlignMol`.
    - **Code**:
    python
    ref_mol = Chem.MolFromSmiles('C1CCCCC1')  # Cyclohexane
    probe_mol = Chem.MolFromSmiles('C1CCCCC1C(=O)O')  # Modified scaffold
    AllChem.AlignMol(probe_mol, ref_mol)  # Align in 3D space
    

    ---

    ### **9. ADMET Property Prediction**
    - **Problem**: Predict solubility (logS).
    - **Solution**: Use descriptors with a pre-trained model.
    - **Code**:
    python
    from rdkit.ML.Descriptors import MoleculeDescriptors
    descriptors = ['MolLogP', 'MolWt', 'NumHDonors']
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptors)
    X = [calc.CalcDescriptors(mol) for mol in dataset]
    logS_predictions = model.predict(X)  # Pre-trained model
    

    ---

    ### **10. Clustering & PCA Visualization**
    - **Problem**: Visualize molecular diversity.
    - **Solution**: Reduce descriptor dimensions with PCA.
    - **Code**:
    python
    from sklearn.decomposition import PCA
    import matplotlib.pyplot as plt
    pca = PCA(n_components=2).fit_transform(descriptor_matrix)
    plt.scatter(pca[:,0], pca[:,1], c=cluster_labels)
    plt.xlabel('PC1'), plt.ylabel('PC2')
    plt.savefig('pca_clusters.png')
    

    ---

    ### **11. 3D Visualization of Molecules**
    - **Problem**: Visualize 3D molecular shapes.
    - **Solution**: Use `py3Dmol`.
    - **Code**:
    python
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import py3Dmol
    mol = Chem.MolFromSmiles("C(C1C(C(C(C(O1)O)O)O)O)O")  # Glucose
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
    view = py3Dmol.view()
    view.addModel(Chem.MolToMolBlock(mol_h),'mol')
    view.setStyle({'stick':{}}).setBackgroundColor('0xeeeeee').zoomTo().show()
    

    ---

    Each example addresses a distinct problem type, leveraging RDKit’s versatility in cheminformatics, visualization, and data analysis."""

class Chemistry:
    def __init__(self):
        pass

    def extract_largest_code_block(self, text):
        code_blocks = re.findall(r'```python(.*?)```', text, re.DOTALL)
        if not code_blocks:
            return None
        return max(code_blocks, key=len)

    def generate(self, question):
        completion = client.chat.completions.create(
            model="google/gemini-2.0-flash-thinking-exp:free",
                messages=[{
                    "role":"system",
                    "content":"You are an rdkit expert that uses the rdkit library to write the necessary code in python to answer the user queries. Make sure all the answers are right and the code should be small and straight to the point.\nMore importantly, generate only required code.\n\n some examples are as follows:\n\n" + chem_cont
                },
                {
                    "role": "user",
                    "content": question
                }
            ]
        ).choices[0].message.content

        # print(completion)

        completion2 = client.chat.completions.create(
        model="google/gemini-2.0-flash-thinking-exp:free",
        messages = [{
                "role":"system",
                "content":"You are an AI assistant that analyzes and corrects the code generated and returns the detailed explanation of the output of the code but not the actual code along with corrected code. Always give me the full code at the top of the response. Also make sure all the ouputs(image/text) of the code is stored in 'chem' folder."
            },
            {
                "role":"user",
                "content": completion
            }]
        ).choices[0].message.content


        code = self.extract_largest_code_block(completion2)
        print(code)

        # class_name= str(uuid.uuid4())
        # filepath = f'/home/agarwalvivek29/Projects/luminosity/backend/uploads/{class_name}.py'
        # with open(filepath, 'w') as file:
        #     file.write(code)

        # out, err = await self.run_manim_code(class_name, code)
        if str(code).find('py3Dmol') != -1:
            code += "\nprint(Chem.MolToMolBlock(mol_h))\nprint(view.uniqueid)"
        print(code)
        out, err = execute_in_jupyter(str(code))
        if str(code).find('py3Dmol') != -1:
            print('Output:', out)
            print('Error:', err)
            a,b = out.split('END\n\n')
            a += "END"
            b = b.strip('\n').strip()
            return (a,b)
            
        print('Output:', out)
        print('Error:', err)