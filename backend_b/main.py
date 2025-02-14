from flask import Flask, request, jsonify
from flask_cors import CORS
import re
from dotenv import load_dotenv
from langchain_openai import AzureChatOpenAI
from langchain.prompts import PromptTemplate
import os
from utils import Utils
import socket
import re
import threading
from openai import OpenAI
import uuid
from gradio_client import Client, handle_file
load_dotenv('.env')


app = Flask(__name__)
CORS(app)





uti = Utils()

llm = AzureChatOpenAI(
        azure_deployment="gpt-4o",
        api_version=os.environ["OPENAI_API_VERSION"],
        temperature=0.15,
        max_tokens=None,
        timeout=None,
        max_retries=2,
        api_key=os.environ["AZURE_OPENAI_KEY"]

    )


def generate_manim_code(question):
    sys_template = """You are a helpful manim coder that first writes all steps involved in creating the video to clarify the user query. steps to be followed are as follows:\n 1)Write down the each and every step involved in solving the problem.\n 2) then based on the steps write a simple manim code step by step without missing any crucial steps like avoiding overlap of multiple texts, missing componeents, wrong calculations etc.. 3) Finally go through the code once again making sure the logic holds and all the rules are followed. 4) Always make sure the code is ready to run state meaning that, do not ommit any lines for breivity. 4) If the process to the solution is lengthy, then make the video also suitably big or small."""
    q_template = """{question}"""

    q_prompt = PromptTemplate(
        template=q_template,
        input_variables=["question"]
    )

    q_prompt = q_prompt.invoke({"question":question})
    voice_over_temp = """
    from manim import *
    from manim_voiceover import VoiceoverScene
    from manim_voiceover.services.gtts import GTTSService


    # Simply inherit from VoiceoverScene instead of Scene to get all the
    # voiceover functionality.
    class RecorderExample(VoiceoverScene):
        def construct(self):
            # You can choose from a multitude of TTS services,
            # or in this example, record your own voice:
            self.set_speech_service(GTTSService())

            circle = Circle()

            # Surround animation sections with with-statements:
            with self.voiceover(text="This circle is drawn as I speak.") as tracker:
                self.play(Create(circle), run_time=tracker.duration)
                # The duration of the animation is received from the audio file
                # and passed to the tracker automatically.

            # This part will not start playing until the previous voiceover is finished.
            with self.voiceover(text="Let's shift it to the left 2 units.") as tracker:
                self.play(circle.animate.shift(2 * LEFT), run_time=tracker.duration)
            .....
    """

    model_response = llm.invoke([("system",str(sys_template)), ("user",str(q_prompt))]).content
    model_response = llm.invoke([("system","You are an AI assistant that helps in analyzing and correcting the provided maxim code to make sure it does not have the following problems:\n\n1) avoiding overlap of multiple texts\n2) missing componeents\n3) wrong calculations.4) making adjustments so that visualizations do not exactly represent represent the output, when the output is very large, short, small, big etc.. to visualize.\n\n Once done, provide me the entire code back to the user."+str(model_response))]).content
    model_response = llm.invoke([("system","Now add voice over to each and every step of the code based on the below example.\n\n"+voice_over_temp+"\n\n Now return the newer code with same functionality but with a clean voice over."),("user",model_response)]).content
    return model_response

# Main function
def math(prompt,class_name):
    # Function to save the generated code to a file
    def save_code_to_file(class_name, code):
        filename = f"{class_name}.py"
        with open(filename, 'w') as file:
            file.write(code)
        return filename

    # Function to activate the environment and run the Manim code
    def run_manim_code(path1, class_name):
        activate_env = f"{path1}"
        manim_command = f"manim -ql {class_name}.py {class_name}"   
        full_command = f"{activate_env} && {manim_command} && python E:\\my_one\\temp.py --class_name {class_name}  && {activate_env.replace('activate.bat','deactivate.bat')}"
        os.system(full_command)
        return "DONE"

    def extract_largest_code_block(text):
        code_blocks = re.findall(r'```python(.*?)```', text, re.DOTALL)
        if not code_blocks:
            return None
        return max(code_blocks, key=len)


    def extract_class_names(code):
        class_pattern = re.compile(r'class\s+(\w+)\s*')
        class_names = class_pattern.findall(code)
        return class_names


    # Generate the Manim code
    manim_code = generate_manim_code(prompt) # make the prev code as this function

    code_cleaned = extract_largest_code_block(manim_code)

    # Extract the class name from the generated code
    class_name1 = extract_class_names(code_cleaned)[0]  # You can parse this from the generated code if needed # make sure to use a regex to obtain it from the generated code
    
    # replace the class name in the generated code
    code_cleaned = code_cleaned.replace(class_name1, class_name)
    # Save the generated code to a file
    save_code_to_file(class_name, code_cleaned)

    # Path to the virtual environment
    path1 = r"C:\Users\BHUVAN\manimations\.venv\Scripts\activate.bat" # replace this with the path to your virtual environment of manim (uv based)

    # Run the Manim code
    return run_manim_code(path1, class_name)


def chemistry(question):


    def extract_largest_code_block(text):
        code_blocks = re.findall(r'```python(.*?)```', text, re.DOTALL)
        if not code_blocks:
            return None
        return max(code_blocks, key=len)

    client = OpenAI(
    base_url="https://openrouter.ai/api/v1",
    api_key=os.environ["OPENAI_API_KEY"],
    )

    chem_cont ="""Here’s a compressed version of the 10 unique RDKit use cases, retaining all essential information while being concise and suitable as context for other models:


### **11. 3D Visualization of Molecules**
- **Problem**: Visualize 3D molecular shapes.
- **Solution**: Use `py3Dmol`.
- **Code**:
  ```python
  from rdkit import Chem
  from rdkit.Chem import AllChem
  import py3Dmol
  mol = Chem.MolFromSmiles("C(C1C(C(C(C(O1)O)O)O)O)O")  # Glucose
  mol_h = Chem.AddHs(mol)
  AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
  view = py3Dmol.view()
  view.addModel(Chem.MolToMolBlock(mol_h),'mol')
  view.setStyle({'stick':{}}).setBackgroundColor('0xeeeeee').zoomTo().show()
  ```

---

Each example addresses a distinct problem type, leveraging RDKit’s versatility in cheminformatics, visualization, and data analysis."""

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

    completion = client.chat.completions.create(
    model="google/gemini-2.0-flash-thinking-exp:free",
    messages = [{
    "role":"system",
    "content":"You are an AI assistant that analyzes and corrects the code generated and returns the detailed explanation of the output of the code but not the actual code along with corrected code. Always give me the full code at the top of the response. Also make sure all the ouputs(image/text) of the code is stored in 'chem' folder."
    },
            {
                "role":"user",
                "content":completion
            }]
    ).choices[0].message.content
    
    code = extract_largest_code_block(completion)

    from new2 import Chemistry
    chem = Chemistry()
    op = chem.generate(code)
    for img in os.listdir('chem'):
        if img.endswith('.png'):
            s3 = uti.save_to_s3(f"chem/{img}",f"{img}","image/png")
            # remove the image from local storage
            os.remove(f"chem/{img}")
        if img.endswith('.svg'):
            s3 = uti.save_to_s3(f"chem/{img}",f"{img}","image/svg+xml")
            # remove the image from local storage
            os.remove(f"chem/{img}")
        if img.endswith('.txt'):
            s3 = uti.save_to_s3(f"chem/{img}",f"{img}","text/plain")
            # remove the image from local storage
            os.remove(f"chem/{img}")
    # if len(op)  == 2 :
    return {"text":{"id":op[1],"rd":op[0]}, "img":0}
    # else:
    #     return {"text":op[0], "img":1}
    # return op

    



@app.route('/api/rag', methods=['POST'])
def get_rag():
    client = Client("bhuvanmdev/QA_document")
    # file_name = request.json.get("file_name")
    file_name = request.form.get('file')
    name = file_name.filename
    file_name.save(name)
    msg = request.json.get("msg")
    result1 = client.predict(
		files=[handle_file(name)],
		chunk_size=1500,
		overlap=300,
		api_name="/process_files"
)
    result2 = client.predict(
		question=msg,
		api_name="/query_streaming"
)
    return jsonify({"result1":result1, "result2":result2})



@app.route('/api/text', methods=['GET'])
def get_text():
    return jsonify({'text': 'Hello, World!'})

@app.route('/api/math', methods=['POST'])
def get_math():
    data = request.json
    vid_flag = data.get("vid",1)
    prompt = data.get("text")
    class_name = "a"+str(uuid.uuid4()).replace("-","")[:10]
    if vid_flag == 1:
        def run_main():
            response = math(prompt,class_name)
        thread = threading.Thread(target=run_main)
        thread.start()
        response = { "text": "Generation Scheduled!!" , "uuid":class_name }
    else:
        response = llm.invoke([("system","You are an AI assistant who helps in solving mathematical problems step by step in detail without any mistakes."),("user",prompt)]).content## normal llm
    return jsonify({"text":response})

@app.route('/api/chemistry', methods=['POST'])
def get_chemistry():
    data = request.json
    prompt = data.get("text")
    img = data.get("img")
    if img == 1:
        res = chemistry(prompt)
        return jsonify(res)
    return jsonify({"text":llm.invoke([("system","You are AI assistant who helps in solving chemistry problems step by step in detail without any mistakes."),("user",prompt)]).content})

if __name__ == '__main__':
    socket.setdefaulttimeout(180)
    app.run(port=5000)