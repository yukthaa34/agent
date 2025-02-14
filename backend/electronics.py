import os
import asyncio
from utils import Utils
from openai import OpenAI
import re
import uuid
from dotenv import load_dotenv
import traceback

load_dotenv()

client = OpenAI(
  base_url="https://openrouter.ai/api/v1",
  api_key="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",
)

basepath = '/home/agarwalvivek29/Projects/luminosity/backend/uploads'
class Electronics:
    def __init__(self):
        pass

    async def verilog_execution(self, module: str, module_content:str, testbench: str):
        # Write code to files
        modulepath = f"{basepath}/{module}.v"
        testbenchpath = f"{basepath}/{module}_tb.v"
        outputpath = f"{basepath}/{module}.out"
        simulationpath = f"{basepath}/{module}.vcd"
        yoysyspath = f"{basepath}/{module}.json"
        rtlviewpath = f"{basepath}/{module}.svg"

        errors = []
        testbench = testbench.replace("module.vcd", simulationpath)

        print(module_content)
        print(testbench)

        with open(modulepath, "w") as file:
            file.write(module_content)
        with open(testbenchpath, "w") as file:
            file.write(testbench)

        # Run Icarus Verilog
        process = await asyncio.create_subprocess_exec(
            "iverilog", "-o", outputpath, modulepath, testbenchpath,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        stdout, stderr = await process.communicate()
        if stderr:
            errors.append(f"Failed to run Icarus Verilog: {stderr.decode()}")
        output = stdout.decode()

        # Simulate the Verilog Code
        process = await asyncio.create_subprocess_exec(
            "vvp", outputpath,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        stdout, stderr = await process.communicate()
        if stderr:
            errors.append(f"Failed to simulate Verilog Code: {stderr.decode()}")

        simulation = None
        rtlview = None

        # Upload Simulation Results
        if os.path.exists(simulationpath):
            utils = Utils()
            object_name = f"{module}.vcd"
            file_path = simulationpath
            simulation = await utils.save_to_github(file_path, object_name) 
            if simulation:
                simulation = f"https://vc.drom.io/?github=agarwalvivek29/exec_files/main/{module}.vcd"

        # Create RTL View and Upload
        process = await asyncio.create_subprocess_shell(
            f'yosys -p "read_verilog {modulepath}; proc; opt; write_json {yoysyspath}" && netlistsvg {yoysyspath} -o {rtlviewpath}',
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        stdout, stderr = await process.communicate()

        if stderr:
            errors.append(f"Failed to generate RTL View: {stderr.decode()}")
            print(stderr.decode())

        if os.path.exists(rtlviewpath):
            utils = Utils()
            object_name = f"{module}.svg"
            rtlview = utils.save_to_s3(rtlviewpath, object_name, 'image/svg+xml')
            if rtlview:
                rtlview = f"https://solace-outputs.s3.ap-south-1.amazonaws.com/innerve/{module}.svg"
            if not rtlview:
                errors.append("Failed to upload RTL View to S3")

        if os.path.exists(modulepath):
            os.remove(modulepath)
        if os.path.exists(testbenchpath):
            os.remove(testbenchpath)
        if os.path.exists(outputpath):
            os.remove(outputpath)
        if os.path.exists(simulationpath):
            os.remove(simulationpath)
        if os.path.exists(yoysyspath):
            os.remove(yoysyspath)
        if os.path.exists(rtlviewpath):
            os.remove(rtlviewpath)

        return {
            "uuid": module,
            "success": True,
            "message": output,
            "simulation": simulation,
            "rtlview": rtlview,
            "errors": errors
        }
    
    async def get_verilog_code(self, prompt: str):
        completion = client.chat.completions.create(
            model="google/gemini-2.0-flash-thinking-exp:free",
            messages = [
                {
                    "role":"system",
                    "content": "Give me a verilog code along with its test bench that is also compatible with ICARUS verilog, Don't mention the steps related I am aware, Just explain the module appropriately and the testbench's working. If the given module / Circuit is a standard circuit kindly give the truth table appropriately. Make sure that the waveforms of the testbench are dumped into module.vcd."
                },
                {
                    "role": "user",
                    "content": prompt
                }
            ]
        )
        response = completion.choices[0].message.content

        # Match Verilog code blocks
        pattern = r"```verilog[\s\S]*?```"
        verilog_blocks = re.findall(pattern, response)
        verilog_blocks = sorted(verilog_blocks, key=lambda x: len(x))
        print(verilog_blocks)
        print(len(verilog_blocks))
        if len(verilog_blocks) >= 2:
            exec = await self.verilog_execution(str(uuid.uuid4()), verilog_blocks[-2][11:-3], verilog_blocks[-1][11:-3])
            exec["message"] = response
            return exec
        
    async def obtain_chip_design(self, codepath, imagepath, module):
        errors = []
        
        process = await asyncio.create_subprocess_shell(
            f'kicad-cli pcb export svg -o {imagepath} -l "F.Cu,B.Cu,F.SilkS,B.SilkS" {codepath}',
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        stdout, stderr = await process.communicate()
        if stderr:
            errors.append(f"Failed to generate PCB Image: {stderr.decode()}")

        utils = Utils()
        object_name = f"{module}.svg"
        image = utils.save_to_s3(imagepath, object_name, 'image/svg+xml')

        if image:
            image = f"https://solace-outputs.s3.ap-south-1.amazonaws.com/innerve/{module}.svg"

        if os.path.exists(codepath):
            os.remove(codepath)
        if os.path.exists(imagepath):
            os.remove(imagepath)

        return {
            "uuid": module,
            "success": True,
            "image": image,
            "errors": errors
        }
        
    async def analyse_chip_design(self, module, filepath):

        with open(filepath, "r") as file:
            code = file.read()

        completion = client.chat.completions.create(
            model="google/gemini-2.0-flash-thinking-exp:free",
            messages = [
                { "role": "system", "content": "You are a Semiconductor Chip design expert" },
                { "role": "user", "content": code + "Given the above code, analyse the circuit and suggest inefficiencies and improvements. Try to provide suggestions on creating better circuits and providing better architectures."},
                { "role": "system", "content": "Now create kicad code for a new PCB which addresses the inefficiencies in the design and provides a better solution" }
            ]
        )
        codepath = filepath
        imagepath = f"{basepath}/{module}.svg"
        response = completion.choices[0].message.content

        res = await self.obtain_chip_design(codepath, imagepath, module)
        res["message"] = response
        return res
        
    async def get_matlab_code(self, prompt: str):
        completion = client.chat.completions.create(
            model="google/gemini-2.0-flash-thinking-exp:free",
            messages = [
                {
                    "role":"system",
                    "content": "Given the following prompt, write a MATLAB code that can be used to solve the problem. The code should be well documented and should be able to solve the problem efficiently. The code should be able to handle edge cases and should be able to provide the correct output for the given input.!!Important: Make sure that the code is compatible with octavius and MATLAB. Make it such that all the graph outputs are saved into module.png"
                },
                {
                    "role": "user",
                    "content": prompt
                }
            ]
        )
        response = completion.choices[0].message.content
        pattern = r"```matlab[\s\S]*?```"
        matlab_blocks = re.findall(pattern, response)
        matlab_blocks = sorted(matlab_blocks, key=lambda x: len(x))
        module = str(uuid.uuid4())

        res = self.matlab_execution(module, matlab_blocks[-1][10:-3])
        res["chat"] = response
        return res
    
    async def matlab_execution(self, module, code: str):
        codepath = f"{basepath}/{module}.m"
        outputpath = f"{basepath}/{module}.png"
        errors = []

        code.replace("module.png", outputpath)

        with open(codepath, "w") as file:
            file.write(code)

        process = await asyncio.create_subprocess_shell(
            f"octave --silent --eval 'source {codepath}'", 
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        stdout, stderr = await process.communicate()

        if stderr:
            errors.append(f"Failed to run MATLAB Code: {stderr.decode()}")
        output = stdout.decode()

        image = None
        if os.path.exists(outputpath):
            utils = Utils()
            image = utils.save_to_s3(outputpath, f"{module}.png", 'image/png')

        if os.path.exists(codepath):
            os.remove(codepath)
        if os.path.exists(outputpath):
            os.remove(outputpath)

        return {
            "success": True,
            "message": output,
            "errors": errors,
            "image": image
        }
    
# el = Electronics()
# output = asyncio.run(el.get_verilog_code("Generate Verilog Code for Full Adder Circuit"))