from fastapi import FastAPI, File, UploadFile
from fastapi.middleware.cors import CORSMiddleware
import uvicorn
import os
from electronics import Electronics
from pydantic import BaseModel
import json
from flowchart import Flowchart
from chemistry import Chemistry
import uuid
from dotenv import load_dotenv
import shutil

load_dotenv()

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
def read_root():
    return {"message": "Hello, World!"}

# @app.post("/chat")
# async def create_chat():
#     chat = await databaseops.create_chat()
#     return {
#         "success": True,
#         "chat": chat
#     }

# @app.get("/chats")
# async def get_all_chats():
#     chats = await databaseops.get_all_chats()
#     return {
#         "success": True,
#         "chats": chats
#     }

# @app.get("/chat/{chatId}")
# async def get_chat_history(chatId: str):
#     messages = await databaseops.get_chat_history(chatId)
#     return {
#         "success": True,
#         messages: messages
#     }

# @app.post("/message")
# async def save_message(message: dict):
#     message = await databaseops.save_message(message)
#     return {
#         "success": True,
#         "message": message
#     }

class ElecReq(BaseModel):
    prompt: str

@app.post("/electronics/matlab")
async def run_matlab_code(elec_req: ElecReq):
    electronics = Electronics()
    result = await electronics.get_matlab_code(elec_req.prompt)
    return {
        "success": True,
        "message": result["message"],
        "output": result["output"],
        "errors": json.dumps(result["errors"]),
        "image": result["image"],
        "chat": result["chat"]
    }

@app.post("/electronics/verilog")
async def run_verilog_code(elec_req: ElecReq):
    electronics = Electronics()
    result = await electronics.get_verilog_code(elec_req.prompt)
    return {
        "success": True,
        "message": result["message"],
        "simulation": str(result["simulation"]),
        "rtlview": result["rtlview"],
        "errors": json.dumps(result["errors"])
    }

@app.post("/electronics/chipdesign")
async def run_chip_design(file: UploadFile = File(...)):
    module = str(uuid.uuid4())

    filepath = f"/home/agarwalvivek29/Projects/luminosity/backend/uploads/{module}.{(file.filename).split('.')[-1]}"

    with open(filepath, "wb") as buffer:
        shutil.copyfileobj(file.file, buffer)

    electronics = Electronics()
    result = await electronics.analyse_chip_design(module, filepath)
    return {
        "success": True,
        "message": result["message"],
        "image": result["image"],
        "uuid": result["uuid"],
        "errors": json.dumps(result["errors"])
    }

@app.post("/mindmaps")
async def generate_flowchart(elecreq: ElecReq):
    prompt = elecreq.prompt
    flowchart = Flowchart()
    flowchart, message = flowchart.generate(prompt)
    return {
        "success": True,
        "flowchart": json.dumps(flowchart),
        "message": message
    }

@app.post("/chemistry/molecule")
def generate_molecule(elecreq: ElecReq):
    prompt = elecreq.prompt
    chemistry = Chemistry()
    try:
        molecule, uniqueId = chemistry.generate(prompt)
        return {
            "success": True,
            "molecule": molecule,
            "uniqueId": uniqueId
        }
    except Exception as e:
        return {
            "success": False,
            "error": str(e)
        }