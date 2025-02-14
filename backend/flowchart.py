from openai import OpenAI
import json
import os
from dotenv import load_dotenv

load_dotenv()

client = OpenAI(
  base_url="https://openrouter.ai/api/v1",
  api_key="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",
)

class Flowchart:
    def __init__(self):
        pass

    def generate(self, prompt: str):
        try:
            completion = client.chat.completions.create(
                model="google/gemini-2.0-flash-thinking-exp:free",
                messages = [
                    {
                        "role":"system",
                        "content": "You are an expert in flowchart diagrams, Given any Information you are supposed to create mindmaps for the same. The mindmaps should be very detailed, Put all the topics in an array, then return the connections between the topics in the form of a adjacency list. The Flowchart should have minimum 15 nodes Give a proper JSON response which has the following format:\n\n{\n  'topics': ['topic1', 'topic2', 'topic3'],\n  'connections': [\n    ['topic1', 'topic2'],\n    ['topic2', 'topic3']\n  ]\n}"
                    },
                    {
                        "role": "user",
                        "content": prompt
                    }
                ]
            )
            response = completion.choices[0].message.content[8:-3]
            response = response.replace("'", "\"")

            completion = client.chat.completions.create(
                model="google/gemini-2.0-flash-thinking-exp:free",
                messages = [
                    # {
                    #     "role":"system",
                    #     "content": "Given the following flowchart explain each and every component in detail as the flowchart proceeds and aims to solve the user's query. The explanation should be detailed but concise."
                    # },
                    {
                        "role": "user",
                        "content": prompt + completion.choices[0].message.content + "Given the following flowchart explain each and every component in detail as the flowchart proceeds and aims to solve the user's query. The explanation should be detailed but concise."
                    }
                ]
            )
            response_message = completion.choices[0].message.content

            return (json.loads(response), response_message)
        except Exception as e:
            return str(e), "An error occurred while generating the flowchart"
    
# se = Flowchart()
# flowchart = se.generate("Draw a flowchart for Types of Machine Learning")