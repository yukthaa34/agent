import boto3
from botocore.exceptions import ClientError
import os
from dotenv import load_dotenv
import requests
import base64
import traceback
import asyncio

load_dotenv()

class Utils:
    def __init__(self):
        self.client = boto3.client('s3', aws_access_key_id=os.getenv('AWS_ACCESS_KEY_ID'), aws_secret_access_key=os.getenv('AWS_SECRET_ACCESS_KEY'))
        self.bucket = os.getenv('S3_BUCKET')

    def check_creds(self):
        s3 = boto3.resource('s3', aws_access_key_id=os.getenv('AWS_ACCESS_KEY_ID'), aws_secret_access_key=os.getenv('AWS_SECRET_ACCESS_KEY'))
        for bucket in s3.buckets.all():
            print(bucket.name)

    def save_to_s3(self, file_path: str, object_name: str, content_type: str):
        # Save file to S3
        try:
            # response = self.client.put_object({
            #     "Body": file_path, 
            #     "Bucket": self.bucket, 
            #     "Key": f"innerve/{object_name}", 
            #     "ContentType": content_type
            # })
            with open(file_path, "rb") as file_data:
                response = self.client.put_object(
                    Body=file_data, 
                    Bucket=self.bucket, 
                    Key=f"innerve/{object_name}", 
                    ContentType=content_type
                )
            return object_name
        except Exception as e:
            print(e)
            print(traceback.format_exc())
            return None
    
    async def save_to_github(self, file_path: str, object_name: str):
        # Save file to GitHub
        try:
            url = f'https://api.github.com/repos/{os.getenv("GITHUB_USER")}/{os.getenv("GITHUB_REPO")}/contents/{object_name}'
            headers = {
                'Authorization': f'token {os.getenv("GITHUB_TOKEN")}'
            }
            with open(file_path, 'rb') as file:
                content = base64.b64encode(file.read()).decode()
            data = {
                'message': 'Add file',
                'content': content,
                'branch': 'main'
            }
            response = requests.put(url, headers=headers, json=data)
            return response
        except Exception as e:
            print(e)
            print(traceback.format_exc())
            return None
            
        
utils = Utils()
# response = asyncio.run(utils.save_to_github(f"{os.getcwd()}/uploads/test.vcd", 'test.vcd'))
# response = utils.save_to_s3(f"{os.getcwd()}/uploads/rtl.svg", 'rtl.svg', 'image/svg+xml')
# response = utils.check_creds()
# print(response)