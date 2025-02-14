from motor.motor_asyncio import AsyncIOMotorClient
from typing import Optional, List, Dict, Any
from bson import ObjectId
import logging
from datetime import datetime

class Database:
    def __init__(self, uri: str, database: str):
        self.client = AsyncIOMotorClient(uri)
        self.db = self.client[database]
        self.messages = self.db.messages
        self.chats = self.db.chats

    async def get_chat_history(self, chat_id: str) -> Optional[List[Dict[str, Any]]]:
        try:
            cursor = self.messages.find({'chatId': chat_id})
            return await cursor.to_list(None)
        except Exception as e:
            logging.error(f"Failed to get chat history: {e}")
            return None

    async def save_message(self, message: Dict[str, Any]) -> Optional[ObjectId]:
        try:
            result = await self.messages.insert_one(message)
            return result.inserted_id
        except Exception as e:
            logging.error(f"Failed to save message: {e}")
            return None

    async def create_chat(self, initial_data: Dict[str, Any] = None) -> Optional[ObjectId]:
        try:
            chat_doc = initial_data if initial_data else {'created_at': datetime.utcnow()}
            result = await self.chats.insert_one(chat_doc)
            return result.inserted_id
        except Exception as e:
            logging.error(f"Failed to create chat: {e}")
            return None

    async def get_all_chats(self) -> Optional[List[Dict[str, Any]]]:
        try:
            return await self.chats.find().to_list(None)
        except Exception as e:
            logging.error(f"Failed to get chats: {e}")
            return None
            
    async def close(self):
        self.client.close()