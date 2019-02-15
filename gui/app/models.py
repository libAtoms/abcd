# from mongoengine import Document, EmbeddedDocument, fields
from app.db import db

from mongoengine import Document, EmbeddedDocument

from mongoengine import fields


class Databases(Document):
    meta = {'collection': 'databases'}
    name = fields.StringField(required=True)
    description = fields.StringField(max_length=300)
    tags = fields.ListField(fields.StringField(max_length=30))
