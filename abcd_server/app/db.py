from abcd.backends import MongoDatabase
from flask_mongoengine import MongoEngine

db = MongoEngine()


class Databases(db.Document):
    meta = {'collection': 'databases'}
    name = db.StringField(required=True)
    description = db.StringField(max_length=300)
    tags = db.ListField(db.StringField(max_length=30))
