from abcd_server.app.db import db


# from mongoengine import Document, EmbeddedDocument
from mongoengine.fields import EmbeddedDocumentField, ListField, DictField, IntField, FloatField, StringField

# class Atoms(Document):
#     meta = {'strict': False}
#
#     # id = ObjectIdField()
#     arrays = DictField()
#     info = DictField()
#     derived = EmbeddedDocumentField(Derived)
#
#     @classmethod
#     def from_ase(cls, atoms):
#         return cls()


class Databases(db.Document):
    meta = {'collection': 'databases'}
    name = db.StringField(required=True)
    description = db.StringField(max_length=300)
    tags = db.ListField(db.StringField(max_length=30))


class Arrays(db.EmbeddedDocument):
    meta = {'strict': False}
    numbers = db.ListField(db.IntField())
    positions = db.ListField(db.ListField(db.FloatField()))


class Info(db.EmbeddedDocument):
    meta = {'strict': False}


class Derived(db.EmbeddedDocument):
    elements = db.DictField(db.IntField())
    n_atoms = db.IntField(default=111)
    arrays_keys = db.ListField(db.StringField())
    info_keys = db.ListField(db.StringField())


class Atoms(db.Document):
    meta = {'strict': False}

    # id = ObjectIdField()
    arrays = db.EmbeddedDocumentField(Arrays)
    info = db.EmbeddedDocumentField(Info)
    derived = db.EmbeddedDocumentField(Derived)


if __name__ == '__main__':
    sample = {
        "pbc": "alias",
        "Lattice": "alias",
        "cell": "alias",
        "positions": "alias",
        "arrays": {
            "numbers":
                [26, 26, 26, "..."],
            "positions":
                [[-0.10204757, -0.24557776, 0.03216429],
                 [1.23818045, 1.64111918, 1.38883472],
                 [-0.01367307, -0.04099747, 3.0774359],
                 "..."],
            "forces":
                [[0.92008082, 1.39988003, -0.1555232],
                 [1.38380678, -1.65299434, 0.55857061],
                 [-0.39713156, -0.11856745, -1.23895032],
                 "..."]},
        "info": {
            "cell":
                [[8.6368128, 0, 0],
                 [0, 8.6368128, 0],
                 [0, 0, 8.6368128]],
            "pbc": [True, True, True],
            "constraints": [],
            "config_name": "bcc_bulk_54_expanded_2_0013",
            "config_type": "bcc_bulk_54_high",
            "degauss": 0.136056917253,
            "ecutwfc": 1224.51225528,
            "energy": -186883.383046,
            "calculator_name": "unknown",
            "calculator_parameters": {},
            "kpoints": [4, 4, 4]
        }
    }

    from mongoengine import connect

    collection_name = 'atoms'

    db = connect('mongodb://localhost:27017/abcd')
    print(db)
