from flask_sqlalchemy import SQLAlchemy
from datetime import datetime

db = SQLAlchemy()


class Atoms(db.Model):
    __tablename__ = 'atoms'
    id = db.Column(db.Integer, primary_key=True)  # Local database id
    unique_id = db.Column(db.String)  # Globally unique hexadecimal id
    timestamp = db.Column(db.DateTime, index=True, default=datetime.utcnow)
    numbers = db.Column(db.ARRAY(db.Integer))  # Atomic numbers shape: (N,)
    pbc = db.Column(db.Boolean)  # Periodic boundary condition flags shape: (3,)
    cell = db.Column(db.Float)  # Unit cell shape: (3, 3)
    positions = db.Column(db.Float)  # Atomic positions shape: (N, 3)
    forces = db.Column(db.Float)  # Atomic forces shape: (N, 3)
    energy = db.Column(db.Float)  # Total energy
    info = db.column(db.JSON)

    # charges = db.Column(db.Float)  # Atomic charges shape: (N,)
    # ctime = db.Column(db.DateTime)  # Creation time
    # mtime = db.Column(db.DateTime)  # Modification time
    # user = db.Coulumn(db.String)  # User name
    # initial_magmoms = db.Column(db.Float)  # Initial atomic magnetic moments shape: (N,)
    # initial_charges = db.Column(db.Float)  # Initial atomic charges shape: (N,)
    # masses = db.Column(db.Float)  # Atomic masses shape: (N,)
    # tags = db.Column(db.Integer)  # Tags shape: (N,)
    # momenta = db.Column(db.Float)  # Atomic momenta shape: (N, 3)
    # stress = db.Column(db.Float)  # Stress tensor shape: (6,)
    # dipole = db.Column(db.Float)  # Electrical dipole shape: (3,)
    # magmom = db.Column(db.Float)  # Magnetic moment
    # magmoms = db.Column(db.Float)  # Atomic magnetic moments shape: (N,)
    # calculator = db.Coulumn(db.String)  # Calculator name
    # constraints	# Constraints	list of dict
    # calculator_parameters	# Calculator parameters	dict

    def __repr__(self):
        return '<Atoms {}>'.format(self.unique_id)

    @classmethod
    def from_json(cls, s):
        return cls()


class Projects(db.Model):
    __tablename__ = 'projects'
    id = db.Column(db.Integer, primary_key=True)


class User(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(80), index=True, unique=True, nullable=False)
    email = db.Column(db.String(120), index=True, unique=True, nullable=False)
    password_hash = db.Column(db.String(128))

    def __repr__(self):
        # return '<User {}>'.format(self.username)
        return f'<User {self.username}>'
