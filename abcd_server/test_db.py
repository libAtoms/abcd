from app import db, create_app
from app.db import Atoms

from ase.io import read

if __name__ == '__main__':
    app = create_app()

    traj = read('../utils/data/bcc_bulk_54_expanded_2_high.xyz')
    print(traj)
    # data = {"key1": "value1", "key2": "value2"}

    u = Atoms(numbers={"key1": "value1", "key2": "value2"})

    # # This will create the database file using SQLAlchemy
    with app.app_context():
        db.drop_all()
        db.create_all()

    with app.app_context():
        db.session.add(u)
        db.session.commit()

    with app.app_context():
        atoms = Atoms.query.all()

    print(atoms)

    for at in atoms:
        print(at.id, at.numbers)
