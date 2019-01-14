from app import db, create_app
from app.db import Atoms

from ase.io import read

if __name__ == '__main__':
    app = create_app()

    traj = read('../data/bcc_bulk_54_expanded_2_high.xyz')
    print(traj)
    # data = {"key1": "value1", "key2": "value2"}

    u = Atoms(numbers={"key1": "value1", "key2": "value2"})

    with app.app_context():
        db.session.add(u)
        db.session.commit()

    with app.app_context():
        atoms = Atoms.query.all()

    print(atoms)

    for at in atoms:
        print(at.id, at.numbers)
