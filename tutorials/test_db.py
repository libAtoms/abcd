from pathlib import Path
from ase.io import iread, read

from abcd import ABCD

if __name__ == "__main__":
    # http requests
    # url = 'http://localhost:5000/api'

    # Pymongo
    # url = 'mongodb://root:example@localhost:27018/?authSource=admin'

    # Mongoengine
    # https://stackoverflow.com/questions/36200288/mongolab-pymongo-connection-error
    url = "mongodb://root:example@localhost:27018/?authSource=admin"

    abcd = ABCD(url, db="abcd", collection="default")
    print(abcd)

    abcd.print_info()

    with abcd as db:
        db.destroy()

    abcd.destroy()

    direcotry = Path("../tutorials/data/")
    file = direcotry / "bcc_bulk_54_expanded_2_high.xyz"
    # file = direcotry / 'GAP_6.xyz'

    traj = read(file.as_posix(), index=slice(None))
    abcd.print_info()

    with abcd as db:
        db.push(traj)

    abcd.print_info()

# from app_old import db, create_app
# from app_old.db import Atoms
#
# from ase.io import read
#
# if __name__ == '__main__':
#     app_old = create_app()
#
#     traj = read('../utils/data/bcc_bulk_54_expanded_2_high.xyz')
#     print(traj)
#     # data = {"key1": "value1", "key2": "value2"}
#
#     u = Atoms(numbers={"key1": "value1", "key2": "value2"})
#
#     # # This will create the database file using SQLAlchemy
#     with app_old.app_context():
#         db.drop_all()
#         db.create_all()
#
#     with app_old.app_context():
#         db.session.add(u)
#         db.session.commit()
#
#     with app_old.app_context():
#         atoms = Atoms.query.all()
#
#     print(atoms)
#
#     for at in atoms:
#         print(at.id, at.numbers)
