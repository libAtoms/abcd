import requests
from ase.io import iread
from pathlib import Path
from abcd import ABCD

# from pymongo import MongoClient
# url = 'mongodb://fekad:qwe123@ds211613.mlab.com:11613/fekad_test'
# client = MongoClient(url)


if __name__ == '__main__':

    count = 0
    # A list to hold our things to do via async
    async_list = []

    direcotry = Path('data/')
    with ABCD(url='http://localhost:5000/api') as db:

        file = direcotry / 'bcc_bulk_54_expanded_2_high.xyz'
        file = direcotry / 'GAP_1.xyz'

        for atoms in iread(file.as_posix(), index=slice(None)):
            # Hack to fix the representation of forces
            atoms.calc.results['forces'] = atoms.arrays['force']

            db.push(atoms)

            count += 1
            print(count)

        # for file in direcotry.glob('*.xyz'):
        #     print(file)
        #
        #
        #     for atoms in iread(file.as_posix()):
        #         print(atoms)
        #         db.push(atoms)
        #
        #     break
