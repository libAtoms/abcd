from ase.io import iread
from pathlib import Path

from pymongo import MongoClient


if __name__ == '__main__':
    url = 'mongodb://fekad:qwe123@ds211613.mlab.com:11613/fekad_test'
    client = MongoClient(url)

    direcotry = Path('data/')

    for file in direcotry.glob('*.xyz'):
        print(file)

        for atoms in iread(file.as_posix()):
            print(atoms)

            pass

        break
