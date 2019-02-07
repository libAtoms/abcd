import sys
from pathlib import Path

sys.path.append('..')

from abcd import ABCD
from utils.ext_xyz import XYZReader

if __name__ == '__main__':

    url = 'mongodb://localhost:27017'
    abcd = ABCD(url)

    for file in Path('GB_alphaFe_001/tilt/').glob('*.xyz'):
        print(file)
        gb_params = {
            'name': 'alphaFe',
            'type': 'tilt',
            'params': file.name[:-4]

        }
        with abcd as db, XYZReader(file) as reader:
            db.push(reader.read_atoms(), extra_info={'GB_params': gb_params})
