# Copyright (c) 2025.
# Authors: Ádám Fekete, Elliott Kasoar, Tamás K. Stenczel
# This program is distributed under the MIT License, see LICENSE.md.

# Copyright (c) 2025.
# Authors: Ádám Fekete, Elliott Kasoar, Tamás K. Stenczel
# This program is distributed under the MIT License, see LICENSE.md.

import sys
from pathlib import Path
import sys

sys.path.append("..")

from utils.ext_xyz import XYZReader

from abcd import ABCD

if __name__ == "__main__":
    url = "mongodb://localhost:27017"
    abcd = ABCD(url)

    for file in Path("GB_alphaFe_001/tilt/").glob("*.xyz"):
        print(file)
        gb_params = {"name": "alphaFe", "type": "tilt", "params": file.name[:-4]}
        with abcd as db, XYZReader(file) as reader:
            db.push(reader.read_atoms(), extra_info={"GB_params": gb_params})
