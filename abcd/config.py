import os
import json
from pathlib import Path


# https://martin-thoma.com/configuration-files-in-python/
class Config(dict):

    def __init__(self, file=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.file = file or os.environ.get('ABCD_CONFIG') or Path('./.abcd')

    @classmethod
    def from_file(cls, file=Path('./.abcd')):
        if not file.is_file():
            return cls(file)

        with open(file) as f:
            data = json.load(f)
        return cls(file, data)

    def save_json(self):
        with open(self.file, 'w') as file:
            json.dump(self, file)

    def __getitem__(self, key):
        return self.get(key) or os.environ.get(key)
