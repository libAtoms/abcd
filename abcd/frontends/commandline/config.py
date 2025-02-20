import json
import logging
import os
from pathlib import Path

logger = logging.getLogger(__name__)


class Config(dict):
    # https://martin-thoma.com/configuration-files-in-python/
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def from_json(cls, filename):
        with open(str(filename)) as json_file:
            obj = json.loads(json_file.read())

        return cls(obj)

    @classmethod
    def load(cls):
        if (
            os.environ.get("ABCD_CONFIG")
            and Path(os.environ.get("ABCD_CONFIG")).is_file()
        ):
            file = Path(os.environ.get("ABCD_CONFIG"))
        elif (Path.home() / ".abcd").is_file():
            file = Path.home() / ".abcd"
        else:
            return cls()

        logger.info(f"Using config file: {file}")

        return cls.from_json(file)

    def save(self):
        file = (
            Path(os.environ.get("ABCD_CONFIG"))
            if os.environ.get("ABCD_CONFIG")
            else Path.home() / ".abcd"
        )

        logger.info(f"The saved config's file: {file}")

        with open(str(file), "w") as file:
            json.dump(self, file)

    def __repr__(self):
        return f"<{self.__class__.__name__} {dict.__repr__(self)}>"
