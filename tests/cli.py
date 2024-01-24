import os
import subprocess
import unittest
import logging
from time import sleep


class CLI(unittest.TestCase):
    """
    Testing OpenSearch database CLI integration.
    """

    @classmethod
    def setUpClass(cls):
        """
        Set up OpenSearch database connection and login with CLI.
        """
        if os.getenv("GITHUB_ACTIONS") != "true":
            raise unittest.SkipTest("Only runs via GitHub Actions")
        cls.security_enabled = os.getenv("security_enabled") == "true"
        cls.port = int(os.environ["port"])
        cls.host = "localhost"

        logging.basicConfig(level=logging.INFO)

        url = f"opensearch://admin:admin@{cls.host}:{cls.port}"
        if not cls.security_enabled:
            url += " --disable_ssl"
        try:
            subprocess.run(f"abcd login {url}", shell=True, check=True)
        except subprocess.CalledProcessError:
            sleep(10)
            subprocess.run(f"abcd login {url}", shell=True, check=True)

    def test_summary(self):
        """
        Test summary output of uploaded data file.
        """
        class_path = os.path.normpath(os.path.abspath(__file__))
        data_file = os.path.dirname(class_path) + "/data/example.xyz"
        subprocess.run(f"abcd upload {data_file}", shell=True, check=True)
        summary = subprocess.run(
            "abcd summary", shell=True, check=True, capture_output=True, text=True
        )
        assert "Total number of configurations: 1" in summary.stdout
