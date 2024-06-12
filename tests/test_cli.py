import logging
import os
from pathlib import Path
import subprocess
from time import sleep

import pytest


DATA_PATH = Path(__file__).parent / "data"

NOT_GTHUB_ACTIONS = True
if os.getenv("GITHUB_ACTIONS") == "true":
    NOT_GTHUB_ACTIONS = False


@pytest.mark.skipif(NOT_GTHUB_ACTIONS, reason="Not running via GitHub Actions")
class TestCli:
    """Testing OpenSearch database CLI integration."""

    @pytest.fixture(autouse=True)
    def abcd(self):
        """Set up OpenSearch database connection and login with CLI."""
        security_enabled = os.getenv("security_enabled") == "true"
        port = int(os.environ["port"])
        host = "localhost"
        if os.environ["opensearch-version"] == "latest":
            credential = "admin:myStrongPassword123!"
        else:
            credential = "admin:admin"

        logging.basicConfig(level=logging.INFO)

        url = f"opensearch://{credential}@{host}:{port}"
        if not security_enabled:
            url += " --disable_ssl"
        try:
            subprocess.run(f"abcd login {url}", shell=True, check=True)
        except subprocess.CalledProcessError:
            sleep(10)
            subprocess.run(f"abcd login {url}", shell=True, check=True)

    def test_summary(self, abcd):
        """
        Test summary output of uploaded data file.
        """
        data_file = DATA_PATH / "example.xyz"

        subprocess.run(
            f"abcd upload {data_file} -i -e 'test_data'", shell=True, check=True
        )
        subprocess.run(f"abcd refresh", shell=True, check=True)

        summary = subprocess.run(
            "abcd summary", shell=True, check=True, capture_output=True, text=True
        )
        assert "Total number of configurations" in summary.stdout
        subprocess.run(f"abcd delete -q 'test_data' -y", shell=True)

    def test_query(self, abcd):
        """
        Test lucene-style query.
        """
        data_file_1 = DATA_PATH / "example.xyz"
        data_file_2 = DATA_PATH / "example_2.xyz"

        subprocess.run(
            f"abcd upload {data_file_1} -i -e 'test_data'", shell=True, check=True
        )
        subprocess.run(
            f"abcd upload {data_file_2} -i -e 'test_data'", shell=True, check=True
        )
        subprocess.run(f"abcd refresh", shell=True, check=True)

        summary = subprocess.run(
            "abcd show -p n_atoms -q 'n_atoms : 2'",
            shell=True,
            check=True,
            capture_output=True,
            text=True,
        )
        assert "2" in summary.stdout and "3" not in summary.stdout
        summary = subprocess.run(
            "abcd show -p n_atoms -q 'n_atoms : 3'",
            shell=True,
            check=True,
            capture_output=True,
            text=True,
        )
        assert "3" in summary.stdout and "2" not in summary.stdout
        subprocess.run(f"abcd delete -q 'test_data' -y", shell=True)

    def test_range_query(self, abcd):
        """
        Test lucene-style ranged query.
        """
        data_file_1 = DATA_PATH / "example.xyz"
        data_file_2 = DATA_PATH / "example_2.xyz"

        subprocess.run(
            f"abcd upload {data_file_1} -i -e 'test_data'", shell=True, check=True
        )
        subprocess.run(
            f"abcd upload {data_file_2} -i -e 'test_data'", shell=True, check=True
        )
        subprocess.run(f"abcd refresh", shell=True, check=True)

        summary = subprocess.run(
            "abcd summary -p energy -q 'energy:[-100 TO -99]'",
            shell=True,
            check=True,
            capture_output=True,
            text=True,
        )
        assert "Total number of configurations: 1" in summary.stdout

        summary = subprocess.run(
            "abcd summary -p energy -q 'energy:[-102 TO -99]'",
            shell=True,
            check=True,
            capture_output=True,
            text=True,
        )
        assert "Total number of configurations: 2" in summary.stdout
        subprocess.run(f"abcd delete -q 'test_data' -y", shell=True)
