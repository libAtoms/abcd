from setuptools import setup, find_packages

setup(
    name='abcd',
    version='0.1',
    packages=['abcd'],
    # packages=find_packages(),
    install_requires=['ase'],
    entry_points={
        'console_scripts': ['abcd-server = abcd_server.cmd:main']},

    # metadata to display on PyPI
    # author="Me",
    # author_email="me@example.com",
    # description="This is an Example Package",
    # license="PSF",
    # keywords="hello world example examples",
    # url="http://example.com/HelloWorld/",  # project home page, if any
    # project_urls={
    #     "Bug Tracker": "https://bugs.example.com/HelloWorld/",
    #     "Documentation": "https://docs.example.com/HelloWorld/",
    #     "Source Code": "https://code.example.com/HelloWorld/",
    # }

)
