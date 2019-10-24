from setuptools import setup, find_packages

setup(
    name='abcd',
    version='0.5',
    author='Adam Fekete, Gabor Csanyi',
    author_email='adam.fekete@kcl.ac.uk',
    description='This is an package witch help to store and share atomistic data.',
    keywords='ase, database, mongo, flask',
    url='https://libatoms.github.io/abcd/',  # project home page, if any
    project_urls={
        'Documentation': 'https://libatoms.github.io/abcd/',
        'Source Code': 'https://github.com/libatoms/abcd',
    },
    packages=find_packages(),
    install_requires=['ase', 'numpy', 'lark-parser', 'pymongo'],
    extras_require={
        'tests': ['mongomock'],
        'mongo': ['pymongo'],
        'http': ['requests'],
        'server-api': ['flask'],
        'server-app': ['flask', 'Flask-Nav', 'Flask-MongoEngine', 'gunicorn', 'flask-paginate'],
    },
    test_require=['mongomock'],
    entry_points={
        'console_scripts': ['abcd=abcd.frontends.commandline:main']
    },
)
