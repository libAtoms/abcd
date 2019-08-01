.. abcd documentation master file, created by
   sphinx-quickstart on Tue Jan 22 19:02:11 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to abcd's documentation!
================================

Goal:

- easy installation (docker support)
- easy to use (python classes for notebooks and command line for scripts)
- secure (api instead of direct db access, https)
- flexible (server-client application, flask blueprints)
- scalable database (nosql like mongodb)

Quickstart/Basic usage
----------------------

``python setup.py install``

::

    with ABCD(url='http://localhost:5000/api') as db:
        results = db.search('formula=Fe3O1;elements=[Fe,*];n_atoms=10,pbc;metadata.collection=iron')
        local_db = [db.get_atoms(id) for id in results]

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install.rst
   usage.rst
   tutorials.rst
   development.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
