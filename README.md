# ABCD

[![Doc](https://img.shields.io/badge/docs-master-green.svg)](https://libatoms.github.io/abcd/)
[![Build Status](https://travis-ci.org/libAtoms/abcd.svg?branch=master)](https://travis-ci.org/libAtoms/abcd)

Database storage and discovery of atomistic data.

Take a look at the `examples.md` file for.. examples!

Main features:

- Configurations that consist of atom positions, elements, forces, and various metadata are stored as a dictionary by a MongoDB backend.
- There is no predefined schema, any combination of keys are allowed for all configurations.
- Two modes: "discovery" and "download". Both use filter-type queries, but in "discovery" mode, summary statistics of the configurations that pass the filter are reported. In "download" mode, the matching configurations are downloaded and exported to a file.
- The "discovery" mode can be used to learn what keys exist in the set of configurations that have passed the current query filter. The user can use this to refine the query.
- Complex queries on dictionary key-value pairs are allowed, and their logical combinations.

## Installation

### General Setup

creating tables and views

```sh
$ pip install git+https://github.com/libAtoms/abcd.git
```

Example Docker installation on Ubuntu:

```sh
sudo apt-get update
sudo apt upgrade
sudo apt install docker.io
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker # or exit and log in
```

Docker can be tested by running:

```sh
docker run hello-world
```

Example Python setup on Ubuntu (pip must be updated for poetry to be used successfully):

```sh
sudo apt install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt install python3.10
sudo apt-get install python3.10-distutils
sudo apt install python3-virtualenv
virtualenv -p /usr/bin/python3.10 venv_10
source venv_10/bin/activate
curl -sS https://bootstrap.pypa.io/get-pip.py | python3.10
```

Building and installing ABCD dependencies via poetry:

```sh
git clone https://github.com/libAtoms/abcd.git
curl -sSL https://install.python-poetry.org | python3 -
export PATH="/home/ubuntu/.local/bin:$PATH"
cd abcd
poetry install
poetry build
```

### MongoDB

If you have an already running MongoDB server, or install your own, then you are ready to go. Alternatively,

```sh
docker run -d --rm --name abcd-mongodb -v <path-on-your-machine-to-store-database>:/data/db -p 27017:27017 mongo
```

will download and install a docker and run a database in it.

To connect to a mongodb that is already running, use

```sh
abcd login mongodb://localhost
```

If you are running `abcd` inside a docker, and want to connect to a mongodb outside that docker use something like this (example is for Mac OS):

```sh
abcd login mongodb://docker.for.mac.localhost
```

The above login command will place create an `~/.abcd` file with the following contents:

```sh
{"url": "mongodb://localhost"}
```

### OpenSearch
If you have an already running OpenSearch server, or install your own, then you are ready to go. Alternatively,

```sh
sudo swapoff -a # optional
sudo sysctl -w vm.swappiness=1 # optional
sudo sysctl -w fs.file-max=262144 # optional
sudo sysctl -w vm.max_map_count=262144
docker run -d --rm --name abcd-opensearch -v <path-on-your-machine-to-store-database>:/data/db -p 9200:9200  --env discovery.type=single-node -it opensearchproject/opensearch:latest
```

will download and install an OpenSearch image and run it. The connection can be tested with:

```sh
curl -vvv -s --insecure -u admin:admin --fail https://localhost:9200
```

To connect to an OpenSearch database that is already running, use

```sh
abcd login opensearch://username:password@localhost
```

## Remote access

You can set up an `abcd` user on your machine where the database is running, and then access it remotely for discovering data. Make sure you have the `~/.abcd` file created for this user, then put this in the `.ssh/authorized_keys` file (substituting your public key for the last part):

```sh
command="/path/to/abcd --remote  ${SSH_ORIGINAL_COMMAND}",no-port-forwarding,no-X11-forwarding,no-agent-forwarding,no-pty ssh-rsa <public-key> your@email
```

Then you'll be able to access the database remotely using, e.g.

```sh
ssh abcd@your.machine summary
```

## GUI through a browser + visualisation

The database has a simple GUI, coupled with a visualiser. Data for now needs to be uploaded on the command line, but query can be done through the browsers. Instructions below (they include running `abcd` from a docker too, but of course you can run it outside the docker as well. )


## Usage in docker
Currently a manual uploaded image is available, that was built on 7/2/2020 by Tamas K. Stenczel.
To access it:
1. pull the image
    ```sh
    docker pull stenczelt/projection-abcd:latest
    ```

2. create a docker network, which enables the containers to communicate with each other and the outside world as well 
    ```sh
    docker network create --driver bridge abcd-network
    ```

3. run the mongo (ABCD) and the visualiser as well
    ```sh
    docker run -d --rm --name abcd-mongodb-net -v <path-on-your-machine-to-store-database>:/data/db -p 27017:27017 --network abcd-network mongo

    docker run -it --rm --name visualiser-dev -p 9999:9999 --network abcd-network stenczelt/projection-abcd
    ```
    NB: You need a a directory where the database files are kept locally and you need to connect this to the mongo 
    container. More info about this can be found in the original ABCD repo

This will start the visualiser with ABCD integration! Have fun!

After usage, for cleanup:

```sh
docker stop visualiser-dev abcd-mongodb-net         # stop the containers
docker rm visualiser-dev abcd-mongodb-net           # remove them if --rm did not
docker network rm abcd-network                      # remove the docker network
```

## Testing

Unit tests are automatically run on push and creation of pull requests. Unit testing using mock databases can also be run in the command line using:

```sh
python -m unittest tests
```
