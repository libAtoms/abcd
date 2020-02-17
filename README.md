# ABCD

[![Doc](https://img.shields.io/badge/docs-master-green.svg)](https://libatoms.github.io/abcd/)
[![Build Status](https://travis-ci.org/libAtoms/abcd.svg?branch=master)](https://travis-ci.org/libAtoms/abcd)

storage and discovery of atomistic data

## Installation

creating tables and views
```
$ pip install git+https://github.com/libAtoms/abcd.git
```

## Setup

If you have an already running mongo server, or install your own, they you are ready to go. Alternatively, 

```
docker run -d --rm --name abcd-mongodb -v <path-on-your-machine-to-store-database>:/data/db -p 27017:27017 mongo
```

will download and install a docker and run a database in it. 

To connect to a mongodb that is already running, use
```
abcd login mongodb://localhost
```

If you are running `abcd` inside a docker, and want to connect to a mongodb outside that docker use something like this (example is for Mac OS):

```
abcd login mongodb://docker.for.mac.localhost
```

The above login command will place create an `~/.abcd` file with the following contents:

```
{"url": "mongodb://localhost"}
```

# Remote access

You can set up an `abcd` user on your machine where the database is running, and then access it remotely for discovering data. Make sure you have the `~/.abcd` file created for this user, then put this in the `.ssh/authorized_keys` file (substituting your public key for the last part):
```
command="/path/to/abcd --remote  ${SSH_ORIGINAL_COMMAND}",no-port-forwarding,no-X11-forwarding,no-agent-forwarding,no-pty ssh-rsa <public-key> your@email
```

Then you'll be able to access the database remotely using, e.g. 
```
ssh abcd@your.machine summary
```

# GUI through a browser + visualisation

The database has a simple GUI, coupled with a visualiser. Data for now needs to be uploaded on the command line, but query can be done through the browsers. Instructions below (they include running `abcd` from a docker too, but of course you can run it outside the docker as well. )


#### Usage in docker
Currently a manual uploaded image is available, that was built on 7/2/2020 by Tamas K. Stenczel.
To access it:
1. pull the image

    `docker pull stenczelt/projection-abcd:latest`

2. create a docker network, which enables the containers to communicate with each other and the outside world as well 

    `docker network create --driver bridge abcd-network`

3. run the mongo (ABCD) and the visualiser as well

    `docker run -d --rm --name abcd-mongodb-net -v <path-on-your-machine-to-store-database>:/data/db -p 27017:27017 --network abcd-network mongo`

    `docker run -it --rm --name visualiser-dev -p 9999:9999 --network abcd-network stenczelt/projection-abcd`

    NB: You need a a directory where the database files are kept locally and you need to connect this to the mongo 
    container. More info about this can be found in the original ABCD repo

This will start the visualiser with ABCD integration! Have fun!
After usage, for cleanup:

    `docker stop visualiser-dev abcd-mongodb-net         # stop the containers`
    
    `docker rm visualiser-dev abcd-mongodb-net           # remove them if --rm did not`
    
    `docker network rm abcd-network                      # remove the docker network`

