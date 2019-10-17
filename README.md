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
