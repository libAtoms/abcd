# ABCD

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

To connect to a mongodb running on the mac from inside the docker that runs abcd, use 

```
abcd login mongodb://docker.for.mac.localhost
```
