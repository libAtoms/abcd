# ABCD

storage and discovery of atomistic data

## usage (to be implemented)

creating tables and views
```
$ psql -d $DB -f sql/create.sql
$ psql -d $DB -f sql/views.sql 
```

importing data
```
$ abcd-import [-db $DB] GAP1.xyz GAP2.xyz 
```

querying data
```
$ abcd [-db $DB] <query> -export stuff.xyz
$ abcd [-db $DB] <query> -hist var
$ abcd [-db $DB] <query> -summary var,var,var
$ abcd [-db $DB] <query> -count
```

using the environment variable `ABCD_DB` if `-db` not specified.

query format
```
query := var <op> value
op := -eq | -gt | -lt
```

## examples

Import Tungsten data from libatoms.org
```
$ abcd-import data/Tungsten/GAP_?.xyz
```

Find what keys are available in the dataset
```
$ abcd --keys
frame keys
----------
     34669 config_type
     34669 total_energy
     12000 virial
       549 virial_not
atom keys
----------
```
which shows that there are 4 keys present. `config_type` and `total_energy` are present for all frames whereas `virial` and `virial_not` are present only for a subset of frames.

Summarise the distribution of a variable
```
$ abcd --summarise total_energy
total_energy : min -144.300776, avg -1509.120000, max -10.822000
-1509.120000 -1414.490000    480
 -589.467000  -510.259000   2966 ■■■
 -510.217000  -507.357000    134
 -134.340000   -10.826900  31083 ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
  -10.822000   -10.822000      6
```
which shows some basic statistics and something which might charitably called a histogram.
