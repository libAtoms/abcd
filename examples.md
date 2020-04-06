# Examples

```
abcd login mongodb://localhost
abcd login mongodb://mongoadmin:secret@localhost:27017/abcd_new
```
```
abcd upload -e cas tungsten_database.xyz
abcd upload -e cas VelocityVerlet_flexible_step0.5.xyz
abcd upload -e cas silicon_database_gp_iter6_sparse9k.xml.xyz
# get around the silly calculator/atoms object storage nightmare of ASE:
# the -i option tells ABCD to ignore the fake ASE calculator that gets created when xyz files are imported
abcd upload -i -e cas Ti_N54_database.xyz 
```
The most often used command is `summary`, it tells you how many configurations match a query. Together with the `-p` option, it gives histograms of those properties.
```
abcd summary
abcd summary -p formula
abcd summary -p formula --all
```
```
abcd summary -q 'formula="H250O125"' -p formula
abcd summary -q 'formula="H250O125"' -q cas -p formula
# regular expressions! the `~` operator does RegExp matching
abcd summary -q 'formula~"H.*O.*"' -q cas -p formula
```
Note how all strings need to be quoted, and those quotes protected from the shell. 

```
abcd rename-key --help
abcd rename-key -q 'formula="H250O125"' cas eszter
```
```
abcd summary
abcd summary -q eszter -p formula
# -p on a real valued property gives the distribution of that property
abcd summary -q eszter -p energy
```
```
abcd summary
abcd summary -p formula
abcd summary -q 'formula~"W"' -p formula
abcd summary -q 'formula~"W"' -p energy
```

downloading configurations:
```
abcd download -q 'formula~"W128"' -f xyz data.xyz
```

executing commands on all configurations that match a query:
```
abcd exec -q 'formula~"W128"' 'print(item.info.cell)'
abcd exec -q 'formula~"W128"' 'print(item.info.cell)' --yes
abcd exec -q 'formula~"W128"' 'at=item.to_atoms(); print(at.cell)' --yes
```
```
abcd exec  'atoms["energy_per_atom"] = atoms["energy"]/atoms["n_atoms"] ; atoms.save()' --yes
abcd summary -p energy_per_atom
```
