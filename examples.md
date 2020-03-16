# Examples

```
abcd login mongodb://localhost
abcd login mongodb://mongoadmin:secret@localhost:27017/abcd_new
```
```
abcd upload -e cas tungsten_database.xyz
abcd upload -e cas VelocityVerlet_flexible_step0.5.xyz
abcd upload -e cas silicon_database_gp_iter6_sparse9k.xml.xyz
abcd upload -e cas Ti_N54_database.xyz
```
```
abcd summary
abcd summary -p formula
abcd summary -p formula --all
```
```
abcd summary -q 'formula="H250O125"' -p formula
abcd summary -q 'formula="H250O125"' -q cas -p formula
```
Note how all strings need to be quoted, and those quotes protected from the shell. 

```
abcd rename-key --help
abcd rename-key -q 'formula="H250O125"' cas eszter
```
```
abcd summary
abcd summary -q eszter -p formula
abcd summary -q eszter -p energy
```
```
abcd summary
abcd summary -p formula
abcd summary -q 'formula~"W"' -p formula
abcd summary -q 'formula~"W"' -p energy
```
the `~` operator does RegExp matching
```
abcd summary -q 'formula~"W128"'
abcd download -q 'formula~"W128"' -f xyz data.xyz
head data.xyz
```
```
abcd exec -q 'formula~"W128"' 'print(item.info.cell)'
abcd exec -q 'formula~"W128"' 'print(item.info.cell)' --yes
abcd exec -q 'formula~"W128"' 'at=item.to_atoms(); print(at.cell)' --yes
```
```
abcd exec  'atoms["energy_per_atom"] = atoms["energy"]/atoms["n_atoms"] ; atoms.save()' --yes
abcd summary -p energy_per_atom
```
