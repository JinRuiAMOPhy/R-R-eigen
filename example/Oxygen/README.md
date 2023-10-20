# expample for O autoionization
# generate targets files and inputs
```
python script.py  
```

This generate a file n2, which contains the target files, and a directory Reigen-inner

```
cd n2/Reigen-inner
```
run inner region calculations, stgb calculations and then
run 
```
python runstgf.py 
```
for interested continuum energies and partialwaves.

run
```
python smooth.py
```
to get a regularized S-matrices data file. 
and 
```
python dfde.py
```
to generate cross sections.

