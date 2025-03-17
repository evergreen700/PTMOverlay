import os
import zipfile
import glob


with zipfile.ZipFile('mass_spec.zip') as zip_ref:
    zip_ref.extractall()

dirs = glob.glob(os.path.join("mass_spec","*"))

for d in dirs:
    bn = os.path.basename(d)
    ptms = bn.split("_")
    if len(ptms) > 1:
        for p in ptms:
            os.symlink(d,os.path.join("mass_spec",p))
