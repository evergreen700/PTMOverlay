import os
import requests
import zipfile
import glob

os.makedirs('mass_spec',exist_ok=True)
response = requests.get('https://byu.box.com/shared/static/seu2e75iw7ihnxqnsp8gi29xmo9pqynr.zip', stream=True)
with open('mass_spec.zip', 'wb') as f:
    for chunk in response.iter_content(chunk_size=8192):
        f.write(chunk)

with zipfile.ZipFile('mass_spec.zip') as zip_ref:
    zip_ref.extractall()

dirs = glob.glob(os.path.join("mass_spec","*"))

for d in dirs:
    bn = os.path.basename(d)
    os.rename(d,os.path.join("mass_spec",bn.replace("_","")))
