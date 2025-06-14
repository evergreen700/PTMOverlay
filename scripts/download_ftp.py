import subprocess
import yaml
import sys
import os

credentials = sys.argv[1]
location = sys.argv[2]
outDir = sys.argv[3]

with open(credentials, "r") as inFile:
    cred = yaml.load(inFile, Loader=yaml.CLoader)

os.makedirs(outDir, exist_ok=True)

secure_flags = []
if "username" in cred:
    secure_flags.append("--user="+cred["username"])
if "password" in cred:
    secure_flags.append("--password="+cred["password"])

subprocess.run(["wget","--mirror","--continue","--no-host-directories","-nd"]+secure_flags+[cred["address"]+"/"+location],cwd=outDir)
