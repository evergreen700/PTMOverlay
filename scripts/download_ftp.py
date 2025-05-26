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

subprocess.run(["wget","--mirror","--continue","--no-host-directories","-nd","--user="+cred["username"],"--password="+cred["password"], cred["address"]+"/"+location],cwd=outDir)
