import sys
import json

partials = sys.argv[1:-2]
kid = sys.argv[-2]
outFile = sys.argv[-1]

complete = dict()

for inFile in partials:
    with open(inFile, "r") as inBuffer:
        inPTM = json.load(inBuffer)
        complete.update(inPTM.get(kid,{}))

with open(outFile, "w") as outBuffer:
    json.dump(complete, outBuffer, indent=4)
