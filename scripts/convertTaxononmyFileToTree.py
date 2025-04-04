import sys
import json

def getTaxonomyTree(inPath):
    taxonomyData = {}
    taxonomyTree = {}
    with open(inPath, "r") as inf:
        header = inf.readline()
        for line in inf:
            lineItems = line.strip().split("\t")
            taxonomyData[lineItems[0]] = lineItems[1:]

            taxonomy = lineItems[3]
            taxonomyItems = taxonomy.split("|")
            prevItem = ""
            prevTree = ""
            for taxItem in taxonomyItems:
                if prevTree == "":
                    if taxItem not in taxonomyTree:
                        taxonomyTree[taxItem] = {}
                    prevTree = taxonomyTree
                    prevItem = taxItem
                else:
                    if taxItem not in prevTree[prevItem]:
                        prevTree[prevItem][taxItem] = {}
                        prevTree = prevTree[prevItem]
                        prevItem = taxItem
                    else:
                        prevTree = prevTree[prevItem]
                        prevItem = taxItem
    return taxonomyTree

def writeTaxonomyTree(outPath, taxonomyTree):
    with open(outPath, "w") as outf:
        outf.write(f"{json.dumps(taxonomyTree)}")

if __name__ == "__main__":
    inPath = sys.argv[1]
    outPath = sys.argv[2]

    taxonomyTree = getTaxonomyTree(inPath)
    writeTaxonomyTree(outPath, taxonomyTree)

    print("done.")
