import sys
import json

def getAccessionIDs(accessionIDsPath):
    accessionIDs = []
    with open(accessionIDsPath) as inf:
        for line in inf:
            accessionIDs.append(line.strip())
    return accessionIDs

def getTaxonomiesDict(taxonomyOBjsPath):
    taxonomiesDict = {}
    with open(taxonomyOBjsPath) as inf:
        for line in inf:
            taxonomyObj = json.loads(line.strip())
            accessionID = taxonomyObj['query']
            lineage = taxonomyObj['lineage']
            # get the taxID
            taxID = lineage[0]['taxid']
            # get the scientific name
            scientificName = lineage[0]['name']
            # get the GenbankCommonName if it exists, otherwise get the first CommonName
            commonName = ""
            for name, nameType in lineage[0]['names'].items():
                if nameType == "GenbankCommonName":
                    commonName = name.title()
                    break
                if commonName == "" and nameType == "CommonName":
                    commonName = name.title()
            # fix weird apostrophy
            commonName = commonName.replace("'S", "'s")
            
            # get the taxonomy from the rest of the lineage objects
            taxonomyRanks = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"]
            taxonomyDict = dict.fromkeys(taxonomyRanks, "")
            for lineageObj in lineage[1:]:
                if lineageObj['rank'] in taxonomyDict:
                    taxonomyDict[lineageObj['rank']] = lineageObj['name']
            # order the taxonomy and put bars between ranks
            taxonomy = []
            for rank in taxonomyRanks:
                rankValue = taxonomyDict[rank]
                if rankValue != "":
                    taxonomy.append(taxonomyDict[rank])
            taxonomy.append(scientificName)
            taxonomy = "|".join(taxonomy)
            # append the entry to the dictionary
            taxonomiesDict[accessionID] = [taxID, commonName, taxonomy]

    return taxonomiesDict

def writeTaxonomiesFile(taxonomyPath, accessionIDs, taxonomiesDict):
    with open(taxonomyPath, "w") as outf:
        outf.write(f"accessionID\ttaxonomy\tcommonName\ttaxonomy\n")
        for accessionID in accessionIDs:
            if accessionID in taxonomiesDict:
                taxonomyData = taxonomiesDict[accessionID]
                outf.write(f"{accessionID}\t{taxonomiesDict[accessionID][0]}\t{taxonomiesDict[accessionID][1]}\t{taxonomiesDict[accessionID][2]}\n")
            else:
                outf.write(f"{accessionID}\tfailed\tfailed\tfailed\n")


if __name__ == '__main__':
    accessionIDsPath = sys.argv[1]
    taxonomyOBjsPath = sys.argv[2]
    taxonomyPath = sys.argv[3]

    accessionIDs = sorted(getAccessionIDs(accessionIDsPath))
    taxonomiesDict = getTaxonomiesDict(taxonomyOBjsPath)
    writeTaxonomiesFile(taxonomyPath, accessionIDs, taxonomiesDict)

    print("done.")