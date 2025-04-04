import json
import re
import sys

def jsonToNewick(data, parentName=''):
    newickStr = ''
    if isinstance(data, dict):
        for childName, childData in data.items():
            if newickStr:
                newickStr += ','
            newickStr += jsonToNewick(childData, childName)
    
    if parentName:
        return f'({newickStr}){parentName}'
    else:
        return newickStr

def jsonStrToNewick(jsonStr):
    newickStr = jsonStr.replace(": {}", "")
    newickStr = re.sub(r"\"[\w\(\).]+( [\w\(\).]+)*\": ", "", newickStr)
    newickStr = re.sub(", ", ",", newickStr)
    # newickStr = newickStr.replace(",{}", "").replace(",}", "}")
    newickStr = newickStr.replace(r"\(", "").replace(r"\)", "")
    # newickStr = newickStr.replace("\"", "")
    newickStr = newickStr.replace("(", "").replace(")", "").replace(",}", "}")
    # replace curly braces with parentheses
    newickStr = newickStr.replace("{", "(").replace("}", ")")
    # remove all colons
    newickStr = newickStr.replace(":", " ")
    # remove all quotes
    newickStr = newickStr.replace("\"", "")
    # add a semi colon at the end if it isn't there
    if newickStr[-1] != ";":
        newickStr = newickStr + ";"
    return newickStr

def readFile(filePath):
    with open(filePath, 'r') as file:
        data = file.readline()
    return data

def writeNewickFile(filePath, newickStr):
    with open(filePath, 'w') as file:
        file.write(newickStr)

    
if __name__ == '__main__':
    jsonFilePath = sys.argv[1]
    newickFilePath = sys.argv[2]

    jsonStr = readFile(jsonFilePath)
    newickFormat = jsonStrToNewick(jsonStr)
    writeNewickFile(newickFilePath, newickFormat)

    print("done.")