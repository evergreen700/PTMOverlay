import os
import sys
from collections import deque
import json
import yaml

inFile = sys.argv[1]
outFile = sys.argv[2]
#window = int(sys.argv[3])
inPTM = sys.argv[3:]

#these are the colors assigned to the PTMs in the order they are presented in.
colors=["#fc1c03","#0324fc","#a903fc","#03fc0f","#fc03e7","#b59a4a"]

os.makedirs(os.path.dirname(outFile), exist_ok=True)

with open("config.yaml","r") as confFile:
    config = yaml.load(confFile, Loader=yaml.CLoader)

ranges = {"starts":set(),"ends":set()}
if "select_intervals" in config:
    for s,e in config["select_intervals"]:
        ranges["starts"].add(s)
        ranges["ends"].add(e)
else:
    ranges["starts"].add(0)

window = config.get("row_length",80)

sequences = []
seqIDs = []
with open(inFile,"r") as readIn:
    seqIDs.append(readIn.readline().strip()[1:])
    seq = []
    line = readIn.readline().strip()
    while line:
        if line[0] == ">":
            sequences.append("".join(seq))
            seqIDs.append(line[1:])
            seq = []
        else:
            seq.append(line)
        line=readIn.readline().strip()
sequences.append("".join(seq))

ptm_types = []
ptms = dict()
starts = dict()
ends = dict()

#read in PTM files
for f in inPTM:
    mod = os.path.basename(f).split("_")[1]
    ptm_types.append(mod)
    with open(f,"r") as readIn:

        #merge ptm sites, read starts, and read ends into dictionaries
        p = json.load(readIn)
        for k in p.keys():
            if k not in ptms:
                ptms[k] = dict()
                starts[k] = dict()
                ends[k] = dict()

            #store ptm sites
            ptms[k][mod] = p[k]['mod_site']

            #store read starts
            for i,j in p[k]['read_start'].items():
                ii = int(i)
                if ii not in starts[k]:
                    starts[k][ii] = 0
                starts[k][ii]+=j

            #store read ends
            for i,j in p[k]['read_end'].items():
                ii = int(i)
                if ii not in ends[k]:
                    ends[k][ii] = 0
                ends[k][ii]+=j


with open(outFile, "w") as writer:

    #write header
    writer.write('''
    <html>
        <head>
                <script src="https://code.jquery.com/jquery-3.1.1.min.js">

                </script><style> pre {       margin-left: 10px;  font-family: "Courier New";     font-size: 10pt;}.spacer{       height: 20px;}.lineHeader{      font-family: "Courier New";     font-size: 10pt;}''')

    #This adds a class for each ptm type
    for i in range(len(ptm_types)):
        writer.write('.'+ptm_types[i]+'{ color: '+colors[i]+'; --hidden-color: '+colors[i]+'; font-weight: bold;}')

    #This completes the header and starts the body
    writer.write('''.secondaryHighlight{ font-weight: bold;      text-decoration: underline;     }.title{    font-weight: bold;    font-size: 14pt;    }
                </style>
        </head><body>
                <h1 class="title">
    ''')
    writer.write(os.path.basename(outFile))
    writer.write('</h1><table>\n')
    for i in range(len(sequences)):
        seq = seqIDs[i]
        writer.write('<tr><td>'+str(i+1)+':</td><td>'+seq+'</td></tr>\n')
    writer.write('''</table><table>
                        <tr>
                                <td class="spacer"></td>
                        </tr>''')

    coverages = {k:0 for k in seqIDs}
    writing = {k:0 for k in seqIDs}
    iters = 0
    while iters < len(sequences[0]):
        for i in range(len(sequences)):
            writer.write('<tr><td><span class="lineHeader" title="'+seqIDs[i]+'">'+str(i+1)+'</span></td><td><pre>')
            seq = sequences[i]
            idx = iters
            writePos = 0
            while writePos < window and idx < len(sequences[0]):
                a = seq[idx]
                c = ""
                if idx in starts[seqIDs[i]]:
                    coverages[seqIDs[i]]+=starts[seqIDs[i]][idx]
                if idx in ends[seqIDs[i]]:
                    coverages[seqIDs[i]]-=ends[seqIDs[i]][idx]
                if idx in ranges["starts"]:
                    writing[seqIDs[i]]+=1
                if idx in ranges["ends"]:
                    writing[seqIDs[i]]-=1
                    if writing[seqIDs[i]] == 0 and idx != max(ranges["ends"]):
                        writer.write('<span title> … </span>')
                        writePos+=3
                if a != "-":
                    for pt in ptm_types:
                        if idx in ptms[seqIDs[i]][pt]:
                            c=pt+' '
                if coverages[seqIDs[i]] > 0:
                    c+='secondaryHighlight'
                if writing[seqIDs[i]] > 0:
                    writer.write('<span class="'+c+'" title>'+a+'</span>')
                    writePos += 1
                idx+=1
            writer.write('</pre></td></tr>')
        iters=idx
        if iters < len(sequences[0]):
            writer.write('''<tr class="spacer"></tr>''')

    writer.write(''' </table><form>
                    <table>''')
    for pt in ptm_types:
        writer.write('''<tr><td><input type="checkbox" name="'''+pt+'''" checked="checked" /><td class="'''+pt+'''">'''+pt+'''</td></td></tr>''')
    writer.write('''</table><script>
                            $(document).ready(function () {    $('input').on('change', function (d) {        $('.' + $(this).attr('name')).css('color', function (d) {            return $(this).css('color') === "rgb(0, 0, 0)" ? 'var(--hidden-color)' : "rgb(0, 0, 0)";        });    });});
                    </script>
            </form>
    </body>''')
