import os
import sys
from collections import deque
import json

inFile = sys.argv[1]
outFile = sys.argv[2]
inPTM = sys.argv[3:]
sequences = []
with open(inFile,"r") as readIn:
    seq = deque([readIn.readline()[1:]])
    line = readIn.readline()
    while line:
        if line[0] == ">":
            sequences.append(seq)
            seq = deque([line[1:]])
        else:
            seq.append(line)
        line=readIn.readline()
sequences.append(seq)

ptm_types = []
ptms = dict()
starts = dict()
ends = dict()
for f in inPTM:
    mod = os.path.basename(f).split("_")[1]
    ptm_types.append(mod)
    with open(f,"r") as readIn:
        p = json.load(readIn)
        for k in p.keys():
            if k not in ptms:
                ptms[k] = dict()
                starts[k] = dict()
                ends[k] = dict()
            ptms[k][mod] = p[k]['mod_site']
            for i,j in p[k]['read_start'].items():
                ii = int(i)
                if ii not in starts[k]:
                    starts[k][ii] = 0
                starts[k][ii]+=j
            for i,j in p[k]['read_end'].items():
                ii = int(i)
                if ii not in ends[k]:
                    ends[k][ii] = 0
                ends[k][ii]+=j

#for k in ptms.keys():
#    ptms[k]['read_start'] = {int(i):j for i,j in ptms[k]['read_start'].items()}
#    ptms[k]['read_end'] = {int(i):j for i,j in ptms[k]['read_end'].items()}

#print(ptm_types)
#print(ptms)
print(starts)
#print(ends)

colors=["red","blue","green","purple","orange"]

with open(outFile, "w") as writer:
    #print header
    writer.write('''
    <html>
        <head>
                <script src="https://code.jquery.com/jquery-3.1.1.min.js">

                </script><style> pre {       margin-left: 10px;  font-family: "Courier New";     font-size: 10pt;}.spacer{       height: 20px;}.lineHeader{      font-family: "Courier New";     font-size: 10pt;}''')
    for i in range(len(ptm_types)):
        writer.write('.'+ptm_types[i]+'{ color: '+colors[i]+'; --hidden-color: '+colors[i]+'; font-weight: bold;}')
    writer.write('''.secondaryHighlight{ font-weight: bold;      text-decoration: underline;     }.title{    font-weight: bold;    font-size: 14pt;    }
                </style>
        </head><body>
                <h1 class="title">
    ''')
    writer.write(outFile+'\n')
    writer.write('</h1><table>\n')
    keyOrder = []
    for i in range(len(sequences)):
        seq = sequences[i].popleft().strip()
        writer.write('<tr><td>'+f"{i:>2}:"+'</td><td>'+seq+'</td></tr>\n')
        keyOrder.append(seq)
    writer.write('''</table><table>
                        <tr>
                                <td class="spacer"></td>
                        </tr>''')

    coverages = {k:0 for k in keyOrder}

    iters = 0

    while len(sequences[0]) != 0:
        for i in range(len(sequences)):
            writer.write('<tr><td><span class="lineHeader" title="'+keyOrder[i]+'">'+str(i)+'</span></td><td><pre>')
            #writer.write(f"{i:>2} ")
            seq = sequences[i].popleft().strip()
            idx = iters
            for a in seq:
                c = ""
                if idx in starts[keyOrder[i]]:
                    coverages[keyOrder[i]]+=starts[keyOrder[i]][idx]
                if idx in ends[keyOrder[i]]:
                    coverages[keyOrder[i]]-=ends[keyOrder[i]][idx]
                if a != "-":
                    for pt in ptm_types:
                        if idx in ptms[keyOrder[i]][pt]:
                            c=pt+' '
                if coverages[keyOrder[i]] > 0:
                    c+='secondaryHighlight'
                writer.write('<span class="'+c+'" title>'+a+'</span>')
                idx+=1
            writer.write('</pre></td></tr>')
            if i == len(sequences)-1:
                iters=idx
        if len(sequences[0]) != 0:
                writer.write('''<tr class="spacer"></tr>''')

    writer.write(''' </table><form>
                    <table>''')
    for pt in ptm_types:
        writer.write('''<tr><td><input type="checkbox" name="'''+pt+'''" checked="checked" /><td>'''+pt+'''</td></td></tr>''')
    writer.write('''</table><script>
                            $(document).ready(function () {    $('input').on('change', function (d) {        $('.' + $(this).attr('name')).css('color', function (d) {            return $(this).css('color') === "rgb(0, 0, 0)" ? 'var(--hidden-color)' : "rgb(0, 0, 0)";        });    });});
                    </script>
            </form>
    </body>''')
