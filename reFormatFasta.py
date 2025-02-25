import os
import sys
from collections import deque
import json

inFile = sys.argv[1]
inPTM = sys.argv[2]
outFile = sys.argv[3]
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

with open(inPTM,"r") as readIn:
    ptms = json.load(readIn)

for k in ptms.keys():
    ptms[k]['read_start'] = {int(i):j for i,j in ptms[k]['read_start'].items()}
    ptms[k]['read_end'] = {int(i):j for i,j in ptms[k]['read_end'].items()}

print(ptms)

with open(outFile, "w") as writer:
    #print header
    writer.write('''
    <html>
        <head>
                <script src="https://code.jquery.com/jquery-3.1.1.min.js">

                </script><style>
                        pre {       margin-left: 10px;  font-family: "Courier New";     font-size: 10pt;}.spacer{       height: 20px;}.lineHeader{      font-family: "Courier New";     font-size: 10pt;}.primaryHighlight{     color: rgb(255, 0, 0);     font-weight: bold;}.secondaryHighlight{ font-weight: bold;      text-decoration: underline;     }.title{    font-weight: bold;    font-size: 14pt;    }
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

    while len(sequences[0]) != 0:
        iters = 0
        for i in range(len(sequences)):
            writer.write('<tr><td><span class="lineHeader" title="'+keyOrder[i]+'">'+str(i)+'</span></td><td><pre>')
            #writer.write(f"{i:>2} ")
            seq = sequences[i].popleft().strip()
            idx = iters
            for a in seq:
                c = ""
                if idx in ptms[keyOrder[i]]['read_start']:
                    coverages[keyOrder[i]]+=ptms[keyOrder[i]]['read_start'][idx]
                if idx in ptms[keyOrder[i]]['read_end']:
                    coverages[keyOrder[i]]-=ptms[keyOrder[i]]['read_end'][idx]
                if a != "-" and idx in ptms[keyOrder[i]]['mod_site']:
                    c='primaryHighlight Phospho '
                if coverages[keyOrder[i]] > 0:
                    c+='secondaryHighlight'
                writer.write('<span class="'+c+'" title>'+a+'</span>')
                idx+=1
            writer.write('</pre></td></tr>')
        if len(sequences[0]) != 0:
                writer.write('''<tr class="spacer"></tr>''')

        iters=idx
    writer.write(''' </table><form>
                    <table>
                            <tr>
                                    <td><input type="checkbox" name="Phospho" checked="checked" /><td>Phospho</td></td>
                            </tr>
                    </table><script>
                            $(document).ready(function () {    $('input').on('change', function (d) {        $('.' + $(this).attr('name')).css('color', function (d) {            return $(this).css('color') === "rgb(255, 0, 0)" ? "rgb(0, 0, 0)" : "rgb(255, 0, 0)";        });    });});
                    </script>
            </form>
    </body>''')
