import svg_stack as ss
import sys

doc = ss.Document()

input_tree = sys.argv[1]
input_graph = sys.argv[2]
output = sys.argv[3]

layout1 = ss.HBoxLayout()
layout1.addSVG(input_tree,alignment=ss.AlignCenter|ss.AlignHCenter)
layout1.addSVG(input_graph,alignment=ss.AlignTop)

layout2 = ss.VBoxLayout()

layout1.addLayout(layout2)

doc.setLayout(layout1)

doc.save(output)