from svgutils.compose import Figure, SVG
import sys

input = sys.argv[1]
output = sys.argv[2]

originalSVG = SVG(input)
originalSVG.scale(1.5)
figure = Figure(float(originalSVG.height) * 1.5, float(originalSVG.width) * 1.5, originalSVG)
figure.save(output)