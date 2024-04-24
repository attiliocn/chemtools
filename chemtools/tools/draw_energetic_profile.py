#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
from energydiagram import ED as EnergyDiagram

parser = ArgumentParser()
parser.add_argument('input')
parser.add_argument('--format', default='svg', choices=['svg','png','jpg'])
args = parser.parse_args()

filename = args.input.split('.',1)[0]

# Below are the same settings use by the default generator 
#fig = plt.figure()
#ax = fig.add_subplot(111, aspect="equal")

fig, ax = plt.subplots(nrows=1, ncols=1)
energy_diagram = EnergyDiagram()

with open(args.input) as input_file:
    parsed_links = []
    for idx, diagram_element in enumerate(input_file):
        diagram_element_settings = diagram_element.strip().split(',')
        label = str(diagram_element_settings[0])
        energy = float(diagram_element_settings[1])
        x_position = str(diagram_element_settings[2])
        try: 
            int(x_position)
        except:
            if x_position == 'None':
                x_position = None
            else:
                x_position = str(x_position)
        requested_links = [int(i) for i in diagram_element_settings[3].split()]

        energy_diagram.add_level(energy=round(energy,1), bottom_text=label, position=x_position, top_text='Energy', linestyle='solid')
        
        for link in requested_links:
            parsed_links.append((idx,link))

for link in parsed_links:
    energy_diagram.add_link(*link)

#energy_diagram.dimension = 10
#energy_diagram.space = 10
energy_diagram.plot(ylabel="Energy / $kcal$ $mol^{-1}$", show_IDs=False, ax=ax)

#plt.savefig(f"{filename}.{args.format}")
plt.show()

