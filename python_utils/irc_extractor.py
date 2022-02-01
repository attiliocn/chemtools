#!/usr/bin/env python

import argparse
from email.policy import default
import re
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(dest='subparser')

parser_single = subparsers.add_parser('single', help='For single output IRC calculation (forward/reverse in the same output')
parser_single.add_argument('output', type=str, help='Gaussian Output File')
parser_single.add_argument('product', choices=['rv','fw'], help='Indicate whether the product is in the forward or reverse direction')
parser_single.add_argument('--downhill', action='store_true', help='Special option for downhill calculations')

parser_double = subparsers.add_parser('double', help='For IRC in which forward and reverse were calculated separately')
parser_double.add_argument('reverse', type=str, help='Gaussian Output File (reverse)')
parser_double.add_argument('forward', type=str, help='Gaussian Output File (forward)')
parser_double.add_argument('product', choices=['rv','fw'], help='Indicate whether the product is in the forward or reverse direction')

parser_single.add_argument('--title', type=str, default=None, help='Filename and Plot Title')
parser_double.add_argument('--title', type=str, default=None, help='Filename and Plot Title')

parser_single.add_argument('--include-title', action='store_true', help='Include plot title in the figure')
parser_double.add_argument('--include-title', action='store_true', help='Include plot title in the figure')

args = parser.parse_args()

def extract_irc_from_output(output_file):
    '''
    This function takes a Gaussian output file and extract the IRC path from it. 
    Returns a 2-column IRC (x= IRC displacement, y= Total energy)
    '''  
    irc_coordinates = [0.0]
    irc_number_of_steps = list()
    irc_energy_all_points = list()

    with open(output_file) as f:
        irc_data = f.readlines()
        for line in irc_data:
            
            intrinsic_coord = re.search('NET REACTION COORDINATE UP TO THIS POINT', line)
            number_of_steps = re.search('# OF STEPS', line)
            scf_energy = re.search("(SCF Done).*(-[0-9]+.[0-9]+)", line)

            if intrinsic_coord:
                irc_coordinates.append(float(line.split()[-1]))
            elif number_of_steps:
                irc_number_of_steps.append(int(line.split()[-1]))
            elif scf_energy:
                irc_energy_all_points.append(float(scf_energy.group(2)))

    irc_energies = list()
    irc_energies.append(irc_energy_all_points[0])

    steps_cumulative = 0 
    for steps in irc_number_of_steps:
        steps_cumulative += steps
        irc_energies.append(irc_energy_all_points[steps_cumulative])

    irc_plot_data = np.column_stack([irc_coordinates, irc_energies])
    return irc_plot_data

def split_singlecalc_irc(irc_calculation):
    for i in range(len(irc_calculation)):
        if irc_calculation[i,0] - irc_calculation[i+1,0] > 0:
            fw_limit = i
            break
    irc_rv = irc_calculation[fw_limit+1:, :]
    irc_fw = irc_calculation[0:fw_limit+1, :]
    return irc_rv, irc_fw

def join_partial_irc(reactant_irc, product_irc):
    '''
    This function takes the processed IRC calculation (as 2-column np_array) for 
    both reactant and product and returns the full IRC calculation
    '''
    # reactant
    processed_reactant_irc = reactant_irc
    processed_reactant_irc[:,0] *= -1
    processed_reactant_irc.sort(axis=0)
    # product
    processed_product_irc = product_irc
    # full irc
    full_irc = np.row_stack([processed_reactant_irc, processed_product_irc])

    where_ts = np.where(full_irc[:,0] == 0)[0]
    if len(where_ts) == 2:
        full_irc = np.delete(full_irc, where_ts[0], axis=0)
    
    return full_irc

def calculate_relative_energies(irc_calculation):
    '''
    This function takes the processed IRC calculation (as 2-column np_array)
    and appends a third column containing the relative energy to the TS (Ets = 0).
    Returns a np_array
    '''
    maximum_energy = np.amax(irc_calculation[:,1])
    relative_energies = list()
    for energy in irc_calculation[:,1]:
        relative_energy = 627.5*(energy - maximum_energy)
        relative_energies.append(relative_energy)
    relative_energies = np.array(relative_energies)
    return np.column_stack([irc_calculation, relative_energies])

def plot_irc(irc_calculation):
    '''
    This function takes a x value (the IRC) and y value (Total or Relative Energy) and returns a 2D scatter plot
    using predefined settings. The expected input is a 2-column numpy array containing the x an y values. 
    '''
    ts_position = int(np.where(irc_calculation[:,0] == 0.0)[0])
    full_irc = irc_calculation
    full_irc_without_ts = np.delete(full_irc, ts_position, axis=0)

    # plot
    fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(9,5))
    ax.scatter(
        x = full_irc_without_ts[:,0],
        y = full_irc_without_ts[:,1],
        marker='o',
        color = 'blue',
        edgecolors='black',
        alpha=0.75
    )

    if not args.downhill:
        ax.scatter(
            x = full_irc[ts_position,0],
            y = full_irc[ts_position,1],
            marker='X',
            s = 96,
            color = 'gold',
            edgecolors='black',
            alpha=0.75
    )
    else:
        ax.scatter(
            x = full_irc[ts_position,0],
            y = full_irc[ts_position,1],
            color = 'blue',
            edgecolors='black',
            alpha=0.75
    )

    ax.set_xlabel('IRC (amu$^{1/2}$-Bohr)')
    ax.set_ylabel('Relative Total Energy (kcal/mol)')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.minorticks_on()

    if args.include_title and args.title:
        ax.set_title(args.title, fontweight='bold', fontsize=16, pad=30)

    if not args.downhill:
        ax.annotate(
            'TS',(0,0), 
            xytext=(10,25), textcoords='offset points', fontsize=12, fontweight='bold', 
            arrowprops=dict(arrowstyle='-',facecolor='black')
            )

    ax.annotate(
        'A',(full_irc[0,0],full_irc[0,1]), 
        xytext=(10,25), textcoords='offset points', fontsize=12, fontweight='bold', 
        arrowprops=dict(arrowstyle='-',facecolor='black')
        )

    ax.annotate(
        'B',(full_irc[-1,0],full_irc[-1,1]), 
        xytext=(10,25), textcoords='offset points', fontsize=12, fontweight='bold', 
        arrowprops=dict(arrowstyle='-',facecolor='black')
        )

    
    plt.tight_layout
    
    if not args.title:
        plt.show()
    else:
        plt.savefig(f"{args.title.replace(' ','_')}.png", dpi=1200, bbox_inches="tight")


if not args.downhill:
    if args.subparser == 'single':
        irc_extracted = extract_irc_from_output(args.output)
        irc_rv, irc_fw = split_singlecalc_irc(irc_extracted)
    elif args.subparser == 'double':
        irc_rv = extract_irc_from_output(args.reverse)
        irc_fw = extract_irc_from_output(args.forward)

    if args.product == 'rv':
        full_irc = join_partial_irc(irc_fw,irc_rv)
    elif args.product == 'fw':
        full_irc = join_partial_irc(irc_rv,irc_fw)
        
    full_irc_relative = calculate_relative_energies(full_irc)
    plot_irc(full_irc_relative[:,[0,2]])
else:
    irc_extracted = extract_irc_from_output(args.output)
    full_irc_relative = calculate_relative_energies(irc_extracted)
    plot_irc(full_irc_relative[:,[0,2]])


