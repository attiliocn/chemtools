#!/usr/bin/env python

import argparse
import re
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(dest='subparser')

parser_single = subparsers.add_parser('single', help='For single output IRC calculation (forward/reverse in the same output')
parser_single.add_argument('output', type=str, help='Gaussian Output File')
parser_single.add_argument('product', choices=['rv','fw'], help='Indicate whether the product is in the forward or reverse direction')
parser_single.add_argument('--downhill', action='store_true', help='Special option for downhill calculations')
parser_single.add_argument('--brute-force', action='store_true', help='Use manual extraction of irc data')

parser_double = subparsers.add_parser('double', help='For IRC in which forward and reverse were calculated separately')
parser_double.add_argument('reverse', type=str, help='Gaussian Output File (reverse)')
parser_double.add_argument('forward', type=str, help='Gaussian Output File (forward)')
parser_double.add_argument('product', choices=['rv','fw'], help='Indicate whether the product is in the forward or reverse direction')
parser_double.add_argument('--downhill', action='store_true', help='Special option for downhill calculations')
parser_double.add_argument('--brute-force', action='store_true', help='Use manual extraction of irc data')

parser_single.add_argument('--title', type=str, default=None, help='Filename and Plot Title')
parser_double.add_argument('--title', type=str, default=None, help='Filename and Plot Title')

parser_single.add_argument('--include-title', action='store_true', help='Include plot title in the figure')
parser_double.add_argument('--include-title', action='store_true', help='Include plot title in the figure')

args = parser.parse_args()

print(f"Running: {args}\n")

def extract_irc_from_output_brute_force(output_file):
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
        if steps_cumulative >= len(irc_energy_all_points):
            irc_energies.append(irc_energy_all_points[-1])
        else:
            irc_energies.append(irc_energy_all_points[steps_cumulative])

    irc_plot_data = np.column_stack([irc_coordinates, irc_energies])
    return irc_plot_data


def extract_irc_from_output(output_file):
    '''
    This function takes a Gaussian output file and extract the IRC path from it. 
    Returns a 2-column IRC (x= IRC displacement, y= Total energy)
    '''  
    irc_coordinates = list()
    irc_energies = list()

    with open(output_file) as f:
        irc_data = f.readlines()

        for i,content in enumerate(irc_data):
            irc_completed_string = re.search('Reaction path calculation complete.', content)
            irc_number_of_points = re.search('Total number of points:', content)
            ts_energy_search = re.search('(Energies reported relative to the TS energy of) *-?([0-9]+.[0-9]+)', content)

            if ts_energy_search:
                ts_energy = float(ts_energy_search.group(2)) * -1
            elif irc_completed_string:
                irc_table_start = i + 7
            elif irc_number_of_points:
                irc_table_end = i - 2
        irc_table = irc_data[irc_table_start:irc_table_end]

        for irc_point in irc_table:
            _, relative_energy, irc_coordinate = irc_point.split()
            relative_energy = float(relative_energy)
            irc_coordinate = abs(float(irc_coordinate))

            irc_coordinates.append(irc_coordinate)
            irc_energies.append(relative_energy + ts_energy)

    irc_plot_data = np.column_stack([irc_coordinates, irc_energies])

    #print("Extracted IRC Table")
    #print(irc_plot_data)
    return irc_plot_data

def split_singlecalc_irc(irc_calculation):
    for i in range(len(irc_calculation)):
        if irc_calculation[i,0] - irc_calculation[i+1,0] < 0:
            fw_limit = i
            break
    irc_fw = irc_calculation[fw_limit+1:, :]
    irc_rv = irc_calculation[0:fw_limit+1, :]
    return irc_rv, irc_fw

def join_partial_irc(reactant_irc, product_irc):
    '''
    This function takes the processed IRC calculation (as 2-column np_array) for 
    both reactant and product and returns the full IRC calculation
    '''
    # reactant
    processed_reactant_irc = reactant_irc
    processed_reactant_irc[:,0] *= -1
    sorted_reactant = np.argsort(processed_reactant_irc[:,0])
    processed_reactant_irc = processed_reactant_irc[sorted_reactant]

    # product
    processed_product_irc = product_irc
    sorted_product = np.argsort(processed_product_irc[:,0])
    processed_product_irc = processed_product_irc[sorted_product]
    # full irc
    full_irc = np.row_stack([processed_reactant_irc, processed_product_irc])

    where_ts = np.where(full_irc[:,0] == 0)[0]
    if len(where_ts) == 2:
        full_irc = np.delete(full_irc, where_ts[0], axis=0)
    
    #print("IRC for plot after processing:")
    #print(full_irc)
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



def call_plot_irc():
    if args.product == 'rv':
        full_irc = join_partial_irc(irc_fw,irc_rv)
    elif args.product == 'fw':
        full_irc = join_partial_irc(irc_rv,irc_fw)
    full_irc_relative = calculate_relative_energies(full_irc)
    plot_irc(full_irc_relative[:,[0,2]])

if args.subparser == 'single':
    if not args.brute_force:
        if not args.downhill:
            irc_extracted = extract_irc_from_output(args.output)
            irc_rv, irc_fw = split_singlecalc_irc(irc_extracted)
            call_plot_irc()
        else:
            irc_extracted = extract_irc_from_output(args.output)
            full_irc_relative = calculate_relative_energies(irc_extracted)
            plot_irc(full_irc_relative[:,[0,2]])
    else:
        if not args.downhill:
            irc_extracted = extract_irc_from_output_brute_force(args.output)
            irc_rv, irc_fw = split_singlecalc_irc(irc_extracted)
            call_plot_irc()
        else:
            irc_extracted = extract_irc_from_output_brute_force(args.output)
            full_irc_relative = calculate_relative_energies(irc_extracted)
            plot_irc(full_irc_relative[:,[0,2]])

elif args.subparser == 'double':
    if not args.brute_force:
        irc_rv = extract_irc_from_output(args.reverse)
        irc_fw = extract_irc_from_output(args.forward)
    else:
        irc_rv = extract_irc_from_output_brute_force(args.reverse)
        irc_fw = extract_irc_from_output_brute_force(args.forward) 
    call_plot_irc()

print(f"\n\n\n")


