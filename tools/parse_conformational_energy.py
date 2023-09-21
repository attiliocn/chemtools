#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
from natsort import natsort_keygen

args = argparse.ArgumentParser()
args.add_argument(
    '-f', '--file', help='.csv file containing the energies to be parsed'
)
args.add_argument(
    '-c', '--energy-column', help='0-indexed column number containing the energy to be parsed', type=int
)
args.add_argument(
    '--sort', help='Sort each subset by energy', action='store_true'
)
args = args.parse_args()

def get_relative_energy_Eh(energies):
    relative_energies = (energies - energies.min())
    return relative_energies

def get_boltzmann_weights(relative_energies, temperature=25):
    GAS_CONSTANT = 8.314
    TEMPERATURE_K = 273.15 + temperature
    RT = (GAS_CONSTANT * TEMPERATURE_K) / 1000 # RT in kJ/mol
    RT = RT / 4.184 # RT in kcal/mol
    relative_energies = relative_energies * 627.5096 # relative energies in kcal/mol
    exponential_terms = np.exp(-relative_energies/RT)
    weights = exponential_terms / exponential_terms.sum()
    return weights

dataset = pd.read_csv(args.file)

basenames = []
for file in dataset.iloc[:,0]:
    basenames.append(f"{file.rsplit('_', 1)[0]}_")
dataset['basename'] = basenames

relative_energies_parsed = {}
boltzmann_weights_parsed = {}
boltzmann_weights_acc_parsed = {}
conformers_up_90_parsed = {}
for basename in dataset['basename'].unique():
    subset = dataset.loc[(dataset['basename'] == basename)]
    energies = subset.iloc[:,args.energy_column]

    relative_energies = get_relative_energy_Eh(energies)
    for i, energy in relative_energies.items():
        relative_energies_parsed[i] = energy * 627.5096

    boltzmann_weights = get_boltzmann_weights(relative_energies)
    for i, weight in boltzmann_weights.items():
        boltzmann_weights_parsed[i] = weight

    boltzmann_weights_acc = boltzmann_weights.sort_values(ascending=False).cumsum()
    for i, weight in boltzmann_weights_acc.items():
        boltzmann_weights_acc_parsed[i] = weight

    conformers_bellow_90 = (boltzmann_weights_acc < 0.9)
    conformers_first_above_90 = conformers_bellow_90.argmin()
    conformers_up_90 = conformers_bellow_90.copy()
    conformers_up_90.iloc[conformers_first_above_90] = True
    for i, weight in conformers_up_90.items():
        conformers_up_90_parsed[i] = weight
        
relative_energies_parsed = pd.Series(relative_energies_parsed, name='Relative Energy kcal/mol')
boltzmann_weights_parsed = pd.Series(boltzmann_weights_parsed, name='Boltzmann weight')
boltzmann_weights_acc_parsed = pd.Series(boltzmann_weights_acc_parsed, name='Boltzmann weight acc')
conformers_up_90_parsed = pd.Series(conformers_up_90_parsed, name='90% Threshold')

parsed_energy  = pd.concat(
    [
        dataset, 
        relative_energies_parsed, 
        boltzmann_weights_parsed, 
        boltzmann_weights_acc_parsed,
        conformers_up_90_parsed,
    ], 
    axis=1
)

if args.sort:
    parsed_energy.sort_values(
        by=['basename', 'Relative Energy kcal/mol'],
        key=natsort_keygen(), 
        inplace=True
    )

parsed_energy.drop('basename', axis=1, inplace=True)

parsed_energy.to_csv('parsed-energy.csv', index=False)