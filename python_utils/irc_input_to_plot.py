#!/usr/bin/env python


with open('plot_input.csv') as f:
    f.readline()
    irc_files_data = f.readlines()
    i = 0
    while i < len(irc_files_data):
        input_data = irc_files_data[i].strip().split(',')
        output_1, irc_type, direction_1, prod_direction, title = irc_files_data[i].strip().split(',')
        
        if "downhill" in output_1:
            print(f"irc_extractor.py single --downhill {output_1} fw --title {title}")
            i += 1
            continue       
 
        if irc_type == 'double':
            i += 1
            output_2, _, direction_2, _, _ = irc_files_data[i].strip().split(',')
            if direction_1 == "fw":
                fw_direction = output_1
                rv_direction = output_2
            elif direction_1 == "rv":
                fw_direction = output_2
                rv_direction = output_1
            print(f"irc_extractor.py {irc_type} {rv_direction} {fw_direction} {prod_direction} --title {title}")

        elif irc_type == 'single':
            print(f"irc_extractor.py {irc_type} {output_1} {prod_direction} --title {title}")
        i += 1
