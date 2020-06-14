#!/usr/bin/env python3

#%%
#Import required modules

import csv, os
import matplotlib.pyplot as plt
import seaborn as sns

#Hacky way of fixing jupyter path weirdness:
def is_interactive():
    import __main__ as main
    return not hasattr(main, '__file__')

if not is_interactive():
    os.chdir('../')
#%%
#Initialize functions for loading and parsing data.
def load_csv_data(csv_location):
    slope_data_dict = dict()
    with open(csv_location) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count != 0:
                slope_data_dict[row[4]] = row[3]
            line_count+=1
    return slope_data_dict

slopes_a_e = load_csv_data('./output/nucl_LSUa_LSUe.csv')
slopes_asg_a = load_csv_data('./output/nucl_LSUASG_LSUa.csv')

common_between_both = list()
for nucl in slopes_asg_a.keys():
    if nucl in slopes_a_e.keys():
        common_between_both.append((slopes_asg_a[nucl]*-1, slopes_a_e[nucl]))
        print(nucl, float(slopes_asg_a[nucl])*-1, slopes_a_e[nucl])

print(len(common_between_both), len(slopes_asg_a))