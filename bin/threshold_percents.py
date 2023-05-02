import pandas as pd
import argparse
import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument('-c', help = 'Conditions table', required = True)
parser.add_argument('-t', help = 'correlation threshold', required = True)
parser.add_argument('-s', help = 'minimum object size (in pixels)', required = False)
parser.add_argument('-o', help = 'output path', required = True)

if __name__ == '__main__':
    args = parser.parse_args()
    conditions = pd.read_csv(args.c)
    output_df = conditions.copy()

    if args.s:
        size_threshold = int(args.s)
    else:
        size_threshold = 0
    cor_threshold = float(args.t)

    percentages = []
    ncells = []
    for i, row in conditions.iterrows():
        try:
            data = pd.read_csv(os.path.join(row['Path'], '_'.join([row['Condition'], row['Marker'], row['Timepoint'], '_correlation.csv'])))
            data_sized = data.loc[data['Area'] >= size_threshold]
            percent_colocalized = 100.*sum(data_sized['Mean correlation'] >= cor_threshold)/len(data_sized['Mean correlation'])
            percentages.append(percent_colocalized)
            ncells.append(len(data_sized['Mean correlation']))
        except FileNotFoundError:
            print('Could not open ' + os.path.join(row['Path'], '_'.join([row['Condition'], row['Marker'], row['Timepoint'], '_correlation.csv'])))
            percentages.append(np.nan)
        except KeyError:
            print('Area data not added to ' + os.path.join(row['Path'], '_'.join([row['Condition'], row['Marker'], row['Timepoint'], '_correlation.csv'])))
            percentages.append(np.nan)
    output_df['Percent colocalized'] = percentages
    output_df['N cells'] = ncells
    output_df.to_csv(args.o)
