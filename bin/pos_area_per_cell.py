import cv2
import argparse
import pandas as pd
import os
import multiprocessing as mp
from functools import partial
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-c', help = 'conditions table',required = True)
parser.add_argument('-t', help = 'threshold', required = True)

def process_dir(indexed_row, args):
    row = indexed_row[1]
    allfiles = os.listdir(row['Path'])
    n_cells = 0
    lysotracker_area = 0.
    for file in allfiles:
        if 'TRED_C' in file:
            prefix = '-'.join(file.split('-')[:-1])
            print(os.path.join(row['Path'], prefix))
            dapifile = prefix + '-DAPI_C.tif'

            red = cv2.imread(os.path.join(row['Path'], file), cv2.IMREAD_GRAYSCALE)
            dapi = cv2.imread(os.path.join(row['Path'], dapifile), cv2.IMREAD_GRAYSCALE)
            # green = cv2.imread(os.path.join(row['Path'], file), cv2.IMREAD_GRAYSCALE)
            # if type(green) == np.ndarray:
            if type(red) == np.ndarray:
                redblur = cv2.GaussianBlur(red, (5,5), 0)
                _, redthresh = cv2.threshold(redblur, int(args.t), 255, cv2.THRESH_BINARY)

                dapiblur = cv2.GaussianBlur(dapi, (5,5), 0)
                _, dapithresh = cv2.threshold(dapiblur, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

                dist = ndi.distance_transform_edt(dapithresh)
                local_max = peak_local_max(dist, indices=False, min_distance=20, labels=dapithresh)
                markers = ndi.label(local_max, structure=np.ones((3, 3)))[0]
                labels = watershed(-dist, markers, mask=dapithresh)
                n_cells += np.amax(labels)
                lysotracker_area += np.sum(redthresh)
            else:
                print("Failed to open " + os.path.join(row['Path'], file))
    return lysotracker_area/n_cells


if __name__ == '__main__':
    args = parser.parse_args()
    cond = pd.read_csv(args.c)
    output_df = cond.copy()
    process_dir_partial = partial(process_dir, args = args)
    with mp.Pool(processes = mp.cpu_count() - 1) as pool:
        output_df['Area per cell'] = pool.map(process_dir_partial, cond.iterrows())
    output_df.to_csv('baf_output.csv')
