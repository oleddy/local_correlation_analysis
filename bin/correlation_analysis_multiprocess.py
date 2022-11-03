import cv2
import argparse
import numpy as np
from scipy.stats import pearsonr
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
import pandas as pd
import os
import sys
import multiprocessing as mp
from functools import partial

# np.set_printoptions(threshold=sys.maxsize)

parser = argparse.ArgumentParser()
parser.add_argument('-c', help = 'input conditions table', required = True)
parser.add_argument('-w', help = 'window size', required = False)
parser.add_argument('-s', help = 'step size', required = False)
parser.add_argument('-t', help = 'threshold', required = True)

def process_dir(indexed_row, args):
    row = indexed_row[1]
    allfiles = os.listdir(row['Path'])
    output = []
    size_output = []
    outfile = '_'.join([row['Condition'], row['Marker'], row['Timepoint'], '_correlation.csv'])
    if outfile not in allfiles:
        for file in allfiles:
            if file.split('-')[-1] == 'FITC_C.tif':
                prefix = '-'.join(file.split('-')[:-1])
                print(os.path.join(row['Path'], prefix))
                redfile = prefix + '-Cy5_C.tif'
                red = cv2.imread(os.path.join(row['Path'], redfile), cv2.IMREAD_GRAYSCALE)
                green = cv2.imread(os.path.join(row['Path'], file), cv2.IMREAD_GRAYSCALE)
                if type(green) == np.ndarray:
                    if args.w:
                        window = int(args.w)
                    else:
                        window = 10

                    if args.s:
                        stepsize = int(args.s)
                    else:
                        stepsize = 1

                    #gaussian smoothing to remove noise before segmentation
                    greenblur = cv2.GaussianBlur(green, (5,5), 0)

                    _, thresh = cv2.threshold(greenblur, int(args.t), 255, cv2.THRESH_BINARY)
                    # cv2.imwrite('example_threshold.tif', thresh)
                    # cv2.imwrite('correlation.tif', cor_array)
                    dist = ndi.distance_transform_edt(thresh)
                    local_max = peak_local_max(dist, indices=False, min_distance=20, labels=thresh)
                    markers = ndi.label(local_max, structure=np.ones((3, 3)))[0]
                    labels = watershed(-dist, markers, mask=thresh)

                    cor_array = np.zeros(red.shape)
                    # p_array = np.zeros(red.shape)
                    indices = labels.nonzero()
                    for i,j in zip(indices[0], indices[1]):
                    # for i in range(red.shape[0]):
                    #     for j in range(red.shape[1]):
                        if (i < window) or (j < window) or (i >= (red.shape[0] - window)) or (j >= (red.shape[0] - window)):
                            pass #leave the value zero if the sliding window would intersect the boundary
                        else:
                            minx = i - window
                            maxx = i + window + 1
                            miny = j - window
                            maxy = j + window + 1
                            r, _ = pearsonr(red[minx:maxx, miny:maxy].flatten(), green[minx:maxx, miny:maxy].flatten())
                            cor_array[i][j] = r
                            # print(r)
                            # p_array[i][j] = p
                    cor_array = np.nan_to_num(cor_array)
                    print(np.amax(cor_array))

                    # print(labels)
                    for label in range(1, np.amax(labels) + 1):
                        label_mean = np.average(cor_array, weights = labels == label) #average only over the pixels that make up the labeled area for a given GFP object (weight = 0 if not in the boundary)
                        # print(label_mean)
                        size_output.append(np.sum(labels==label))
                        output.append(label_mean)
                else:
                    print("Failed to open " + os.path.join(row['Path'], file))
        output_df = pd.DataFrame({'Mean correlation' : output, 'Area' : size_output})
        output_df.to_csv(os.path.join(row['Path'],outfile))

if __name__ == '__main__':
    args = parser.parse_args()
    # red = np.array(Image.open(args.r))
    # green = np.array(Image.open(args.g))
    cond = pd.read_csv(args.c)
    process_dir_partial = partial(process_dir, args = args)
    with mp.Pool(processes = mp.cpu_count() - 1) as pool:
        pool.map(process_dir_partial, cond.iterrows())


    # cor_array = cor_array*255
    # cor_array = cor_array.astype(int)
    # cv2.imwrite(args.o + '_correlation.tif', cor_array)
