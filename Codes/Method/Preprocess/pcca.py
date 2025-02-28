# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2024-12

import os
from math import fabs, sin, cos, sqrt, pi
from osgeo import gdal
import matplotlib.pyplot as plt
import numpy as np
from math import log
from tqdm import tqdm
import scipy
import Parameters.consts as const
import Parameters.plot_utils as utils

"""
Module generating the percolation graph for each country/region.

This module generates the percolation graph for each country/region based on the ``Percolation-based City Clustering Algorithm``.
"""

def CountSpheroidArea(lat, dLat, dLng):
    """
    Calculate the area of the grid cell.

    Args:
        lat (float): The central latitude of the grid cell, in the unit of radians.
        dLat (float): The latitude interval of the grid cell, in the unit of radians.
        dLng (float): The longitude interval of the grid cell, in the unit of radians.

    Returns:
        float: The area of the grid cell, in the unit of square kilometer.
    """
    a = 6378.137
    b = 6356.752314245179
    e = sqrt(1 - (b / a) * (b / a))
    tmp = 1 - e * e * sin(lat) * sin(lat)
    m = a * (1 - e * e) / (tmp * sqrt(tmp))
    n = a / sqrt(tmp)
    area = dLat * dLng * m * n * cos(lat)
    return fabs(area)


def BinarizeTiff(NTLArray, thresholds):
    """
    Binarize the raw data to the urban/non-urban data under each potential threshold.

    Args:
        NTLArray (ndarray): The array representation of the NTL raster.
        thresholds (ndarray): The range of thresholds to be processed.

    Returns:
        dict: A dictionary where keys are the thresholds and values are 2D numpy arrays
            representing the binarized NTL data. The value is 1 for urban and 0 for non-urban.
    """

    units = {}  # < thres: binarized array >
    for thres in thresholds:
        units[thres] = np.zeros(NTLArray.shape, dtype='float32')

    # Binarize the tiff data to the urban/non-urban data under each potential threshold
    for i in range(NTLArray.shape[0]):
        for j in range(NTLArray.shape[1]):
            for thres in thresholds:
                if NTLArray[i][j] <= thres:    # thres is in ascending order
                    break
                units[thres][(i, j)] = 1
    return units


def Aggregate(units, thresholds):
    """
    Merge urban cells into 8-connected urban clusters.

    Args:
        units (dict): A dictionary where keys are the thresholds and values are 2D numpy arrays
            representing the binarized NTL data. The value is 1 for urban and 0 for non-urban.
        thresholds (ndarray): The range of thresholds to be processed.

    Returns:
        dict: A dictionary where keys are thresholds and values are dictionaries representing the
            collections of grid cells for each cluster.
    """
    structure = np.array([[1, 1, 1],
                          [1, 1, 1],
                          [1, 1, 1]])
    clusters = {}  # <thres: <cluster_id: list[cluster_units]> >
    for thres in thresholds:
        _unit = units[thres]
        labelArray, labelNum = scipy.ndimage.label(_unit, structure=structure)
        _cluster = { i : [] for i in np.arange(1, labelNum+1)}  # <cluster_id: list[cluster_units,...]>
        for i in range(_unit.shape[0]):
            for j in range(_unit.shape[1]):
                if labelArray[i][j] != 0:
                    _cluster[labelArray[i][j]].append((i, j))
        clusters[thres] = _cluster

    return clusters

def CalculateAreas(clusters, thresholds):
    """
    Calculate the spherical areas for the clusters under each potential threshold.

     Args:
         clusters (dict): A dictionary where keys are thresholds and values are dictionaries representing the
             collections of grid cells for each cluster.
        thresholds (ndarray): The range of thresholds to be processed.

     Returns:
        dict: A dictionary where keys are thresholds and values are dictionaries representing the area of each cluster.
        dict: A dictionary where keys are thresholds and values are the proportion of the area occupied by the largest cluster.
        dict: A dictionary where keys are thresholds and values are the entropy of the area distribution of clusters.
     """

    clustersArea = {
        thres: {idx: sum(CountSpheroidArea(lat=(topX + (i + 0.5) * pixelHeight) * pi / 180.0, dLat=dLat, dLng=dLng)
                   for i, j in cluster) for idx, cluster in clusters[thres].items()}
        for thres in thresholds} # <thres: <cluster_id, cluster area>>

    largestArea = {
        thres: 0 if (total := sum(clustersArea[thres].values())) == 0
        else max(clustersArea[thres].values()) / total
        for thres in thresholds}  # <thres: largest cluster area portion)>

    entropies = {
        thres: 0 if (total := sum(clustersArea[thres].values())) == 0
        else sum([- (area / total) * log(area / total) for area in clustersArea[thres].values()])
        for thres in thresholds}  # <thres: cluster areas entropy>

    return clustersArea, largestArea, entropies

def PercolatioinGraph(largestArea, thresholds, optThres=0.5):
    """
    Plot the percolation graph for each country/region.

    Args:
        largestArea (dict): A dictionary where keys are thresholds and values are the proportion of the area
            occupied by the largest cluster.
        optThres (float, optional): The manually identified optimal threshold for the country/region.
        thresholds (ndarray): The range of thresholds to be processed.

    Returns:
         None. Save the figure to specified path.
    """
    if max(thresholds) == 10:
        fig, ax = plt.subplots(1, 1, figsize=(11.0 /2.54, 8.0 /2.54), dpi=600)
    else:
        fig, ax = plt.subplots(1, 1, figsize=(11.0 * 3/ 2.54, 8.0 / 2.54), dpi=600)


    ax.plot(const.THRESHOLDS, largestArea.values(), color='#305F92', linewidth=0.8, marker='o', markersize=5, markeredgecolor='grey',
            markeredgewidth=0.3)
    ax.axvline(x=optThres, color='k', linewidth=1)

    ax.set_xlabel('Threshold', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_ylabel('Largest cluster proportion', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_title(f'{country} : {const.NAMES[country]}', fontsize=utils.TITLE_SIZE, pad=utils.TITLE_PAD)
    ax.set_xticks(ticks=range(11))
    ax.yaxis.grid(color='#cccccc', linewidth=0.25, alpha=1)

    utils.FormatAxis(ax)

    plt.subplots_adjust(top=0.90, right=0.98, left=0.10, bottom=0.10, wspace=0, hspace=0)
    plt.savefig(const.FILE_PATH + fileItems[1])

    plt.cla()
    plt.close("all")

def ReadRaster(path):
    """
    Read the raster specified by file path.

    Args:
        path (str): The file path of the raster.

    Returns:.
        array (ndarray): The array representation of the raster.
        topX (float): The top x coordinates.
        pixelHeight (float): The height of each pixel.
        dLat (float): The latitude interval of the grid cell, in the unit of radians.
        dLng (float): The longitude interval of the grid cell, in the unit of radians.
    """
    if not os.path.exists(path):
        print(f'{country} donnot have NTL data!')
        return

    tiff = gdal.Open(path)
    band = tiff.GetRasterBand(1)
    leftY, pixelWidth, _, topX, _, pixelHeight = tiff.GetGeoTransform()
    dLat = pixelHeight * pi / 180.0
    dLng = pixelWidth * pi / 180.0
    array = band.ReadAsArray(0, 0, tiff.RasterXSize, tiff.RasterYSize)

    return array, topX, pixelHeight, dLat, dLng

def ExecutebyOneThreshold():
    """
    Iterate one threshold in each time to generate the percolation data.

    For large countries/regions, we can not process all thresholds simultaneously because of the insufficient memory.
    So in each time, we only process single threshold and storage its results.

    Args:
        None.

    Returns:
        dict: A dictionary where keys are thresholds and values are the proportion of the area occupied by the largest cluster.
    """
    TclustersArea = {}
    TlargestArea = {}
    Tentropies = {}  # Temporary variables to storage the result of each iteration
    for time in tqdm(np.arange(len(const.THRESHOLDS))):
        Tthresholds = np.arange(const.THRESHOLDS_STEP * time, const.THRESHOLDS_STEP * (time + 1),
                                const.THRESHOLDS_STEP)
        units = BinarizeTiff(NTLArray, Tthresholds)
        clusters = Aggregate(units, Tthresholds)
        clustersArea, largestArea, entropy = CalculateAreas(clusters, Tthresholds)

        # Update data
        TclustersArea.update(clustersArea)
        TlargestArea.update(largestArea)
        Tentropies.update(entropy)

    # Sort data
    clustersArea = dict(sorted(TclustersArea.items(), key=lambda x: x[0]))
    largestArea = dict(sorted(TlargestArea.items(), key=lambda x: x[0]))
    entropies = dict(sorted(Tentropies.items(), key=lambda x: x[0]))

    return largestArea

def ExecutebyAllThreholds():
    """
    Process the thresholds range in one time.

    Args:
        None.

    Returns:
        dict: A dictionary where keys are thresholds and values are the proportion of the area occupied by the largest cluster.
    """

    units = BinarizeTiff(NTLArray, const.THRESHOLDS)
    clusters = Aggregate(units, const.THRESHOLDS)
    clustersArea, largestArea, entropy = CalculateAreas(clusters, const.THRESHOLDS)

    return largestArea


if __name__ == '__main__':
    for country in range(233):
        print('Country : ', country)

        fileItems = [f'Data\\Bright\\Country\\Light{country}.tif',
                     f'Result\\PCCA\\Final\\pcca{country}',
                     f'Temp\\portionR{country}']

        NTLArray, topX, pixelHeight, dLat, dLng = ReadRaster(const.FILE_PATH + fileItems[0])
        clustersArea = ExecutebyOneThreshold()
        PercolatioinGraph(clustersArea, const.THRESHOLDS, const.OPTIMAL_THRESHOLD[country] if const.OPTIMAL_THRESHOLD[country] != None else 0.5)

