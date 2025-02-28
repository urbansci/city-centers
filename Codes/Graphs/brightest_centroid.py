# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2024-12

import geopandas as gpd
import pandas as pd
import rasterio
import numpy as np
from shapely.geometry import Point
import rasterio.mask
from geopy.distance import geodesic
from tqdm import tqdm
import Parameters.plot_utils as utils
import seaborn as sns
import matplotlib.pyplot as plt
import Parameters.consts as const

"""
Module for the generation of the ``Comparison between the identified centers with the brightest grids and centroids`` (Extended Data Figure 6).
"""


def GetBrightestAndCentroid(areas, raster):
    """
    Extract the brightest point and centroid of NTL within each urban area.

    Args:
        areas (DataFrame): The dataframe containing the extent ('geometry') of urban areas.
        raster (rasterio.DatasetReader): The NTL raster.

    Returns:
        None. Append two geometry columns to the urban areas dataframe, 'Brightest' : position of brightest points;
        'Centroid' : position of the centroids.
    """
    for idx, row in tqdm(areas.iterrows(), total=areas.shape[0], desc="Processing areas"):
        polygon = row['geometry']
        out_image, out_transform = rasterio.mask.mask(raster, [polygon], crop=True, all_touched=True)
        out_image = out_image[0]
        mask = out_image != raster.nodata

        if mask.sum() == 0:
            print(idx, "Error!")
            return None, None

        # Brightest
        brightest_index = np.unravel_index(np.argmax(out_image, axis=None), out_image.shape)
        brightest_coords = out_transform * (brightest_index[1], brightest_index[0])
        areas.loc[idx, 'Brightest'] = Point(brightest_coords[0], brightest_coords[1])

        # Centroid
        rows, cols = np.where(mask)
        coords = np.array([out_transform * (col, row) for row, col in zip(rows, cols)])
        pixel_values = out_image[rows, cols]
        weighted_x = np.sum(coords[:, 0] * pixel_values) / np.sum(pixel_values)
        weighted_y = np.sum(coords[:, 1] * pixel_values) / np.sum(pixel_values)
        areas.loc[idx, 'Centroid'] = Point(weighted_x, weighted_y)


def CalculateGeodesicDistance(row, des):
    """
    Calculate the geodesic distance under WGS84 between main center and des in row.
    """
    center = (row['geometry_'].y, row['geometry_'].x)  # (lat, lon)
    destination = (row[des].y, row[des].x)  # (lat, lon)

    return geodesic(center, destination).km


def CalculateDistance():
    """
    Calculate the distance (km) from the identified main center to the brightest point or centroid of the NTL for each urban area.

    Args:
        None.

    Returns:
        None. Append two columns to the original city shapefile, 'Dist_Brig': distance to the brightest point;
        'Dist_Cent': distance to the centroid.
    """
    areas = gpd.read_file(const.FILE_PATH + fileItems[0])
    centers = gpd.read_file(const.FILE_PATH + fileItems[1])
    raster = rasterio.open(const.FILE_PATH + fileItems[2])
    centers = (centers[centers['Is_MC'] == '1'])[['Country', 'Cluster_ID', 'geometry']]

    areas = pd.merge(areas, centers, left_on=['Country', 'ID'], right_on=['Country', 'Cluster_ID'], how='left',
                     suffixes=("", "_"))

    # Append the geometry columns
    GetBrightestAndCentroid(areas, raster)

    # Calculate the distance
    areas['Dist_Brig'] = areas.apply(CalculateGeodesicDistance, axis=1, des='Brightest')
    areas['Dist_Cent'] = areas.apply(CalculateGeodesicDistance, axis=1, des='Centroid')

    gpd.GeoDataFrame(areas[['Centroid', 'ID', 'Country']], geometry='Centroid', crs="EPSG:4326").to_file(
        const.FILE_PATH + fileItems[4])
    gpd.GeoDataFrame(areas[['Brightest', 'ID', 'Country']], geometry='Brightest', crs="EPSG:4326").to_file(
        const.FILE_PATH + fileItems[5])

    areas.drop(columns=['geometry_', 'Centroid', 'Brightest'], inplace=True)
    gpd.GeoDataFrame(areas, geometry='geometry', crs="EPSG:4326").to_file(const.FILE_PATH + fileItems[3])


def KDE():
    """
    Plot the kernel density estimation of the distance.
    """

    plt.rcParams['font.family'] = 'Arial'
    f, ax = plt.subplots(1, 2, figsize=(11 / 2.54, 5 / 2.54), dpi=600)

    areas = gpd.read_file(const.FILE_PATH + fileItems[3])
    dist = {'Brightest grid': areas['Dist_Brig'], 'Centroid': areas['Dist_Cent']}
    colors = ['#31657E', '#9B3524']

    for idx, (key, series) in enumerate(dist.items()):

        data = [series, series[areas[areas['Area'] > 200].index]]
        df = pd.DataFrame({
            1: series,
            2: series[areas[areas['Area'] > 200].index]
        })
        # Plot the kde of distance
        for i in range(len(data)):
            sns.kdeplot(ax=ax[idx], data=data[i], fill=False, color=colors[i], clip=(0, 10), linewidth=0.8, zorder=5)
            ax[idx].plot([data[i].mean()] * 2, [0, 0.02], color=colors[i], linestyle='solid', linewidth=0.8, zorder=10)
            print(key, i, data[i].mean())

        ax[idx].set_xlabel("")
        ax[idx].set_ylabel("")

        ax[idx].set_title(key, fontsize=utils.TITLE_SIZE, pad=utils.TITLE_PAD)
        utils.FormatAxis(ax[idx])
        ax[idx].set_xlim(0, 10)
        ax[idx].set_ylim(0, 0.3)

        # Inset axes to plot the bar plot
        axins = ax[idx].inset_axes([0.6, 0.3, 0.4, 0.6])
        sns.boxplot(
            data=df,
            ax=axins,
            width=0.7,
            linewidth=0.5,
            showfliers=False,
            fill=False,
            palette=colors,
            zorder=5
        )

        axins.grid(axis='y', which='major', color='#C9CACA', linewidth=0.5, zorder=1)
        axins.set_ylim(0, 25)
        axins.set_yticks([0, 10, 20])
        axins.set_xlabel("")

        utils.FormatAxis(axins, positions=[], x_ticks=False, y_ticks=False, x_labels=False)

    f.text(0.02, 0.5, f'Density', fontsize=utils.LABEL_SIZE, rotation='vertical', ha='center',
           va='center')
    f.text(0.5, 0.04, 'Distance to the main center (km)', fontsize=utils.LABEL_SIZE, ha='center', va='center')

    ax[0].legend(loc='best', fontsize=utils.LEGEND_SIZE, frameon=False)
    ax[1].legend().remove()

    plt.subplots_adjust(top=0.9, right=0.96, left=0.1, bottom=0.15, wspace=0.2, hspace=0.2)
    plt.savefig("./Graph/Distance.png")

fileItems = ["Result/Shp/Global/GlobalAreas.shp", "Result/Shp/Global/GlobalCenters.shp",
             "Data/Bright/World/VNL_v2_npp_2020_global_vcmslcfg_c202102150000.average_masked.tif",
             "Validation/Brig_Cent/Distance.shp", "Validation/Brig_Cent/Centroid.shp",
             "Validation/Brig_Cent/Brighest.shp"]

if __name__ == '__main__':
    CalculateDistance()
    KDE()
