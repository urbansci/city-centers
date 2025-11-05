# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2025-11

import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rasterio
import Parameters.plot_utils as utils
import Parameters.consts as const
from pyproj import Geod
from rasterio.features import geometry_mask
from tqdm import tqdm
from scipy.optimize import curve_fit
from rasterio.features import geometry_window

"""
Module for the generation of the ``Population-weighted average distance to the center'' plot (Figure 5b).
"""

def average_distance_of_city(pop_raster, poly_geom, points_geoms):
    """
       Calculate the weighted average distance from population raster pixels to the city center points.

       Args:
           pop_raster (rasterio.io.DatasetReader): The population raster representing the population count of each pixel.
           poly_geom (shapely.geometry.Polygon): The polygon geometry representing the boundary of the urban area.
           points_geoms (iterable of shapely.geometry.Point): An iterable of `shapely.geometry.Point` objects representing
                                                              the coordinates of the city center points.

       Returns:
           float: The weighted average distance to the center points, in kilometers.
    """
    site_coords = np.array([(p.x, p.y) for p in points_geoms])

    window = geometry_window(pop_raster, [poly_geom])
    vals = pop_raster.read(1, window=window, masked=True)
    mask = geometry_mask([poly_geom], transform=pop_raster.window_transform(window), out_shape=vals.shape)
    vals.mask |= mask
    rows, cols = np.where(~vals.mask)

    xs, ys = pop_raster.xy(rows + window.row_off, cols + window.col_off)
    pixel_coords = np.column_stack([xs, ys])
    pixel_vals = vals[~vals.mask]
    if np.sum(pixel_vals) == 0 or len(site_coords) == 0:
        return None

    geod = Geod(ellps="WGS84")
    if len(site_coords) == 1:
        lon, lat = site_coords[0]
        dist = geod.inv(
            np.full(len(pixel_coords), lon),
            np.full(len(pixel_coords), lat),
            pixel_coords[:, 0],
            pixel_coords[:, 1]
        )[2]
    else:
        dist = np.min([
            geod.inv(
                np.full(len(pixel_coords), lon),
                np.full(len(pixel_coords), lat),
                pixel_coords[:, 0],
                pixel_coords[:, 1]
            )[2]
            for lon, lat in site_coords
        ], axis=0)

    return np.sum(pixel_vals * dist) / (1e3 * np.sum(pixel_vals))

def plot_actuality_vs_mono(ax):
    """
    Plot the average distance relative to the center number of cities.
    """
    def pop_weighted_distance():
        """
        Calculate the population-weighted distance to the nearest center/to the main center for each urban area.
        """
        points = gpd.read_file(const.FILE_PATH + fileItems[0])
        polys = gpd.read_file(const.FILE_PATH + fileItems[1])
        raster = rasterio.open(const.FILE_PATH + fileItems[2])
        points.to_crs(4326)
        polys.to_crs(4326)

        results_all = []
        results_main = []
        for _, poly in tqdm(polys.iterrows(), total=len(polys)):
            sites = points[points['ID_UC_G0'] == poly['ID_UC_G0']]
            main_site = sites[sites['Main'] == 1].iloc[0]

            results_all.append(average_distance_of_city(raster, poly.geometry, sites.geometry))
            results_main.append(average_distance_of_city(raster, poly.geometry, [main_site.geometry]))

        polys['dist_all'] = results_all
        polys['dist_main'] = results_main
        # polys.to_file(const.FILE_PATH + fileItems[3])
        return polys

    def errorbar_plt(ax, data, marker, color, label):
        """
        Errorbar plot to show the distribution of distance for urban areas grouped by the center number.
        """
        ax.errorbar(data.mean().index, data.mean().values, yerr=data.std().values, fmt=marker, markersize=4,
                    markerfacecolor=color, color='black', markeredgewidth=0.3, capsize=2, capthick=0.25, elinewidth=0.25,
                    label=label, zorder=10)

    polys = pop_weighted_distance().dropna(subset='dist_all')
    # group cities by center number, and aggregate those with 20+ center.
    bins = list(range(1, 21, 1))
    bins.append(np.inf)
    labels = list(range(1, 21))
    polys['Center_Num_bin'] = pd.cut(polys['Center_Num'], bins=bins, labels=labels, right=False)
    group = polys.groupby('Center_Num_bin', observed=False)
    dist = group['dist_all']
    dist_main = group['dist_main']

    errorbar_plt(ax, dist, 's', colors[0], 'Nearest center')
    errorbar_plt(ax, dist_main, '^', colors[1], 'Main center')

    print(f"Average distance to nearest center: {dist.mean().mean()} \n "
          f"Average distance to main center: {dist_main.mean().mean()} \n")

    # Statistics
    ax.hlines(dist.mean().mean(), *ax.get_xlim(), color=colors[0], lw=0.5, zorder=20)

    def power_law(x, a, b):
        return a * x ** b

    x = dist_main.mean().index.astype(float)
    y = dist_main.mean().values
    params, cov = curve_fit(power_law, x, y)
    a, b = params
    perr = np.sqrt(np.diag(cov))
    y_pred = power_law(x, a, b)
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot

    print(f"a = {a:.4f} ± {perr[0]:.4f}\n"
          f"b = {b:.4f} ± {perr[1]:.4f}\n"
          f"r2 = {r2:.4f}")
    x_fit = np.linspace(0, 21, 200)
    y_fit = power_law(x_fit, a, b)
    ax.plot(x_fit, y_fit, color=colors[1], lw=0.5, zorder=20)

    ax.set_xlabel('Center number', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_ylabel(f'Population-weighted average distance (km)', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)

    ax.set_xlim((0, 21))
    ax.set_xticks(ticks=[1, 5, 10, 15, 20], labels=[1, 5, 10, 15, '20+'])
    ax.set_ylim((0, 30))
    ax.minorticks_off()
    ax.legend().remove()
    utils.FormatAxis(ax)


fileItems = [r'GHSL\Global\globalCenters_edited.shp', r'GHSL\Global\globalClusters.shp',
             r'Data\WorldPop\World\ppp_2020_1km_Aggregated.tif']
colors = ['#CF82D9', '#1DA0CC']

if __name__ == '__main__':
    fig, ax = plt.subplots(1, 2, figsize=(18 / 2.54, 6 / 2.54), dpi=600)
    plot_actuality_vs_mono(ax[1])
    ax[0].axis('off')

    ax[0].annotate("a", xy=(-0.1, 1.05), xycoords="axes fraction", fontsize=utils.ORDER_SIZE, weight="bold")
    ax[1].annotate("b", xy=(-0.1, 1.05), xycoords="axes fraction", fontsize=utils.ORDER_SIZE, weight="bold")

    plt.subplots_adjust(top=0.9, bottom=0.1, left=0.05, right=0.98, wspace=0.1, hspace=0.2)

    plt.savefig(r'D:\Research\CityCenter\Statistics\Accessibility\dist_bin.pdf')