# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2024-12

import geopandas as gpd
import rasterio
import numpy as np
import matplotlib.pyplot as plt
from rasterio.mask import mask
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from geopy.distance import geodesic
from geographiclib.geodesic import Geodesic
from shapely.geometry import Polygon
import Parameters.plot_utils as utils
import Parameters.consts as const

"""
Module for the generation of the ``Population density decay curve`` plot (Extended Data Figure 5).
"""

def GeodesicArea(polygon):
    """
    Calculate the geodesic area of the given polygon.

    Args:
        polygon (geometry.Polygon): The input polygon object.

    Returns:
        float: The geodesic area, with the unit of square kilometer.
    """
    geod = Geodesic.WGS84
    poly = geod.Polygon()

    for lon, lat in polygon.exterior.coords:
        poly.AddPoint(lat=lat, lon=lon)

    num, perimeter, area = poly.Compute()
    return abs(area/1e6)

def GeodesicBuffer(point, distance):
    """
    Generate the buffer zone to the given point under the given geodesic distance.

    Args:
        point (geometry.Point): The input central point.
        distance (float): The given geodesic distance, with the unit of kilometer.

    Returns:
        geometry.Polygon: The buffer polygon.
    """
    buffer = []

    for angle in range(0, 360):
        destination = geodesic(kilometers=distance).destination((point.y, point.x), angle)
        buffer.append((destination.longitude, destination.latitude))

    buffer_polygon = Polygon(buffer)
    return buffer_polygon

def GetRingsDensity(point, width, radii):
    """
    Divide a set of concentric rings around the central point and calculate the average population density of each ring.

    Args:
        point (geometry.Point): The central point.
        width (float): The width of each ring, in the unit of km.
        radii (list[float]): The radius of the outer circle for each ring, in the unit of km

    Returns:
        list[float]: The average population density of each ring.
    """

    densities = []
    raster = rasterio.open(const.FILE_PATH + fileItems[0]) # pop density raster
    transform = raster.transform
    pixel_size_x = transform[0]
    pixel_size_y = -transform[4]

    for r in radii:
        # Generate Ring
        outerCircle = GeodesicBuffer(point, r)
        innerCircle = GeodesicBuffer(point, r-width)
        ring = outerCircle.difference(innerCircle)
        ring_gdf = gpd.GeoDataFrame({'geometry': [ring]}, crs='WGS84')

        # Get Mask Array
        maskArray, maskTransform = mask(raster, ring_gdf.geometry, crop=True, all_touched=True)
        maskArray = maskArray[0]

        # Calculate Average Density
        P = 0  # Sum of population
        A = 0  # Sum of Area
        for i in range(maskArray.shape[0]):
            for j in range(maskArray.shape[1]):
                if maskArray[i, j] == raster.nodata:
                    continue

                # Calculate the intersection area of the pixel and ring, the pixel's weight is given by the intersection / pixel area
                x = maskTransform[2] + j * maskTransform[0] + i * maskTransform[1]
                y = maskTransform[5] + j * maskTransform[3] + i * maskTransform[4]
                pixel = Polygon([(x, y), (x+pixel_size_x, y), (x+pixel_size_x, y+pixel_size_y), (x, y+pixel_size_y)])

                intersection = ring.intersection(pixel)

                if intersection.geom_type == "Polygon":
                    intersectionArea = GeodesicArea(intersection)
                elif intersection.geom_type == "MultiPolygon":
                    intersectionArea = sum(GeodesicArea(poly) for poly in intersection.geoms)
                else:
                    intersectionArea = 0

                P += maskArray[i, j] * intersectionArea
                A += intersectionArea

        if P == 0:
            densities.append(1) # avoid log 0
        else:
            densities.append(P/A)
    print(densities)
    return densities

def GenerateData():
    """
    Generate the density decay data of the main centers and subcenters for the specified cities.

    Args:
        None.

    Returns:
        pd.dataframe: The dataframe storing the density data, where the row representing the city and column representing the center.
    """

    width = 2
    radii = np.arange(width, 20 + width, width)

    data = pd.DataFrame(columns=['main', 'sub1', 'sub2'], index=cities)
    centers = gpd.read_file(const.FILE_PATH + fileItems[1])
    for city in cities:
        _centers = centers[centers['City'] == city]
        data.loc[city, 'main'] = GetRingsDensity(_centers[_centers['Type'] == 'Main'].iloc[0].geometry, width, radii)
        data.loc[city, 'sub1'] = GetRingsDensity(_centers[_centers['Type'] == 'Sub1'].iloc[0].geometry, width, radii)
        data.loc[city, 'sub2'] = GetRingsDensity(_centers[_centers['Type'] == 'Sub2'].iloc[0].geometry, width, radii)

    return data

# Linear regression to avoid the effect caused by the large numbers(short distance)
def func(x, a, b):
    return a-b*x


def PlotDecayCurve(axes, data):
    """
    Draw the density decay curve and fit it by the exponential model for each city.

    Args:
        axes (plt.Axes): The axes to plot.
        data (pd.Dataframe): The density data.

    Returns:
        None.
    """

    x = np.arange(1, 21, 2)
    axes = axes.flatten()
    for i, city in enumerate(cities):
        ax = axes[i]
        print(city)
        for j, center in enumerate(data.columns):
            y = data.loc[city, center]
            ax.scatter(x, y, marker='^', facecolor='white', edgecolor=utils.COUNTRY_COLORS[j], linewidth=0.3,
                       s=10,  alpha=1, zorder=12-j)

            # Fit Curve
            popt, pcov = curve_fit(func, x, np.log10(y))
            a, b = popt
            y_pred = [func(xx, a, b) for xx in x]
            r2 = r2_score(np.log10(y), y_pred)

            print(f'{center}\n   R2:{r2:.3f}, A:{10**a:.3f}, B:{b/np.log10(np.e):.3f}')
            ax.plot([x[0], x[-1]], [10**func(x[0], a, b), 10**func(x[-1], a, b)], linestyle='solid', color=utils.COUNTRY_COLORS[j],
                    linewidth=0.8, alpha=1, zorder=10-i)

        ticks = range(0, max(x)+5, 5)
        ax.set_xticks(ticks=ticks)
        ax.set_yscale('log')

        if i==0:
            ax.set_yticks(ticks=[10 ** 3, 10 ** 3.5, 1e4])
            ax.set_yticklabels(labels=[3, 3.5, 4])
        elif i==1:
            ax.set_yticks(ticks=[10 ** 3.5, 10 ** 4, 10 ** 4.5])
            ax.set_yticklabels(labels=[3.5, 4, 4.5])
        elif i==2:
            ax.set_yticks(ticks=[10 ** 2, 10 ** 3, 10 ** 4])
            ax.set_yticklabels(labels=[2, 3, 4])
        else:
            ax.set_yticks(ticks=[10 ** 3, 10 ** 3.5, 1e4])
            ax.set_yticklabels(labels=[3, 3.5, 4])
        ax.minorticks_off()
        utils.FormatAxis(ax)

fileItems = ['Data\WorldPop\World\ppp_2020_1km_Aggregated.tif',  r'Data\\Population Density\\Centers.shp']
cities = ['Los Angeles', 'Beijing', 'Berlin', 'London']

if __name__ == '__main__':

    plt.rcParams['font.family'] = 'Arial'
    fig, axes = plt.subplots(2, 2, figsize=(11/2.54, 9/2.54), dpi=600)

    PlotDecayCurve(axes, GenerateData())

    fig.text(0.02, 0.5, r'$\log_{10}$Population density (people/km$^{2}$)', fontsize=utils.LABEL_SIZE + 1, rotation='vertical', ha='center', va='center')
    fig.text(0.5, 0.03, 'Distance to the center (km)', fontsize=utils.LABEL_SIZE + 1, ha='center', va='center')

    plt.subplots_adjust(top=0.95, right=0.98, left=0.08, bottom=0.10, wspace=0.2, hspace=0.3)
    plt.savefig(r'.\Graph\Density\DensityFunction.png')