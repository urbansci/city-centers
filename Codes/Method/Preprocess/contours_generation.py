# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2024-12

import os.path
from osgeo import gdal, ogr
import geopandas as gpd
from shapely.geometry import Polygon
from pyproj import Geod
from tqdm import tqdm
import numpy as np
from scipy.ndimage import gaussian_filter
import Parameters.consts as const

"""
Module generating the contours map of the smoothed NTL for each country/region.

For the NTL of each country/region, this module first conduct one Gaussian convolution, then generate the contours map 
and calculate the enclosed area of contours. Note that the contours are generated and winnowed according to 
the "contour interval" and "minimum contour area" in this module, and are further winnowed based on 
the urban area-specified "starting contour value" in the 'centers_generator.py'.  
"""


def GaussianConvolution(radius, sigma):
    """Conduct gaussian smooth to the NTL for each country/region.

    Note:
        We do not play filter to the nodata pixels, which are commonly water bodies or located beyond the national boundary.

    Args:
        radius (int) : Radius of the Gaussian kernel.
        sigma (float) : Standard Deviation for the Gaussian kernel.

    Returns:
        None: Save the smoothed NTl to shapefile.
    """

    input_ds = gdal.Open(const.FILE_PATH + fileItems[0])
    if input_ds is None:
        raise RuntimeError("Cannot open the raster!")

    band = input_ds.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    raster_array = band.ReadAsArray().astype(np.float64)

    mask = (raster_array == nodata)
    smoothed_array = gaussian_filter(raster_array, sigma, radius=radius)
    smoothed_array[mask] = nodata  # Re-apply the no-data mask to the smoothed array.

    # Create output raster.
    driver = gdal.GetDriverByName('GTiff')
    output_ds = driver.Create(const.FILE_PATH + fileItems[1], input_ds.RasterXSize, input_ds.RasterYSize, 1,
                              gdal.GDT_Float32)
    output_ds.SetProjection(input_ds.GetProjection())
    output_ds.SetGeoTransform(input_ds.GetGeoTransform())
    output_band = output_ds.GetRasterBand(1)
    output_band.WriteArray(smoothed_array)
    output_band.SetNoDataValue(nodata)

    # Clean up
    input_ds = None
    output_ds = None


def GenerateContours(interval):
    """Generate the contours map of the smoothed NTL.

    Args:
        interval(float): The interval of the contour value.

    Returns:
        None: Save the contours to shapefile.
    """
    input_ds = gdal.Open(const.FILE_PATH + fileItems[1])
    if input_ds is None:
        raise RuntimeError("Cannot open the raster!")

    driver = ogr.GetDriverByName("ESRI Shapefile")
    output_ds = driver.CreateDataSource(const.FILE_PATH + fileItems[2])
    srs = ogr.osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    layer = output_ds.CreateLayer("Contours", srs=srs, geom_type=ogr.wkbLineString)
    layer.CreateField(ogr.FieldDefn("Value", ogr.OFTReal))

    band = input_ds.GetRasterBand(1)
    gdal.ContourGenerate(band, interval, 0, [], 0, 0, layer, -1, 0)

    input_ds = None
    output_ds = None


def CalculateArea():
    """Calculate the enclosed area of each contour under the geographical coordinate system.

    Contours smaller than the minimum contour area are eliminated.

    Args:
        None.

    Returns:
        None: Append one 'Area' field to the shapefile.
    """

    contours = gpd.read_file(const.FILE_PATH + fileItems[2])
    if len(contours) == 0:
        return
    geod = Geod(ellps='WGS84')
    contours['Area'] = contours['geometry'].apply(
        lambda x: np.nan if x is None or len(list(x.coords)) < 3
        else abs(geod.geometry_area_perimeter(Polygon(x))[0]) / 1e+6)

    contours = contours[contours['Area'] >= const.MINIMUN_CONTOUR_AREA]

    contours.to_file(const.FILE_PATH + fileItems[2], driver='ESRI Shapefile')


if __name__ == "__main__":

    for country in tqdm(range(233)):
        fileItems = [f'\Data\Bright\Country\Light{country}.TIF',
                     f'\Data\Bright\Gaussian\Gaussian{country}.TIF',
                     f'\Data\Bright\Contour\Contour{country}.shp']

        if not os.path.exists(const.FILE_PATH + fileItems[0]):
            print(f'Country {country} do not have NTL data!')
            continue

        GaussianConvolution(radius=1, sigma=5)
        GenerateContours(interval=3)
        CalculateArea()
