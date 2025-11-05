# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2025-11

import os.path
from osgeo import gdal, ogr
import geopandas as gpd
from shapely.geometry import Polygon, LineString
from pyproj import Geod
from tqdm import tqdm
import numpy as np
from scipy.ndimage import gaussian_filter
from geopandas import GeoDataFrame
from ..Parameters import consts as const

"""
Module generating the contours map of the smoothed NTL for each country/region.

For the NTL of each country/region, this module first conduct the Gaussian convolution, then generate the contours map 
and calculate the enclosed area of contours.
"""

class ContoursGenerator():
    """
    A class generating the contours map of the country/region.
    """
    def __init__(self):
        return

    def GaussianConvolution(self, NTLPath, smoothNTLPath, radius, sigma):
        """Conduct gaussian smooth to the NTL for each country/region.

        Note:
            We do not play filter to the nodata pixels, which are commonly water bodies or located beyond the national boundary.

        Args:
            NTLPath (str): The file path of the NTL of the country/region.
            smoothNTLPath (str): The file path of the smoothed NTL of the country/region.
            radius (int) : Radius of the Gaussian kernel.
            sigma (float) : Standard deviation for the Gaussian kernel.

        Returns:
            None: Save the smoothed NTl as a shapefile.
        """

        input_ds = gdal.Open(NTLPath)
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
        output_ds = driver.Create(smoothNTLPath, input_ds.RasterXSize, input_ds.RasterYSize, 1,
                                  gdal.GDT_Float32)
        output_ds.SetProjection(input_ds.GetProjection())
        output_ds.SetGeoTransform(input_ds.GetGeoTransform())
        output_band = output_ds.GetRasterBand(1)
        output_band.WriteArray(smoothed_array)
        output_band.SetNoDataValue(nodata)

        # Clean up
        input_ds = None
        output_ds = None

    def GenerateContours(self, smoothNTLPath, contourPath):
        """Generate the contours map of the smoothed NTL.

        Args:
            smoothNTLPath (str): The file path of the smoothed NTL of the country/region.
            contourPath (str): The file path of the contours of the country/region.

        Returns:
            None: Save the contours to shapefile.
        """
        input_ds = gdal.Open(smoothNTLPath)
        if input_ds is None:
            raise RuntimeError("Cannot open the raster!")

        driver = ogr.GetDriverByName("ESRI Shapefile")
        output_ds = driver.CreateDataSource(contourPath)
        srs = ogr.osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        layer = output_ds.CreateLayer("Contours", srs=srs, geom_type=ogr.wkbLineString)
        layer.CreateField(ogr.FieldDefn("Value", ogr.OFTReal))

        band = input_ds.GetRasterBand(1)
        gdal.ContourGenerate(band, 1, 0, [], 0, 0, layer, -1, 0)

        input_ds = None
        output_ds = None

    def CalculateArea(self, contourPath):
        """Calculate the enclosed area of each contour under the geographical coordinate system.

        Args:
            contourPath (str): The file path of the smoothed NTL of the country/region.

        Returns:
            None: Append one 'Area' field to the shapefile.
        """

        contours = gpd.read_file(contourPath)
        if len(contours) == 0:
            return
        geod = Geod(ellps='WGS84')

        contours_filtered = contours[
            contours['geometry'].apply(
                lambda geom: isinstance(geom, LineString) and
                             len(geom.coords) >= 4 and
                             geom.coords[0] == geom.coords[-1])].copy()

        contours_filtered['Area'] = contours_filtered['geometry'].apply(
            lambda x: abs(geod.geometry_area_perimeter(Polygon(x))[0]) / 1e+6)

        if contours_filtered.empty:
            empty_gdf = GeoDataFrame(columns=contours.columns, geometry=[], crs=contours.crs)
            empty_gdf.set_geometry('geometry', inplace=True)
            empty_gdf.to_file(contourPath, driver='ESRI Shapefile')
        else:
            contours_filtered.to_file(contourPath, driver='ESRI Shapefile')

    def Execute(self, NTLPath, smoothNTLPath, contourPath):
        """
        Generate the contours map.

        Args:
            NTLPath (str): The file path of the NTL of the country/region.
            smoothNTLPath (str): The file path of the smoothed NTL of the country/region.
            contourPath (str): The file path of the contours of the country/region.

        Return:
            None. Save the contours map as a shapefile.
        """

        self.GaussianConvolution(NTLPath, smoothNTLPath, radius=1, sigma=5)
        self.GenerateContours(smoothNTLPath, contourPath)
        self.CalculateArea(contourPath)