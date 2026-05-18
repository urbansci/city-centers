# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2026-04

import os.path
from osgeo import gdal, ogr
import geopandas as gpd
from shapely.geometry import Polygon, LineString
from pyproj import Geod
from tqdm import tqdm
import numpy as np
from scipy.ndimage import gaussian_filter
from geopandas import GeoDataFrame
import Parameters.consts as const

"""
Module generating the contours map of the smoothed NTL.

This module first conducts the Gaussian convolution, then generates the contours map
and calculates the enclosed area of contours.
"""

class ContoursGenerator():
    """
    A class generating the contours map of the country/region.
    """
    def __init__(self):
        return

    def GaussianConvolution(self, NTLPath, smoothNTLPath, radius, sigma, tile_size=4096):
        """Conduct gaussian smooth to the NTL using tiled processing for large rasters.

        Note:
            We do not apply filter to the nodata pixels, which are commonly water bodies or
            located beyond the national boundary.

        Args:
            NTLPath (str): The file path of the NTL raster.
            smoothNTLPath (str): The file path of the smoothed NTL raster.
            radius (int): Radius of the Gaussian kernel.
            sigma (float): Standard deviation for the Gaussian kernel.
            tile_size (int): Size of each processing tile (default 4096).

        Returns:
            None: Save the smoothed NTL as a GeoTIFF.
        """
        import rasterio
        from rasterio.windows import Window

        # Padding size: the kernel extends radius*sigma pixels in each direction
        pad = int(np.ceil(radius * sigma)) + 1

        with rasterio.open(NTLPath) as src:
            meta = src.meta.copy()
            meta.update(dtype='float32', driver='GTiff')
            nodata = src.nodata
            height, width = src.height, src.width

            with rasterio.open(smoothNTLPath, 'w', **meta) as dst:
                for row_off in tqdm(range(0, height, tile_size), desc="Gaussian smoothing"):
                    for col_off in range(0, width, tile_size):
                        # Determine the padded read window (clamped to raster bounds)
                        read_row = max(row_off - pad, 0)
                        read_col = max(col_off - pad, 0)
                        read_h = min(row_off + tile_size + pad, height) - read_row
                        read_w = min(col_off + tile_size + pad, width) - read_col
                        read_window = Window(read_col, read_row, read_w, read_h)

                        tile = src.read(1, window=read_window).astype(np.float64)

                        # Apply Gaussian filter on the padded tile
                        nodata_mask = (tile == nodata) if nodata is not None else np.zeros_like(tile, dtype=bool)
                        smoothed = gaussian_filter(tile, sigma, radius=radius)
                        smoothed[nodata_mask] = nodata if nodata is not None else 0

                        # Trim padding to extract the core tile
                        top_pad = row_off - read_row
                        left_pad = col_off - read_col
                        core_h = min(tile_size, height - row_off)
                        core_w = min(tile_size, width - col_off)
                        core = smoothed[top_pad:top_pad + core_h, left_pad:left_pad + core_w]

                        # Write the core tile
                        write_window = Window(col_off, row_off, core_w, core_h)
                        dst.write(core.astype(np.float32), 1, window=write_window)

    def GenerateContours(self, smoothNTLPath, contourPath):
        """Generate the contours map of the smoothed NTL.

        Args:
            smoothNTLPath (str): The file path of the smoothed NTL raster.
            contourPath (str): The file path of the output contours (.gpkg).

        Returns:
            None: Save the contours to GeoPackage.
        """
        input_ds = gdal.Open(smoothNTLPath)
        if input_ds is None:
            raise RuntimeError("Cannot open the raster!")

        driver = ogr.GetDriverByName("GPKG")
        output_ds = driver.CreateDataSource(contourPath)
        srs = ogr.osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        layer = output_ds.CreateLayer("contours", srs=srs, geom_type=ogr.wkbLineString)
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
            empty_gdf.to_file(contourPath, driver='GPKG', layer='contours', OVERWRITE='YES')
        else:
            contours_filtered.to_file(contourPath, driver='GPKG', layer='contours', OVERWRITE='YES')

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