# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2024-12

import numpy as np
from shapely.geometry import Polygon, MultiPolygon
import scipy
from rasterstats import zonal_stats
import Parameters.consts as const
from shapely.affinity import affine_transform
from pyproj import Geod
import geopandas as gpd

"""
Module generating the urban areas system under the optimal threshold for each country/region.

This module binarizes and aggregates the NTL raster to urban areas under the optimal threshold obtained in ``pcca.py``, 
calculates their attributes including area and population, and conducts a filter based on the minimum population size 
and population density.
"""

class Cluster:
    """
    A class representing an urban area.

    This class is used to define a derived urban area that holds attributes including its
    geographical area, brightness, population, and other relevant properties.

    Args:
        ID (int): The unique identifier for the cluster.
        Area (float, optional): The area of the cluster. Defaults to 0.0.
        Brightness (float, optional): The total sum of brightness value of cells within the cluster. Defaults to 0.0.
        Population (float, optional): The total population of the cluster. Defaults to 0.0.

    Attributes:
        geometry (shapely.geometry): A geometry object to hold the geometry data of the cluster.
        ID (int): The unique identifier for the cluster.
        Area (float): The area of the cluster.
        Brightness (float, optional): The total sum of brightness value of cells within the cluster.
        Population (float): The population of the cluster.
        CenterNumber (int): The number assigned to the cluster center.
        Pixels (list): A list to hold the indices of pixels of the NTL raster within the cluster.
        Name (str): The geoname of the cluster.
    """
    def __init__(self, ID, Area=0.0, Brightness=0.0, Population=0.0):
        self.geometry = None
        self.Pixels = []
        self.ID = ID
        self.Area = Area
        self.Brightness = Brightness
        self.Population = Population
        self.CenterNumber = 0
        self.Name = ""

class UrbanAreasGenerator:
    """
    A class generating the urban areas system of a country/region.

    This class binarizes and aggregates the NTL raster to urban areas under the optimal threshold obtained in ``pcca.py``,
    calculates their attributes including area and population, and conducts a filter based on the minimum population size.

    Args:
        None.

    Attributes:
        clusters (list[Cluster]): The list of the identified urban areas.
    """

    def __init__(self):
        self.clusters = []

    def Binarize(self, NTLArray, optThres):
        """
        Binarize the raw data to the urban/non-urban cells under the optimal threshold.

        Cells with brightness values higher than the threshold will be marked as urban cells, otherwise non-urban cells.

        Args:
            NTLArray (ndarray): The array representation of the NTL raster.
            optThres (float): The optimal brightness threshold.

        Returns:
            ndarray: The binarized array under the optimal threshold with the same shape of input ``NTLArray``.
        """
        binaryArray = (NTLArray > optThres).astype(int)

        return binaryArray

    def Aggregate(self, binaryArray, transform):
        """
        Merge urban cells into urban clusters.

        Extract the 8-connected regions, i.e., the urban areas, in the binarzied raster.

        Args:
            binaryArray (ndarray): The binarized array denoting the urban/non-urban cells.
            transform (affine.Affine): The affine matrix from the raster coordinates to the geographic coordinates.

        Returns:
            None. Construct the global variable ``clusters`` with geometry information.
        """

        structure = np.array([[1, 1, 1],
                              [1, 1, 1],
                              [1, 1, 1]])
        labelArray, labelNum = scipy.ndimage.label(binaryArray, structure=structure)
        # Obtain the pixels list for each urban area
        self.clusters = [Cluster(i) for i in range(labelNum)]
        for i in range(labelArray.shape[0]):
            for j in range(labelArray.shape[1]):
                if labelArray[i][j] != 0:
                    self.clusters[labelArray[i][j]-1].Pixels.append((i, j))
        # Extract the boundary for each urban area
        for _cluster in self.clusters:
            polygons = [
                Polygon([
                (i, j),
                (i + 1, j),
                (i + 1, j + 1),
                (i, j + 1) ])
                for (j, i) in _cluster.Pixels ]

            _boundary = MultiPolygon(polygons).buffer(0)
            boundary = affine_transform(
                _boundary,
                [transform.a, transform.b, transform.d, transform.e, transform.xoff, transform.yoff]
            )  # Transform to raster coors
            _cluster.geometry = boundary

    def CalAttributes(self, NTLArray, NTLTransform, PopArray, PopTransform):
        """
        Obtain the attributes of urban areas, including area, total brightness value and total population.

        Urban areas with total population less than 2,000 are eliminated.

        Args:
            NTLArray (ndarray): The array representation of the NTL raster.
            NTLTransform (affine.Affine): The affine matrix of the NTL raster.
            PopArray (ndarray): The array representation of the Pop raster.
            PopTransform (affine.Affine): The affine matrix of the Pop raster.

        Returns:
            None. Append attributes to the global variable ``clusters``.
        """

        geod = Geod(ellps='WGS84')
        for _cluster in self.clusters:
            _cluster.Area = abs(geod.geometry_area_perimeter(_cluster.geometry)[0])/1E6  # square kilometer
            _cluster.Brightness = zonal_stats(_cluster.geometry, NTLArray, affine=NTLTransform, stats='sum', all_touched=True)[0]['sum']
            _cluster.Population = zonal_stats(_cluster.geometry, PopArray, affine=PopTransform, stats='sum', all_touched=True)[0]['sum']

        filteredClusters = [_cluster for _cluster in self.clusters if _cluster.Population >= const.MINIMUN_CITY_POPULATION
                            and _cluster.Population/_cluster.Area >= const.MINIMUN_CITY_DENSITY]
        for i in range(len(filteredClusters)):
            filteredClusters[i].ID = i  # reindex

        self.clusters = filteredClusters

    def Execute(self, optThres, NTLArray, NTLTransform, PopArray, PopTransform):
        """
        Generate the urban areas under the optimal threshold by the 'Percolation-based city cluster algorithm'.

        Args:
            optThres (float): The optimal brightness threshold.
            NTLArray (ndarray): The array representation of the NTL raster.
            NTLTransform (affine.Affine): The affine matrix of the NTL raster.
            PopArray (ndarray): The array representation of the Pop raster.
            PopTransform (affine.Affine): The affine matrix of the Pop raster.

        Returns:
            list[Cluster]: The list of ``Cluster`` objects representing the identified urban areas.
        """

        self.clusters = []

        binaryArray = self.Binarize(NTLArray, optThres)
        self.Aggregate(binaryArray, NTLTransform)
        self.CalAttributes(NTLArray, NTLTransform, PopArray, PopTransform)

        return self.clusters