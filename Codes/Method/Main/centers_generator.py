# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2024-12

import os
import numpy as np
from tqdm import tqdm
import geopandas as gpd
from shapely.geometry import Polygon, Point
from rasterstats import zonal_stats
import rasterio
import rasterio.mask
import Parameters.consts as const
from geopandas.tools import sjoin

"""
Module generating the centers within the identified urban areas for each country/region.

This module applies the ``localized contour tree method`` and ``iterative main tree selection method`` to detect 
the centers and identify the main center within the urban area.
"""

class Contour:
    """
    A class representing a contour upon the NTL surface.

    This class is used to define a contour that holds attributes including its ID, geographical extent, area,
    and other relevant properties.

    Args:
        ID (int): The unique identifier for the contour.
        Value (int): The brightness value of the contour.
        Area (float): The enclosed area of the contour.

    Attributes:
        geometry (shapely.geometry): A geometry object to hold the geometry data of the contour.
        ID (int): The unique identifier for the contour.
        Value (int): The brightness value of the contour.
        Area (float): The enclosed area of the contour.
        CenterID (int): The ID of the center corresponding to the contour. Only used for seed contours.
        IsMainCenter (int): Whether the contour corresponds to the main center, 0 for false, 1 for true. Only used for
            seed contours.
        Childs (list): The list of ID of the childs of the contour, only direct-connected childs.
        Parent (int): The ID of the contour's parent. -1 for root contour.

    """
    def __init__(self, ID, Value, geometry, Area=0):
        self.ID = ID
        self.Value = Value
        self.geometry = geometry
        self.Area = Area
        self.Childs = []
        self.Parent = -1
        self.CenterID = -1
        self.IsMainCenter = 0


class Center:
    """
    A class representing a city center.

    This class is used to define a center that holds attributes including its geographical location, corresponding cluster,
    and whether a main center.

    Args:
        ID (int): The unique identifier for the center.
        geometry (shapely.geometry): A geometry object to hold the geometry location the center.
        Latitude (float): The latitude of the center location.
        Longitude (float): The longitude of the center location.
        ClusterID (int): The ID/order of the urban area possessing the center.

    Attributes:
        ID (int): The unique identifier for the center.
        geometry (shapely.geometry): A geometry object to hold the geometry location the center.
        Latitude (float): The latitude of the center location.
        Longitude (float): The longitude of the center location.
        ClusterID (int): The ID of the urban area possessing the center.
        IsMainCenter (int): Whether the center is the main center of the urban area, 0 for false, 1 for true.

    """
    def __init__(self, ID, geometry, Latitude, Longitude, ClusterID):
        self.ID = ID
        self.geometry = geometry
        self.Latitude = Latitude
        self.Longitude = Longitude
        self.ClusterID = ClusterID
        self.IsMainCenter = 0
        self.Name = ""


class CentersGenerator():
    """
    A class generating the centers within the identified urban areas for each country/region.

    This class applies the ``localized contour tree method`` and ``iterative main tree selection method`` to detect
    the centers and identify the main center within each identified urban area.

    Args:
        None.

    Attributes:
        clusters (list[Cluster]): A list of the input urban areas.
        seedContours (list[Contour]): A list of the identified seed contours.
        centers (list[Center]): A list of the identified centers.

    Methods:
        GetCentroid(self, geometry, NTL): Extract the centroid weighted by the NTL of the contour.
        DetectCenter(self, contourPath, NTL, NTLArray, PopArray, PopTransform): Detect the centers within the urban area and identify the main center.
        IdentifyBasecontour(self, alterBasecontours, PopArray, PopTransform): Identify the base contour from alter basecontours.
        IdentifyMainCenter(self, centers, contours, PopArray, PopTransform): Identify the main center among centers.
        AssignGeoname(self, geonamePath): Assign geonames to urban areas and its centers.
        Execute(self, Clusters, contourPath, geonamePath, NTL, NTLArray, PopArray, PopTransform): Generate the centers within the urban areas.
    """

    def __init__(self):

        self.clusters = None
        self.centers = None
        self.seedContours = None

    def GetCentroid(self, geometry, NTL):
        """Extract the centroid weighted by the NTL of the contour.

        Args:
            geometry (shapely.geometry): The geometry object of the seed contour.
            NTL (rasterio.DataReader): The NTL raster object.
        Returns:
            float: The latitude of the center.
            float: The longitude of the center.
        """

        out_image, out_transform = rasterio.mask.mask(NTL, [Polygon(geometry)], crop=True, all_touched=False)
        out_image = out_image[0]
        mask = out_image != NTL.nodata

        if mask.sum() == 0:
            raise RuntimeError("Error: No pixel within the seed contour.")

        rows, cols = np.where(mask)
        coords = np.array([out_transform * (col + 0.5, row + 0.5) for row, col in zip(rows, cols)]) # central location of pixels
        pixelValues = out_image[rows, cols]
        x = np.sum(coords[:, 0] * pixelValues) / np.sum(pixelValues)
        y = np.sum(coords[:, 1] * pixelValues) / np.sum(pixelValues)

        return y, x


    def DetectCenter(self, contourPath, NTL, NTLArray, PopArray, PopTransform):
        """
        Detect the centers within the urban area and identify the main center.

        Args:
            contourPath (str): The file path of the contours of the country/region.
            NTL (rasterio.DataReader): The NTL raster object.
            NTLArray (ndarray): The array representation of the NTL raster.
            PopArray (ndarray): The array representation of the Pop raster.
            PopTransform (affine.Affine): The affine matrix of the Pop raster.

        Returns:
            None. Initialize the global variable ``centers`` and ``seedContours`` and filter the ``clusters``.
        """

        gdf = gpd.read_file(contourPath)  # Contour
        self.clusters = sorted(self.clusters, key=lambda x: x.Area, reverse=True)
        filteredClusters = []

        for _cluster in tqdm(self.clusters):
            # Step 1: Extract all contours within the cluster
            lines = gdf[gdf.geometry.apply(lambda x: x.intersects(_cluster.geometry))]
            gdf.drop(lines.index, inplace=True)
            if len(lines) == 0:
                continue

            # Step 2: Eliminate the unexpected contours according to the "start contour value"
            brightThres = np.percentile(np.array([NTLArray[row][col] for row, col in _cluster.Pixels]), const.START_CONTOUR_VALUE_PERCENTILE)
            lines = lines[lines['Value'] >= brightThres]  # winnowing by value
            lines.sort_values(by=['Area'], ascending=False, inplace=True)
            contours = [ Contour(i, line['Value'], line['geometry'], line['Area']) for i, (_, line) in enumerate(lines.iterrows()) ] # <ID: :obj:`Contour`>

            # Step 3: Construct topology tree
            # Identify all childs for each contour, in the descending order of area.
            for i in range(len(contours)):
                contour = contours[i]
                for j in range(i + 1, len(contours)):
                    _contour = contours[j]
                    # If _contour has different parent compared to contour, then contour will not include _contour.
                    if _contour.Parent != contour.Parent:
                        continue
                    # If _contour has one point located within contour, then it must be fully included by contour.
                    coords = list(_contour.geometry.coords)
                    if Point(coords[int(len(coords) / 2)]).within(Polygon(contour.geometry)):
                        contour.Childs.append(j)
                        # Update the _contour's parent and child
                        if _contour.Parent != -1:
                            contours[_contour.Parent].Childs.remove(j)
                        _contour.Parent = i

            # Step 4: Extract seed contours
            LCs = [contour for contour in contours if len(contour.Childs) == 0]  # Leaf contours
            _seedContours = [LC for LC in LCs if LC.Parent == -1 or contours[LC.Parent].Value < LC.Value]  # Seed contours, i.e., peak leaf contours

            # If no elementary contours, then drop this urban area
            if len(_seedContours) == 0:
               continue
            _cluster.CenterNumber = len(_seedContours)
            _cluster.ID = len(filteredClusters)
            filteredClusters.append(_cluster)

            # Step 5: Extract centers of gravity of seed contours
            _centers = []  # Current cluster
            for i, _sc in enumerate(_seedContours):
                _sc.CenterID = i + len(self.centers)

                lat, lng = self.GetCentroid(_sc.geometry, NTL)
                _centers.append( Center(i + len(self.centers), Point(lng, lat), lat, lng, _cluster.ID) )

            # Step 6: Identify the main center
            self.IdentifyMainCenter(_centers, contours, PopArray, PopTransform)

            self.centers.extend(_centers)
            self.seedContours.extend(_seedContours)

        self.clusters = filteredClusters


    def IdentifyBasecontour(self, alterBasecontours, PopArray, PopTransform):
        """
        Identify the base contour from alter basecontours.

        The contour with the highest population density among those whose area is greater than half of the largest area is selected.

        Args:
            alterBasecontours (list): List of contour objects with the basecontour value.
            PopArray (ndarray): The array representation of the Pop raster.
            PopTransform (affine.Affine): The affine matrix of the Pop raster.

        Returns:
            Contour: The identified basecontour.

        """
        if len(alterBasecontours) == 1:
            basecontour = alterBasecontours[0]
        else:
            maxArea = max([contour.Area for contour in alterBasecontours])
            potentialContours = [contour for contour in alterBasecontours if contour.Area > maxArea * const.BASECONTOUR_AREA_THRESHOLD]

            if len(potentialContours) == 1:
                basecontour = potentialContours[0]
            else:
                basecontour = max(potentialContours,
                key=lambda contour: zonal_stats(Polygon(contour.geometry), PopArray, affine=PopTransform, stats='sum')[0]['sum'] / contour.Area)
        return basecontour

    def IdentifyMainCenter(self, centers, contours, PopArray, PopTransform):
        """
        Identify the main center among centers.

        Args:
            centers (lisf of :obj:`Center`): A dictionary where keys are IDs and values are the ``Center`` objects of the urban area.
            contours (list of :obj:`Contour`): The contours within the urban area.
            PopArray (ndarray): The array representation of the Pop raster.
            PopTransform (affine.Affine): The affine matrix of the Pop raster.

        Returns:
            None. Set the ``IsMainCenter`` attribute of the identified center to be true.
        """

        if len(centers) == 1:
            centers[0].IsMainCenter = 1

        else:
            baseValue = min([_contour.Value for _contour in contours])
            alterBasecontours = [contour for contour in contours if contour.Value == baseValue]
            basecontour = self.IdentifyBasecontour(alterBasecontours, PopArray, PopTransform)

            while len(basecontour.Childs) != 0:
                baseValue += const.CONTOUR_INTERVAL
                alterBasecontours = [contours[i] for i in basecontour.Childs]
                basecontour = self.IdentifyBasecontour(alterBasecontours, PopArray, PopTransform)

            basecontour.IsMainCenter = 1
            _center = next((_center for _center in centers if _center.ID == basecontour.CenterID), None)
            _center.IsMainCenter = 1
    def AssignGeoname(self, geonamePath):
        """
        Assign geonames to urban areas and its centers.

        If the country/region don't have geonames dataset, then ignore this step.

        Args:
            geonamePath (str): The file pathe of the GeoNames of the country/region.

        Returns:
            Append the ``name`` attribute to global variables ``clusters`` and ``centers``.
        """

        if not os.path.exists(geonamePath):
            return

        names = gpd.read_file(geonamePath)  # point features having geoname
        names = names.to_crs("EPSG:4326")

        # Assign geoname to the urban areas
        for _cluster in self.clusters:
            results = sjoin(names, gpd.GeoDataFrame(geometry=[_cluster.geometry], crs="EPSG:4326"), how="inner", predicate="within")
            if results.empty:
                continue

            _cluster.Name = results.loc[results['population'].idxmax(), 'asciiname']

        # Assign geoname to the centers
        for _sc in self.seedContours:
            results = sjoin(names, gpd.GeoDataFrame(geometry=[Polygon(_sc.geometry)], crs="EPSG:4326"), how="inner", predicate="within")
            if results.empty:
                continue
            _center = next((_center for _center in self.centers if _center.ID == _sc.CenterID), None)
            _center.Name = results.loc[results['population'].idxmax(), 'asciiname']



    def Execute(self, Clusters, contourPath, geonamePath, NTL, NTLArray, PopArray, PopTransform):
        """
        Generate the centers within the urban areas.

        The centers are generated by the ``Localized Contour Tree Method`` and the main centers are identified by the
        ``Iterative Main Tree Selection Method``.

        Args:
            Clusters (list): The list of ``Cluster`` objects representing the identified urban areas.
            contourPath (str): The file path of the contours of the country/region.
            geonamePath (str): The file pathe of the GeoNames of the country/region.
            NTL (rasterio.DataReader): The NTL raster object.
            NTLArray (ndarray): The array representation of the NTL raster.
            PopArray (ndarray): The array representation of the Pop raster.
            PopTransform (affine.Affine): The affine matrix of the Pop raster.

        Returns:
            None. Construct the global variables ``clusters`` and ``seedContours``.
        """

        self.centers = []
        self.seedContours = []
        self.clusters = Clusters

        self.DetectCenter(contourPath, NTL, NTLArray, PopArray, PopTransform)
        self.AssignGeoname(geonamePath)
