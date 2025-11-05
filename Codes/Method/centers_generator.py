# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2025-11

import os
from itertools import islice
import fiona
from shapely.geometry import shape, Polygon, Point
from collections import defaultdict
import rasterio
import rasterio.mask
import pandas as pd
from geopandas.tools import sjoin
from tqdm import tqdm
from rasterstats import zonal_stats
from shapely.strtree import STRtree
from .contours_generator import *

"""
Module generating centers within urban areas for each country/region.

This module identifies the economic centers and labels the main center within the urban area.
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
        CenterID (int): The ID of the center corresponding to this contour. Only used for seed contours (contours correspond
            to the economic centers).
        ClusterID (int): The ID of the city containing this contour.
        Value (int): The brightness value of the contour.
        Area (float): The enclosed area of the contour.
        Parent (int): The ID of the contour's parent. -1 for root contour.
        Childs (list): The list of ID of the childs of this contour, only direct-connected childs.
        IsMainCenter (int): Whether the corresponding center is the main center of the urban area, 0 for false, 1 for true.
            Only used for seed contours.
        POICount (int): The number of POI points locating within this contour.
        POIDensity (int): The number density of the POI points locating within this contour.
    """

    def __init__(self, geometry, ID, ClusterID,  Value,  Area=0):
        self.geometry = geometry
        self.ID = ID
        self.CenterID = -1
        self.ClusterID = ClusterID
        self.Value = Value
        self.Area = Area
        self.Parent = -1
        self.Childs = []
        self.IsMainCenter = 0
        self.POICount = 0
        self.POIDensity = 0

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

    Args:
        None.

    Attributes:
        contoursGenerator (:obj:`ContourGenerator`): A class object to generate the contours map.
        clusters (geopandas.GeoDataFrame): A object containing the input urban areas.
        seedContours (list[Contour]): A list of the identified seed contours.
        centers (list[Center]): A list of the identified centers.
    """

    def __init__(self):

        self.contoursGenerator = None
        self.clusters = None
        self.centers = None
        self.seedContours = None

    def GetBrightestPoint(self, geometry, NTL):
        """
        Extract the brightest point (maximum NTL) within the given geometry.

        Args:
            geometry (shapely.geometry): The geometry object of the seed contour.
            NTL (rasterio.DatasetReader): The NTL raster object.

        Returns:
            float: The latitude of the brightest point.
            float: The longitude of the brightest point.
        """

        out_image, out_transform = rasterio.mask.mask(NTL, [Polygon(geometry)], crop=True, all_touched=True)
        out_image = out_image[0]
        mask_valid = out_image != NTL.nodata

        max_val = out_image[mask_valid].max()
        rows, cols = np.where((out_image == max_val) & mask_valid)

        row, col = rows[0], cols[0]
        lon, lat = out_transform * (col + 0.5, row + 0.5)

        return lat, lon

    def DetectCenter(self, NTLPath, contourPath):
        """
        Detect the centers within this urban area .

        Args:
            NTLPath (str): The file path of the NTL data.
            contourPath (str): The file path of the contours of the country/region.

        Returns:
            None. Initialize the global variable ``centers'' and ``seedContours'' and filter the ``clusters''.
        """

        NTL, NTLArray, NTLTransform = self.ReadRaster(NTLPath)
        gdf = gpd.read_file(contourPath)

        # Construct spatial Rtree for contours
        if const.contour_tree is None:
            const.contour_tree = STRtree(gdf.geometry.values)
        tree = const.contour_tree

        self.clusters.sort_values(by='Area', ascending=False, inplace=True)
        for idx, _cluster in tqdm(self.clusters.iterrows(), desc=f"cluster", position=2, leave=False, total=self.clusters.shape[0]):
            # Step 1: Extract all contours within the cluster
            indices = tree.query(_cluster.geometry)
            lines = gdf.iloc[[idx for idx in indices if gdf.iloc[idx]['Area'] >= const.MINIMUM_CONTOUR_AREA and gdf.iloc[idx].geometry.intersects(_cluster.geometry)]]
            if len(lines) == 0:
                continue

            lines = lines.sort_values(by='Area', ascending=False)
            contours = [Contour(line['geometry'], i, idx, line['Value'], line['Area']) for i, (_, line) in
                        enumerate(lines.iterrows())]  # [ID: :obj:`Contour`]

            # Step 2: Construct topology tree
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

            # Step 3: Extract leaf contours corresponding to the peak terrain
            lcs = [contour for contour in contours if len(contour.Childs) == 0]  # Leaf contours
            pscs = [] # Potential Seed contours, i.e., peak leaf contours
            for lc in lcs:
                if lc.Parent == -1: # for isolated/root contours, determine whether the terrain is a peak or a valley
                    mean_val = zonal_stats(
                        Polygon(lc.geometry),
                        NTLArray,
                        affine = NTLTransform,
                        stats=['mean'],
                        all_touched=True,
                        nodata=None
                    )
                    if lc.Value < mean_val[0]['mean']:
                        pscs.append(lc)
                elif contours[lc.Parent].Value < lc.Value: # ensuring the corresponding terrain represents a peak rather than a valley
                    pscs.append(lc)

            # Step 4: Extract the brightest pixel of the seed contour
            _centers = []
            _seedContours = []
            i = 0
            for psc in pscs:
                lat, lng = self.GetBrightestPoint(psc.geometry, NTL)
                # Ensure the center locate within the city boundary.
                if Point(lng, lat).within(_cluster.geometry):
                    _centers.append(Center(i + len(self.centers), Point(lng, lat), lat, lng, idx))
                    _seedContours.append(psc)
                    psc.CenterID = i + len(self.centers)
                    i += 1

            # If no elementary contours, then drop this urban area
            if len(_seedContours) == 0:
                continue

            self.centers.extend(_centers)
            self.seedContours.extend(_seedContours)

    def IdentifyMainCenter(self, POIPath):
        """
        Identify the main center among centers of each urban area.

        Args:
            POIPath (str): The file path of the POI Data.

        Returns:
            None. Set the ``IsMainCenter'' attribute of the identified main centers to be true.
        """

        # Calculate the count of poi points within each seed contour
        BATCH_SIZE = 100_000
        trees = []
        with fiona.open(POIPath, layer=0) as poi_src:
            batch_index = 0
            iterator = iter(poi_src)

            with tqdm(desc="Batch", position=2, leave=False) as pbar:
                while True:
                    # Read next batch
                    batch_feats = list(islice(iterator, BATCH_SIZE))
                    if not batch_feats:
                        break

                    # Build geometries
                    batch = [shape(feat['geometry']) for feat in batch_feats]

                    # Build STRtree spatial index for the current batch
                    if const.poi_tree is None:
                        tree = STRtree(batch)
                        trees.append(tree)
                    else:
                        tree = const.poi_tree[batch_index]

                    # Inner loop: update each contour
                    for contour in self.seedContours:
                        poly = Polygon(contour.geometry)
                        indices = tree.query(poly)
                        contour.POICount += sum(1 for idx in indices if poly.contains(batch[idx]))

                    # Update progress bar by batch
                    pbar.update(1)
                    batch_index += 1

        if const.poi_tree is None:
            const.poi_tree = trees

        # Group seed contours by urban area
        groups = defaultdict(list)
        for contour in self.seedContours:
            groups[contour.ClusterID].append(contour)
            contour.POIDensity = contour.POICount / contour.Area if contour.Area > 0 else 0

        # Identify the main center for each urban area
        for contours in groups.values():
            # Mark the contour with the highest POI density as the main center
            mainContour = max(contours, key=lambda x: x.POIDensity)
            mainContour.IsMainCenter = 1
            self.centers[mainContour.CenterID].IsMainCenter = 1


    def ReadRaster(self, path):
        """
        Read the raster specified by file path.

        Args:
            path (str): The file path of the raster.

        Returns:
            rasterio.DatasetReader: The rasterio object of the raster.
            ndarray: The array representation of the raster.
            affine.Affine: The transform matrix of the raster.
        """
        if not os.path.exists(path):
            print(f'{path} do not exist!')
            return

        raster = rasterio.open(path)
        array = raster.read(1)
        transform = raster.transform

        return raster, array, transform


    def Execute(self, clusters, NTLPath, smoothNTLPath, POIPath, contourPath):
        """
        Generate the centers within the urban areas.

        Args:
            clusters (geopandas.GeoDataFrame): The object representing the urban areas.
            NTLPath (str): The file path of the NTL of the country/region.
            smoothNTLPath (str): The file path of the smoothed NTL of the country/region.
            POIPath (str): The file path of the POI data of the country/region.
            contourPath (str): The file path of the contours of the country/region.

        Returns:
            None. Construct the global variables ``clusters`` and ``seedContours``.
        """

        if not os.path.exists(contourPath):
            self.contoursGenerator = ContoursGenerator()
            self.contoursGenerator.Execute(NTLPath, smoothNTLPath, contourPath)

        self.centers = []
        self.seedContours = []
        self.clusters = clusters

        self.DetectCenter(NTLPath, contourPath)
        self.IdentifyMainCenter(POIPath)