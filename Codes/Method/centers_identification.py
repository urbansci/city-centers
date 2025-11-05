# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2024-12
import os.path

from .centers_generator import *

"""
Module generating the economic centers for each country/region under different contour parameters.

This module calls the functions defined in ``contours_generator.py'' and ``centers_generator.py'' to generate 
the economic centers.
"""

class CityCentersDatasetGenerator():
    """
    A class generating the economic centers for the specific country/region.

    This class integrates the class ``CentersGenerator`` to export the center files.

    Args:
        None.

    Attributes:
        country (int): The index of the country/region.
        centersGenerator (:obj:`CentersGenerator`): A class object to generator the centers.
        fileItems (list[str]): A list containing the paths of all input datasets and output files.
   """
    def __init__(self):

        self.country = None
        self.centersGenerator = CentersGenerator()
        self.fileItems = None

    def Export(self):
        """
        Export the identified centers and their corresponding contours as shapefiles.

        Args:
            None.

        Returns:
            None. Results are saved to shapefiles.
        """
        # Centers
        points = gpd.GeoDataFrame({
            'geometry': [_center.geometry for _center in self.centersGenerator.centers],
            **{col: [_center.__getattribute__(attr) for _center in self.centersGenerator.centers] for col, attr in
               [('ID', 'ID'), ('Cluster_ID', 'ClusterID'), ('Latitude', 'Latitude'),  ('Longitude', 'Longitude'),
                ('Is_MC', 'IsMainCenter'), ("Name", "Name")]}}, crs="EPSG:4326")

        #  Contours
        lines = gpd.GeoDataFrame({
            'geometry': [_sc.geometry for _sc in self.centersGenerator.seedContours],
            **{col: [_sc.__getattribute__(col) for _sc in self.centersGenerator.seedContours] for col in
               ['CenterID', 'Area', 'POIDensity']}}, crs="EPSG:4326")

        points.to_file(const.FILE_PATH + self.fileItems[5])
        lines.to_file(const.FILE_PATH + self.fileItems[6])

    def Execute(self, country):
        """
        Generator the centers dataset for the specific country/region.

        Args:
            country (int): The index of the country/region to be processed.

        Returns:
            None. Save the dataset, including centers and urban areas, to shapefiles.
        """
        self.country = country
        self.fileItems = [f'Data\\Bright\\Country\\Light{country}.tif',
                          f'Data\\Bright\\Gaussian\\Gaussian{country}.tif',
                          f'Data\\POI\\Country\\{country}.gpkg',
                          f'GHSL\\Contour\\Contour{country}.shp',
                          f'GHSL\\Cluster\\{country}.shp',
                          f'GHSL\\Center\\Center{country}_a{const.MINIMUM_CONTOUR_AREA}.shp',
                          f'GHSL\\Seed\\Seed{country}_a{const.MINIMUM_CONTOUR_AREA}.shp']
        if not os.path.exists(const.FILE_PATH + self.fileItems[4]):
            return

        clusters = gpd.read_file(const.FILE_PATH + self.fileItems[4])
        self.centersGenerator.Execute(clusters,
                                      const.FILE_PATH + self.fileItems[0],
                                      const.FILE_PATH + self.fileItems[1],
                                      const.FILE_PATH + self.fileItems[2],
                                      const.FILE_PATH + self.fileItems[3],)
        self.Export()


if __name__ == '__main__':
    parms = [1, 2, 3, 4, 5, 6, 7, 8] # potential minimum contour area values

    generator = CityCentersDatasetGenerator()
    for country in tqdm(range(200, 234), desc="Country (200 - 233)", position=0, leave=True):
        for area in tqdm(parms, desc=f"{country}", position=1, leave=False):
            const.MINIMUM_CONTOUR_AREA = area
            generator.Execute(country)

        # Clear Rtree
        const.contour_tree = None
        const.poi_tree = None