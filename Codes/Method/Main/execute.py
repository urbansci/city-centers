# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2024-12

from Method.Main.centers_generator import *
from Method.Main.urban_areas_generator import *

"""
Module generating the city centers dataset for each country/region.

This module calls the functions defined in ``urban_areas_generator.py`` and ``centers_generator.py`` to generate 
the urban areas and centers in sequence. Note that this module should be executed after the ``pcca.py`` and
 ``contours_generation.py``.
"""

class CityCentersDatasetGenerator():
    """
    A class generating the centers dataset for the specific country/region.

    This class integrates two classes, the ``UrbanAreasGenerator``and ``CentersGenerator`` to generate the whole centers
    dataset for the specific country/region.

    Args:
        None.

    Attributes:
        country (int): The index of the country/region.
        urbanAreasGenerator (:obj:`UrbanAreasGenerator`): A class object to generate the urban areas.
        centersGenerator (:obj:`CentersGenerator`): A class object to generator the centers.
        fileItems (list[str]): A list containing the paths of all input datasets and output files.
   """
    def __init__(self):

        self.country = None
        self.urbanAreasGenerator = UrbanAreasGenerator()
        self.centersGenerator = CentersGenerator()
        self.fileItems = None

    def Export(self):
        """
        Export the global variables to shapefiles.

        Args:
            None.

        Returns:
            None. Results are saved to shapefiles.
        """
        # Urban areas
        polygons = gpd.GeoDataFrame({
            'geometry': [_cluster.geometry for _cluster in self.centersGenerator.clusters],
            **{col: [_cluster.__getattribute__(attr) for _cluster in self.centersGenerator.clusters] for col, attr in
               [('ID', 'ID'), ('Area', 'Area'), ('Bright', 'Brightness'), ('Pop', 'Population'),
                ('Center_Num', 'CenterNumber'), ("Name", "Name") ]}}, crs="EPSG:4326")

        # Centers
        points = gpd.GeoDataFrame({
            'geometry': [_center.geometry for _center in self.centersGenerator.centers],
            **{col: [_center.__getattribute__(attr) for _center in self.centersGenerator.centers] for col, attr in
               [('ID', 'ID'), ( 'Cluster_ID', 'ClusterID'), ('Latitude', 'Latitude'),  ('Longitude', 'Longitude'),
                ('Is_MC', 'IsMainCenter'), ("Name", "Name")]}}, crs="EPSG:4326")

        #  Elementary contours
        lines = gpd.GeoDataFrame({
            'geometry': [_sc.geometry for _sc in self.centersGenerator.seedContours],
            **{col: [_sc.__getattribute__(col) for _sc in self.centersGenerator.seedContours] for col in
               ['ID', 'Area']}}, crs="EPSG:4326")

        polygons.to_file(const.FILE_PATH + self.fileItems[4])
        points.to_file(const.FILE_PATH + self.fileItems[5])
        lines.to_file(const.FILE_PATH + self.fileItems[6])


    def ReadRaster(self, path):
        """
        Read the raster specified by file path.

        Args:
            path (str): The file path of the raster.

        Returns:
            raster (rasterio.DatasetReader): The rasterio object of the raster.
            array (ndarray): The array representation of the raster.
            transform (affine.Affine): The transform matrix of the raster.
        """
        if not os.path.exists(path):
            print(f'{path} do not exist!')
            return

        raster = rasterio.open(path)
        array = raster.read(1)
        transform = raster.transform

        return raster, array, transform

    def Execute(self, country):
        """
        Generator the centers dataset for the specific country/region.

        Args:
            country (int): The index of the country/region to be processed.

        Returns:
            None. Save the dataset, including centers and urban areas, to shapefiles.
        """
        print(f"Country : {country}")
        self.country = country
        self.fileItems = [f'Data\\Bright\\Country\\Light{country}.tif',
                          f'Data\\WorldPop\\Country\\Pop{country}.tif',
                          f'Data\\GeoNames\\{country}.shp',
                          f'Data\\Bright\\Contour\\Contour{country}.shp',
                          f'Result\\Shp\\Cluster\\{country}.shp',
                          f'Result\\Shp\\Center\\{country}.shp',
                          f'Result\\Shp\\Contour\\{country}.shp']

        NTL, NTLArray, NTLTransform = self.ReadRaster(const.FILE_PATH + self.fileItems[0])
        Pop, PopArray, PopTransform = self.ReadRaster(const.FILE_PATH + self.fileItems[1])

        clusters = self.urbanAreasGenerator.Execute(const.OPTIMAL_THRESHOLD[self.country] if const.OPTIMAL_THRESHOLD[self.country] != None else 0.5, NTLArray, NTLTransform, PopArray,
                                               PopTransform)
        self.centersGenerator.Execute(clusters, const.FILE_PATH + self.fileItems[3], const.FILE_PATH + self.fileItems[2], NTL,
                                 NTLArray, PopArray, PopTransform)
        self.Export()


if __name__ == '__main__':

    generator = CityCentersDatasetGenerator()
    for country in range(128,129):
        generator.Execute(country)