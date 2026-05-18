# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2026-04

import os
import time
from multiprocessing import Pool
from functools import partial
import numpy as np
import rasterio
import geopandas as gpd
from tqdm import tqdm
from GHSL.centers_generator import CentersGenerator
from GHSL.contours_generator import ContoursGenerator
import Parameters.consts as const

"""
Module generating the economic centers for all urban areas globally under different contour parameters.

This module calls the functions defined in ``centers_generator.py'' to generate the economic centers.
"""


def export(generator, area):
    """
    Export the identified centers and their corresponding contours as shapefiles.

    Args:
        generator (CentersGenerator): The generator object containing identified centers and seed contours.
        area (int): The minimum contour area parameter value used in this run.
    """
    # Centers
    points = gpd.GeoDataFrame({
        'geometry': [_center.geometry for _center in generator.centers],
        **{col: [_center.__getattribute__(attr) for _center in generator.centers] for col, attr in
           [('ID', 'ID'), ('ID_UC_G0', 'ClusterID'), ('Latitude', 'Latitude'), ('Longitude', 'Longitude'),
            ('Main', 'IsMainCenter')]}}, crs="EPSG:4326")

    # Contours
    lines = gpd.GeoDataFrame({
        'geometry': [_sc.geometry for _sc in generator.seedContours],
        'ID_UC_G0': [_sc.ClusterID for _sc in generator.seedContours],
        **{col: [_sc.__getattribute__(col) for _sc in generator.seedContours] for col in
           ['CenterID', 'Area', 'POIDensity']}}, crs="EPSG:4326")

    points.to_file(fileItems[6].format(area=area))
    lines.to_file(fileItems[7].format(area=area))


def run_single(area, clustersPath, NTLPath, smoothNTLPath, POIPath, contourPath,
               mmap_path, mmap_shape, mmap_dtype, NTLTransform):
    """
    Run center identification and export for a single area parameter.

    Args:
        area (int): The minimum contour area parameter value.
        clustersPath (str): The file path of the urban area clusters.
        NTLPath (str): The file path of the NTL raster.
        smoothNTLPath (str): The file path of the smoothed NTL raster.
        POIPath (str): The file path of the POI data.
        contourPath (str): The file path of the contours.
        mmap_path (str): Path to the memory-mapped NTL array file.
        mmap_shape (tuple): Shape of the NTL array.
        mmap_dtype (str): Dtype of the NTL array.
        NTLTransform (affine.Affine): The transform matrix of the NTL raster.
    """
    t0 = time.time()

    clusters = gpd.read_file(clustersPath)
    NTLArray = np.memmap(mmap_path, dtype=mmap_dtype, mode='r', shape=mmap_shape)

    generator = CentersGenerator()
    generator.Execute(clusters, NTLPath, smoothNTLPath, POIPath, contourPath,
                      minContourArea=area, NTLArray=NTLArray, NTLTransform=NTLTransform, verbose=False)

    export(generator, area)
    return area, len(generator.centers), len(generator.seedContours), time.time() - t0

fileItems = [
    r'D:\Research\CityCenter\2020\Input\VNL_v2_npp_2015_global_vcmslcfg_c202102150000.average_masked.tif',  # 0: NTL
    r'D:\Research\CityCenter\2020\Input\smoothed_ntl.tif',          # 1: Smoothed NTL
    r'D:\Research\CityCenter\2020\Input\contours.gpkg',             # 2: Contours
    r'D:\Research\CityCenter\GHS-FUA\GHS_FUA_UCDB2015_GLOBE_R2019A_54009_1K_V1_0.gpkg',                    # 3: Clusters
    r'D:\Research\CityCenter\2020\Input\poi.gpkg',                  # 4: POI
    r'D:\Research\CityCenter\2020\Input\_ntl_shared.dat',           # 5: NTL memmap temp
    r'D:\Research\CityCenter\GHS-FUA\Center\center_a{area}.shp',       # 6: Center output
    r'D:\Research\CityCenter\GHS-FUA\Seed\seed_a{area}.shp',           # 7: Seed output
]


if __name__ == '__main__':
    parms = range(1, 11)

    # Generate smoothed NTL and contours from the global NTL raster if they don't exist
    if not os.path.exists(fileItems[1]):
        contours_gen = ContoursGenerator()
        contours_gen.Execute(fileItems[0], fileItems[1], fileItems[2])

    # Load NTL array once and write to a temporary memmap file for cross-process sharing
    with rasterio.open(fileItems[0]) as src:
        NTLArray = src.read(1)
        NTLTransform = src.transform

    mmap_array = np.memmap(fileItems[5], dtype=NTLArray.dtype, mode='w+', shape=NTLArray.shape)
    mmap_array[:] = NTLArray
    mmap_array.flush()
    NTLShape = NTLArray.shape
    NTLDtype = str(NTLArray.dtype)
    del NTLArray, mmap_array

    worker = partial(
        run_single,
        clustersPath=fileItems[3],
        NTLPath=fileItems[0],
        smoothNTLPath=fileItems[1],
        POIPath=fileItems[4],
        contourPath=fileItems[2],
        mmap_path=fileItems[5],
        mmap_shape=NTLShape,
        mmap_dtype=NTLDtype,
        NTLTransform=NTLTransform,
    )

    t_start = time.time()

    with Pool(processes=2) as pool:
        with tqdm(total=len(parms), desc="Overall") as pbar:
            for area, n_centers, n_contours, elapsed in pool.imap_unordered(worker, parms):
                pbar.set_postfix(last=f"a{area}({n_centers}c,{elapsed:.0f}s)")
                pbar.update(1)

    print(f"All tasks finished in {time.time() - t_start:.1f}s")
    os.remove(fileItems[5])
