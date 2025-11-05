# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2025-11

import glob
import pandas as pd
from tqdm import tqdm
import numpy as np
import geopandas as gpd
from shapely.geometry import Point, Polygon
from shapely.strtree import STRtree
import fiona
from fiona.env import Env
from geopy.distance import geodesic
import Parameters.consts as const

"""
Module providing additional attributes to the identified centers, including the geoname based on GeoNames and the functional
category based on Foursquare OS Places.

For each center, generate a 1-km radius buffer, classify its functional category based on the categorydistribution of 
poi points within it, and assign its geoname based on the GeoNames points within it.
"""

def generate_buffers():
    """
    Generate a 1km-radius buffer for each center for subsequent processing.
    """

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

    centers = gpd.read_file(const.FILE_PATH + fileItems[0])
    buffers = gpd.GeoDataFrame({'Center_ID': centers.index.tolist()},
                               geometry=[GeodesicBuffer(center.geometry, 1) for _, center in centers.iterrows()],
                               crs="EPSG:4326")
    del centers
    return buffers


def classify_based_on_buffers(buffers):
    """
    Classify the centers' functional categories.
    """

    def category_mapping():
        """
        Get the reclassified mapping of POI points from the fsq_category_id to the center category.
        """

        cat_df = pd.read_csv(const.FILE_PATH + r'Data\POI\personalization-apis-movement-sdk-categories_.csv')
        cat_df.rename(columns=lambda x: x.strip(), inplace=True)
        category_mapping = pd.Series(cat_df['Center Type'].values, index=cat_df["Category ID"]).dropna()
        return category_mapping

    def query_poi(parquet_file, polygons):
        """
        Process for each parquet file, retrieve the poi points within each contour.
        """

        df = pd.read_parquet(parquet_file)
        if 'geometry' not in df.columns:
            df = df.dropna(subset=['latitude', 'longitude'])
            geometry = [Point(xy) for xy in zip(df['longitude'], df['latitude'])]
            poi_gdf = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:4326')
        else:
            poi_gdf = gpd.GeoDataFrame(df, geometry=gpd.GeoSeries.from_wkt(df['geometry']), crs='EPSG:4326')

        # Construct spatial rtree
        geom_list = poi_gdf.geometry.values
        rtree = STRtree(geom_list)

        match_pois = []
        for idx, row in polygons.iterrows():
            polygon = Polygon(row.geometry)
            candidates = rtree.query(polygon)

            for pt in candidates:
                if polygon.contains(geom_list[pt]):
                    poi = poi_gdf.iloc[pt]
                    match_pois.append({
                        'Center_ID': row['Center_ID'],
                        'fsq_category_ids': list(poi['fsq_category_ids'])
                        if isinstance(poi["fsq_category_ids"], (list, tuple, np.ndarray)) and len(poi['fsq_category_ids']) > 0
                        else None
                    })

        return pd.DataFrame(match_pois)

    def classify_category(pois_df):
        """
        Classify the seed contours based on the category distribution of poi points within it.

        First, if the distribution entropy is higher than 1.5, label the contour as a 'Mixed' type; otherwise, label it
        as the most common type of poi points. For cases where more than one category tied for the first place,
        label it as a 'Mixed' type.
        """

        def cal_entropy(counts):
            probs = counts / counts.sum()
            return -np.sum(probs * np.log2(probs + 1e-10))

        results = []
        for center_id, group in pois_df.groupby('Center_ID'):

            all_categories = pd.Series([cat for cats in group['category'] for cat in cats])

            category_counts = all_categories.value_counts()
            entropy = cal_entropy(category_counts)

            if entropy > 1.5:
                category = 'Mixed'
            else:
                max_count = category_counts.max()
                top_categories = category_counts[category_counts == max_count]

                if len(top_categories) > 1:
                    category = 'Mixed'
                else:
                    category = top_categories.idxmax()

            results.append({
                'Center_ID': center_id,
                'Entropy': entropy,
                'Category': category
            })

        return pd.DataFrame(results)

    parquet_dir = const.FILE_PATH + 'Data\POI\parquet'
    parquet_files = sorted(glob.glob(f"{parquet_dir}\\places-*.parquet"))
    results = []

    mapping = category_mapping()
    # Chunked processing for each parquet file
    for file in tqdm(parquet_files, desc="Processing Parquet Files"):

        matched_df = query_poi(file, buffers)
        if matched_df.empty:
            continue

        matched_df['category'] = matched_df['fsq_category_ids'].apply(
            lambda ids: [mapping[i] for i in ids if i in mapping] if isinstance(ids, list) and any(
                i in mapping for i in ids) else None)
        matched_df = matched_df.dropna(subset=['category'])

        results.append(matched_df)
    pois = pd.concat(results, ignore_index=True)

    categories = classify_category(pois)
    categories.set_index(['Center_ID'], inplace=True)
    return categories

def names_assignment(buffers):
    """
    For each center, assign its name as the name of the most populous points within its buffer.
    """

    with Env():
        with fiona.open(const.FILE_PATH + fileItems[3], encoding='utf-8') as src:
            geonames = gpd.GeoDataFrame.from_features(src, crs=src.crs)
            geonames = geonames.to_crs(buffers.crs)

    joined = gpd.sjoin(buffers, geonames, how='left', predicate='contains').reset_index(drop=True)
    names = joined.groupby('Center_ID').apply(lambda x: x.loc[x['population'].idxmax(), 'asciiname']if x['population'].notna().any() else None)
    names.name = 'Name'

    return names

fileItems = [r'GHSL/Global/globalCenters_edited_v2.shp', r'GHSL\Global\buffers.shp', r'GHSL/Global/dataset.shp', r'Data/GeoNames/global.shp']
if __name__ == '__main__':
    centers = gpd.read_file(const.FILE_PATH + fileItems[0])

    buffers = generate_buffers()
    buffers.to_file(const.FILE_PATH + fileItems[1])

    categories = classify_based_on_buffers(buffers)
    names = names_assignment(buffers)

    centers['Entropy'] = -1
    centers['Category'] = 'No-POI'
    centers['Name'] = None
    centers.update(categories[['Entropy', 'Category']])
    centers.update(names)
    centers.to_file(const.FILE_PATH + fileItems[2])
