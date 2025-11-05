# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2025-11

import re
import geopandas as gpd
import pandas as pd
import os
from ..Parameters import consts as const

"""
Module integrating the original results under varying thresholds to produce the final centers dataset. 

Retain centers that are identified under at least half of the threshold settings (majority voting).
"""
def majority_voting(country, parms):
    """
    Obtained the voted results among different contour parameter values.

    Match the centers of a specific country/region across varying area thresholds, recording the occurrence count and
    source files of each center. Designate the one labeled most frequently as the main center to be the final main center.

    Args:
        country (int): The ID of the country.
        parms (list[int]): The parameter range.

    Returns:
        None. Save the matching result as a shapefile.
    """
    center_path = const.FILE_PATH + f'GHSL\\Center'
    contour_path = const.FILE_PATH + f'GHSL\\Seed'
    cluster_path = const.FILE_PATH + f'GHSL\\Cluster'

    # Concatenate the centers of different area thresholds
    points = []
    for area in parms:
        if not os.path.exists(os.path.join(center_path, f'Center{country}_a{area}.shp')):
            continue
        centers = gpd.read_file(os.path.join(center_path, f'Center{country}_a{area}.shp'))
        contours = gpd.read_file(os.path.join(contour_path, f'Seed{country}_a{area}.shp'))
        centers = centers.merge(
            contours[['CenterID', 'POIDensity']],
            left_on='ID',
            right_on='CenterID',
            how='left'
        )
        centers = centers.drop(columns=['CenterID'])
        centers["source"] = area

        if not centers.empty:
            points.append(centers)
    if len(points) == 0:
        return

    combined = gpd.GeoDataFrame(pd.concat(points, ignore_index=True), crs=points[0].crs)
    combined['x'] = combined.geometry.x.round(6)
    combined['y'] = combined.geometry.y.round(6)
    combined['pos'] = combined.apply(lambda row: (row['x'], row['y']), axis=1)

    # Append the 'ID_UC_G0' column
    clusters = gpd.read_file(os.path.join(cluster_path, f'{country}.shp'))
    combined = combined.merge(clusters['ID_UC_G0'], left_on='Cluster_ID', right_index=True, how='left')

    grouped = combined.groupby(['Cluster_ID', 'pos'])
    result_rows = []
    for (cluster_id, pos), group in grouped:
        sources = group['source'].unique()
        main_prob = group['Is_MC'].sum() / len(parms)
        avg_density = group[group['Is_MC'] == 1]['POIDensity'].mean() if main_prob > 0 else 0

        result_rows.append({
            "Cluster_ID": cluster_id,
            "ID_UC_G0": group.iloc[0]['ID_UC_G0'],
            "geometry": group.iloc[0].geometry,
            "Latitude": group.iloc[0]['Latitude'],
            "Longitude": group.iloc[0]['Longitude'],
            "Count": len(group),
            "Source": ", ".join(map(str, sorted(sources))),
            "Main_Prob": main_prob,
            "Avg_Dens": avg_density,
        })
    if len(result_rows) == 0:
        return

    result = gpd.GeoDataFrame(result_rows, crs=combined.crs)
    result['Main'] = 0

    slice = result[result['Count'] >= len(parms) / 2]
    slice = slice.sort_values(['Main_Prob', 'Avg_Dens'], ascending=[False, False])
    # If multiple centers are tied for the highest probability, compute the average POI density for each when serving as
    # the main center, and select the one with the highest value as the final main center.
    main_idx = slice.groupby('Cluster_ID').head(1).index
    result.loc[main_idx, 'Main'] = 1

    result = result[[
        "Cluster_ID",
        "ID_UC_G0",
        "Latitude",
        "Longitude",
        "Count",
        "Main",
        "geometry",
    ]]
    output_path = os.path.join(center_path, f"Matched\\{country}.shp")
    result.to_file(output_path, driver="ESRI Shapefile")

if __name__ == '__main__':
    parms = [1, 2, 3, 4, 5, 6, 7, 8]
    for country in range(234):
        majority_voting(country, parms)