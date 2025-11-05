# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2025-11

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import Parameters.consts as const
import Parameters.plot_utils as utils
from GHSL.attribute_enhancement import fileItems
import pandas as pd
from scipy.stats import gaussian_kde
from geovoronoi import voronoi_regions_from_coords
from pyproj import Geod
from tqdm import tqdm
from shapely.geometry import LineString
from shapely.ops import split
from shapely.ops import unary_union

"""
Module for the generation of the ``Spatial distribution of centers within the city'' plot (Figure 4b, 4c and 4d).
"""

def world_scale_area(ax, gdf):
    """
    Plot kernal density estimation of city areas grouped by the number of centers (Figure 4b).

    Args:
        ax (plt.Axes): The axes to be plotted on.
        gdf (gpd.GeoDataFrame): The dataframe containing the center number and area of cities.
    Return:
        None.
    """
    center_nums = [1, 2, 3, 4, 5]
    colors = ['#52AECC', '#CFCED4','#AFA9F2', '#EB98F2', '#B669BF']

    for i, cnum in enumerate(center_nums):
        # Filter cities with specific center number
        subset = gdf[gdf['Center_Num'] == cnum]
        data = subset['Area']

        # Plot KDE
        sns.kdeplot(ax=ax, data=data, bw_adjust=0.5, color=colors[i], fill=True, alpha=0.3)

        # Plot vertical line for mean
        mean_val = data.mean()
        ax.axvline(mean_val, color=colors[i], linestyle=(-3, (6, 6)), linewidth=0.8)

    ax.set_xlim(0, 400)
    ax.set_xticks([0, 100, 200, 300, 400])
    ax.set_ylim(0, 0.05)
    ax.set_yticks([0, 0.01, 0.02, 0.03, 0.04, 0.05])

    ax.set_xlabel(rf'Area of the city (km²)', fontname='Arial')
    ax.set_ylabel('Kernel density estimation')

    ax.legend()
    utils.FormatAxis(ax)

def country_scale_area(ax, gdf):
    """
    Violin plot to show the area distribution of centers' voronoi polygons for multiple countries (Figure 4c).

    Args:
        ax (plt.Axes): The axes to be plotted on.
        gdf (gpd.GeoDataFrame): The dataframe containing the center data, with columns including country, city id, and area.

    Return:
        None.
    """

    group_sizes = gdf.groupby(['Country', 'ID_UC_G0'])['ID_UC_G0'].transform('size')
    gdf['Type'] = group_sizes.apply(lambda x: 'Mono' if x == 1 else 'Poly')

    colors = {'Mono': '#1DA0CC', 'Poly': '#C061CC'}
    violin_width = 25
    violin_gap = 0.25

    scatter_size = 15
    lw = 0.5

    for i, country in enumerate(countries):
        for j, t in enumerate(['Mono', 'Poly']):
            vals = gdf[(gdf['Country'] == country) & (gdf['Type'] == t)]['Area'].values

            kde = gaussian_kde(vals, bw_method=0.4)
            y = np.linspace(0, 200, 200)
            density = kde(y)
            density = density * (violin_width / 2)

            if t == 'Mono':
                pos = i - violin_gap/2
                ax.fill_betweenx(y, pos - density, pos, facecolor=colors[t], alpha=0.6, linewidth=0)
            else:
                pos = i + violin_gap/2
                ax.fill_betweenx(y, pos, pos + density, facecolor=colors[t], alpha=0.6, linewidth=0)

            box_pos = i + (j - 0.5) * (violin_gap - 0.1)
            mean_val = np.mean(vals)
            std = np.std(vals)

            ax.plot(
                [box_pos, box_pos],
                [mean_val-std, mean_val+std],
                color='black',
                linewidth=lw,
                solid_capstyle='round',
                zorder=2
            )
            ax.scatter(
                box_pos,
                mean_val,
                marker='s',
                s=scatter_size,
                color=colors[t],
                edgecolor='none',
                lw=0,
                zorder=3
            )
            if t == 'Poly':
                print(country, mean_val)

    ax.legend().remove()

    ax.set_xticks(ticks=np.arange(len(countries)), labels=[const.NAMES[country] for country in countries])

    ax.set_yticks([0, 50, 100, 150])
    ax.set_xlabel(None)
    ax.set_ylabel(rf"Area of the center's voronoi polygon (km²)", fontname="Arial")

    ax.set_ylim(0, 150)
    utils.FormatAxis(ax, positions=[])
    ax.spines['bottom'].set_position(('outward', 8))
    ax.tick_params(axis='x', which='major', pad=3)

def city_scale_area(ax, gdf):
    """
    Violin plot to show the area distribution of centers' voronoi polygons for multiple cities (Figure 4d).

    Args:
        ax (plt.Axes): The axes to be plotted on.
        gdf (gpd.GeoDataFrame): The dataframe containing the center data, with columns including country, city id, and area.

    Return:
        None.
    """

    city_dict = {232:[10756, 10394], 150:[2292, 2824], 144:[1781, 1915], 2:[17, 930], 219:[11386, 11287], }
    city_name = {17:'Los Angeles', 930:'New York', 10394:'Guangzhou', 10756:'Shanghai', 2292:'Essen', 2824:'Berlin',
                 1915:'London', 1781:'Manchester', 11386:'Tokyo', 11287:'Osaka' }

    group_sizes = gdf.groupby(['Country', 'ID_UC_G0'])['ID_UC_G0'].transform('size')
    gdf['size'] = group_sizes

    gap = 0.4
    x_offset = 60
    y_offset = 0.04
    lw = 0.5
    violin_width = 10
    color = '#AF58BA'
    means = []

    for i, country in enumerate(city_dict.keys()):
        for j, city in enumerate(city_dict[country]):
            values = gdf[(gdf['Country'] == country) & (gdf['ID_UC_G0'] == city)]['Area'].values

            # KDE
            kde = gaussian_kde(values, bw_method=0.3)
            x = np.linspace(0, 150, 200)
            density = kde(x) * violin_width

            y_pos = i + j * gap
            ax.fill_between(
                x + x_offset,
                y_pos,
                y_pos + density,
                facecolor='#C061CC',
                alpha=0.6,
                edgecolor='none',
                lw=0,
                zorder=5
            )

            mean_val = np.mean(values)
            means.append(mean_val)
            y_center = y_pos - y_offset
            ax.scatter(
                mean_val + x_offset,
                y_center,
                marker='s',
                s=10,
                color=color,
                edgecolor='none',
                lw=0,
                zorder=10
            )
            std = np.std(values)
            ax.plot([x_offset + mean_val - std, x_offset + mean_val + std], [y_center, y_center], color='black', lw=lw, solid_capstyle='round')

            # baseline
            ax.hlines(y=y_pos, xmin=x_offset, xmax=x_offset+150, color=color, lw=0.2, alpha=0.6)

            ax.text(x_offset-10, y_pos, city_name[city], va='center', ha='right', color='black')
            # print(city_name[city], mean_val, std)

    for i, country in enumerate(city_dict.keys()):
        ax.text(x_offset-60, i + gap/2,const.NAMES[country], va='center', ha='center', color='black')
    ax.axvline(np.mean(means)+x_offset, color=color, lw=0.8, linestyle=(-3, (6, 6)), zorder=1)

    ax.set_xlim(0, 150+x_offset)
    ax.set_ylim(-0.2, 4.8)
    ax.set_xticks([x_offset, 50+x_offset, 100+x_offset, 150+x_offset], labels=[0,50,100,150])
    ax.set_yticks([])
    ax.set_xlabel(rf"Area of the center's voronoi polygon (km²)", fontname="Arial")
    ax.set_ylabel(None)

    utils.FormatAxis(ax, positions=[])

def generate_voronois():
    """
    Generate the voronoi partitions of the centers within each city.
    """

    def points_to_voronoi(points_gdf, clip_poly, idx, country):
        """
        Compute Voronoi polygons of centers within the urban area.

        Args:
            points_gdf (gpd.GeoDataFrame): The GeoDataFrame containing point geometries representing the centers.
            clip_poly (shapely.geometry.Polygon): The boundary of the urban area.
            idx (int): The unique identifier for the urban area for which the Voronoi polygons are computed.
            country (str): The unique identifier of the country to which the urban area belongs.

        Returns:
            gpd.GeoDataFrame: A GeoDataFrame containing the Voronoi polygons for the centers within the specified
                               urban area.
        """
        if len(points_gdf) == 1: # one center
            return gpd.GeoDataFrame({
                'ID_UC_G0': [idx],
                'Country': [country],
                'Main': [points_gdf.iloc[0]['Main']],
                'geometry': [clip_poly]
            }, crs=points_gdf.crs)

        elif len(points_gdf) == 2: # two centers
            coords = np.array(list(zip(points_gdf.geometry.x, points_gdf.geometry.y)))
            mid = coords.mean(axis=0)
            vec = coords[1] - coords[0]

            perp_vec = np.array([-vec[1], vec[0]])
            factor = 1e5
            p1 = mid - perp_vec * factor
            p2 = mid + perp_vec * factor
            split_line = LineString([p1, p2])

            split_result = split(clip_poly, split_line)
            split_polys = [p for p in split_result.geoms if not p.is_empty]

            # Merge to two multipolygons corresponding to two sides
            poly_left = []
            poly_right = []
            def side(poly, line_point, line_vec):
                # Determine the position of the polygon relative to the split line (left or right).
                c = np.array([poly.centroid.x, poly.centroid.y])
                v = c - line_point
                cross = line_vec[0] * v[1] - line_vec[1] * v[0]
                return cross

            for poly in split_polys:
                if side(poly, mid, perp_vec) > 0:
                    poly_left.append(poly)
                else:
                    poly_right.append(poly)

            voronois = [unary_union(poly_left), unary_union(poly_right)]

            poly_to_main = []
            if voronois[0].contains(points_gdf.geometry.iloc[0]):
                poly_to_main = [points_gdf.loc[0, 'Main'], points_gdf.loc[1, 'Main']]
            else:
                poly_to_main = [points_gdf.loc[1, 'Main'], points_gdf.loc[0, 'Main']]

            return gpd.GeoDataFrame({
                'ID_UC_G0': [idx] * 2,
                'Country': [country] * 2,
                'Main': poly_to_main,
                'geometry': voronois
            }, crs=points_gdf.crs)

        else: # more than two centers
            coords = np.vstack(points_gdf.geometry.apply(lambda p: (p.x, p.y)).values)
            region_polys, idx_map = voronoi_regions_from_coords(coords, clip_poly)
            out = gpd.GeoDataFrame({'ID_UC_G0': [idx] * len(region_polys),
                                    'Country': [country] * len(region_polys),
                                    'Main': [points_gdf.iloc[idx_map[i][0]]['Main'] for i in idx_map.keys()]},
                                   geometry=list(region_polys.values()), crs=points_gdf.crs)
            return out

    # Read data
    polys = gpd.read_file(const.FILE_PATH + fileItems[1])
    pts = gpd.read_file(const.FILE_PATH + fileItems[2])
    polys = polys[polys['Country'].isin(countries)]

    # To web mercator
    pts = pts.to_crs(epsg=3857)
    polys = polys.to_crs(epsg=3857)

    results = []

    for poly_idx, poly_row in tqdm(polys.iterrows(), total=len(polys)):
        poly_geom = poly_row.geometry
        idx = poly_row['ID_UC_G0']
        country = poly_row['Country']
        pts_in = pts[(pts['ID_UC_G0'] == idx) & (pts['Country'] == country)]
        vor_gdf = points_to_voronoi(pts_in.reset_index(drop=True), poly_geom.buffer(0), idx, country)

        results.append(vor_gdf)

    merged = gpd.GeoDataFrame(pd.concat(results, ignore_index=True), crs=polys.crs)
    merged = merged.to_crs(epsg=4326)
    geod = Geod(ellps='WGS84')
    merged['Area'] = merged.geometry.apply(lambda geom: abs(geod.geometry_area_perimeter(geom)[0]) / 1e6)
    merged.to_file(const.FILE_PATH + fileItems[3])

fileItems = [r'GHSL\Global\globalClusters.shp', r'GHSL\Global\splitClusters.shp', r'GHSL\Global\globalCenters_edited.shp',
             r'Statistics\Voronoi\voronoi.shp']
countries = [232, 150, 144, 2, 219]

if __name__ == '__main__':
    # generate_voronois()
    _, ax = plt.subplots(2, 2, figsize = (17/2.54, 15/2.54), dpi=600,  gridspec_kw={
        'width_ratios': [1, 1],
        'height_ratios': [3, 4]})

    world_scale_area(ax[0, 1], gpd.read_file(const.FILE_PATH + fileItems[0]))
    country_scale_area(ax[1, 0], gpd.read_file(const.FILE_PATH + fileItems[3]))
    city_scale_area(ax[1, 1], gpd.read_file(const.FILE_PATH + fileItems[3]))

    ax[0,0].axis('off')
    ax[0,0].annotate("a", xy=(-0.12, 1.05), xycoords="axes fraction", fontsize=utils.ORDER_SIZE, weight="bold")
    ax[0,1].annotate("b", xy=(-0.12, 1.05), xycoords="axes fraction", fontsize=utils.ORDER_SIZE, weight="bold")
    ax[1,0].annotate("c", xy=(-0.12, 1.05), xycoords="axes fraction", fontsize=utils.ORDER_SIZE, weight="bold")
    ax[1,1].annotate("d", xy=(-0.12, 1.05), xycoords="axes fraction", fontsize=utils.ORDER_SIZE, weight="bold")

    plt.subplots_adjust(left=0.06, bottom=0.08, right=0.98, top=0.96, wspace=0.1)
    plt.savefig(r'D:\Research\Graph\GHSL\panel\spatial_pattern.png')

