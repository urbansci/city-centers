# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2026-04

"""
Population-weighted average distance analysis (Figure 5 & Supplementary Figure 12).

Figure 5: Average distance to nearest center vs main center across center numbers.
Supplementary Figure 12: Robustness checks across population datasets and
    minimum contour area thresholds.
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from pyproj import Geod
from rasterio.features import geometry_mask, geometry_window
from scipy.optimize import curve_fit
from tqdm import tqdm

import Parameters.consts as const
import Parameters.plot_utils as utils

# ── Constants ─────────────────────────────────────────────────────────────────

CENTERS_PATH = const.FILE_PATH + r'2020\Center\center_a4.shp'
CLUSTERS_PATH = const.FILE_PATH + r'2020\Input\GHSL_UCDB_MTUC_2020_GLOBE_R2024_WGS84.shp'
POP_RASTER = const.FILE_PATH + r'Data\Pop\WorldPop\ppp_2020_1km_Aggregated.tif'
OUTPUT_FIG5 = r'D:\Research\Graph\GHSL\panel\accessibility.png'
OUTPUT_SUPP = r'D:\Research\Graph\GHSL\panel\extended\accessibility robustness.png'

POP_DATASETS = {
    'LandScan': (r'Data\Pop\LandScan\landscan-global-2020.tif', 'o', '#6a994e'),
    'GPW':      (r'Data\Pop\GPW\ciesen_nasa_gpw_v4_population_count_2020.tif', '^', '#1DA0CC'),
    'GHS-POP':  (r'Data\Pop\GHS-POP\GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif', 'v', '#e6b516'),
    'WorldPop': (r'Data\Pop\WorldPop\ppp_2020_1km_Aggregated.tif', 's', '#be52cc'),
}

THRESHOLDS = {
    '2': ('o', '#6a994e'),
    '3': ('^', '#1DA0CC'),
    '4': ('s', '#be52cc'),
    '5': ('v', '#e6b516'),
}

COLORS_MAIN = ['#CF82D9', '#1DA0CC']


# ── Core computation ──────────────────────────────────────────────────────────

def _average_distance_of_city(pop_raster, poly_geom, points_geoms):
    """Population-weighted average distance from raster pixels to nearest center (km)."""
    site_coords = np.array([(p.x, p.y) for p in points_geoms])

    window = geometry_window(pop_raster, [poly_geom])
    vals = pop_raster.read(1, window=window, masked=True)
    mask = geometry_mask([poly_geom], transform=pop_raster.window_transform(window), out_shape=vals.shape)
    vals.mask |= mask
    rows, cols = np.where(~vals.mask)

    xs, ys = pop_raster.xy(rows + window.row_off, cols + window.col_off)
    pixel_coords = np.column_stack([xs, ys])
    pixel_vals = vals[~vals.mask]
    if np.sum(pixel_vals) == 0 or len(pixel_vals) == 0 or len(site_coords) == 0:
        return None

    geod = Geod(ellps="WGS84")
    if len(site_coords) == 1:
        lon, lat = site_coords[0]
        dist = geod.inv(
            np.full(len(pixel_coords), lon),
            np.full(len(pixel_coords), lat),
            pixel_coords[:, 0],
            pixel_coords[:, 1]
        )[2]
    else:
        dist = np.min([
            geod.inv(
                np.full(len(pixel_coords), lon),
                np.full(len(pixel_coords), lat),
                pixel_coords[:, 0],
                pixel_coords[:, 1]
            )[2]
            for lon, lat in site_coords
        ], axis=0)

    return np.sum(pixel_vals * dist) / (1e3 * np.sum(pixel_vals))


def _compute_distances(centers_path, clusters_path, pop_path):
    """Compute polycentric and monocentric distances for all cities."""
    points = gpd.read_file(centers_path).to_crs(4326)
    polys = gpd.read_file(clusters_path).to_crs(4326)
    raster = rasterio.open(pop_path)

    center_counts = points.groupby('ID_UC_G0').size().reset_index(name='Center_Num')
    polys = polys.merge(center_counts, on='ID_UC_G0', how='left')
    polys['Center_Num'] = polys['Center_Num'].fillna(0).astype(int)
    polys = polys[polys['Center_Num'] > 0]

    results_all, results_main = [], []
    for _, poly in tqdm(polys.iterrows(), total=len(polys)):
        sites = points[points['ID_UC_G0'] == poly['ID_UC_G0']]
        main_candidates = sites[sites['Main'] == 1]
        if len(main_candidates) == 0:
            results_all.append(np.nan)
            results_main.append(np.nan)
            continue
        results_all.append(_average_distance_of_city(raster, poly.geometry, sites.geometry))
        results_main.append(_average_distance_of_city(raster, poly.geometry, [main_candidates.iloc[0].geometry]))

    polys['dist_nearest'] = results_all
    polys['dist_main'] = results_main
    return polys


def _compute_distances_single(centers_path, clusters_path, pop_path):
    """Compute nearest-center distance for all cities (single distance column)."""
    points = gpd.read_file(centers_path)
    polys = gpd.read_file(clusters_path)
    raster = rasterio.open(pop_path)

    rows = []
    for city_id, sites in points.groupby('ID_UC_G0'):
        poly = polys[polys['ID_UC_G0'] == city_id].iloc[0]
        rows.append({
            'Center_Num': len(sites),
            'dist': _average_distance_of_city(raster, poly.geometry, sites.geometry),
        })
    return pd.DataFrame(rows)


def _bin_by_center_num(df, col):
    """Bin cities by center number (1-19 individual, 20+ grouped), return grouped series."""
    bins = list(range(1, 21)) + [np.inf]
    labels = list(range(1, 21))
    df['bin'] = pd.cut(df['Center_Num'], bins=bins, labels=labels, right=False)
    return df.groupby('bin', observed=False)[col]


# ── Panel plot functions ──────────────────────────────────────────────────────

def _plot_actuality_vs_mono(ax):
    """Fig. 5: polycentric vs monocentric average distance."""
    polys = _compute_distances(CENTERS_PATH, CLUSTERS_PATH, POP_RASTER).dropna(subset='dist_nearest')

    dist = _bin_by_center_num(polys, 'dist_nearest')
    dist_main = _bin_by_center_num(polys, 'dist_main')

    x_vals = dist.mean().index.astype(float)

    _, caps1, _ = ax.errorbar(x_vals, dist.mean().values, yerr=dist.std().values,
                              fmt='s', markersize=4, markerfacecolor='white', markeredgecolor=COLORS_MAIN[0],
                              markeredgewidth=0.4, color=COLORS_MAIN[0], capsize=1.5, elinewidth=0.3, zorder=15)
    for cap in caps1:
        cap.set_markeredgewidth(0.3)

    _, caps2, _ = ax.errorbar(x_vals, dist_main.mean().values, yerr=dist_main.std().values,
                              fmt='^', markersize=4.5, markerfacecolor='white', markeredgecolor=COLORS_MAIN[1],
                              markeredgewidth=0.4, color=COLORS_MAIN[1], capsize=1.5, elinewidth=0.3, zorder=10)
    for cap in caps2:
        cap.set_markeredgewidth(0.3)

    print(f"Average distance to nearest center: {dist.mean().mean()}")
    print(f"Average distance to main center: {dist_main.mean().mean()}")

    ax.axhline(dist.mean().mean(), color=COLORS_MAIN[0], lw=0.6, zorder=5)

    def power_law(x, a, b):
        return a * x ** b

    x = dist_main.mean().index.astype(float)
    y = dist_main.mean().values
    params, cov = curve_fit(power_law, x, y)
    a, b = params
    perr = np.sqrt(np.diag(cov))
    y_pred = power_law(x, a, b)
    r2 = 1 - np.sum((y - y_pred) ** 2) / np.sum((y - np.mean(y)) ** 2)

    print(f"a = {a:.4f} +/- {perr[0]:.4f}")
    print(f"b = {b:.4f} +/- {perr[1]:.4f}")
    print(f"R2 = {r2:.4f}")

    x_fit = np.linspace(1, 21, 200)
    ax.plot(x_fit, power_law(x_fit, a, b), color=COLORS_MAIN[1], lw=0.6, zorder=5)

    ax.set_xlabel('Center number', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_ylabel('Population-weighted average\ndistance to the center (km)', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_xlim(0, 21)
    ax.set_xticks([1, 5, 10, 15, 20], [1, 5, 10, 15, '20+'])
    ax.set_ylim(0, 35)
    ax.set_yticks([0, 5, 10, 15, 20, 25, 30, 35])
    ax.minorticks_off()
    utils.FormatAxis(ax)


def _plot_robustness_series(ax, thres, pop_path, marker, color):
    """Single robustness series: average distance for a given threshold and pop dataset."""
    centers_path = const.FILE_PATH + rf'2020\Center\center_a{thres}.shp'
    df = _compute_distances_single(centers_path, CLUSTERS_PATH, pop_path).dropna(subset='dist')

    dist = _bin_by_center_num(df, 'dist')
    ax.errorbar(dist.mean().index, dist.mean().values, fmt=marker, markersize=6,
                markerfacecolor='none', markeredgecolor=color, color='black',
                markeredgewidth=0.5, capsize=2, capthick=0.25, elinewidth=0.25, zorder=10)

    mean_val = dist.mean().mean()
    std_val = dist.mean().std()
    ax.hlines(mean_val, *ax.get_xlim(), color=color, lw=0.5, zorder=5)

    print(f"Threshold: {thres}; pop: {pop_path}")
    print(f"Average: {mean_val:.4f}; STD: {std_val:.4f}")
    return mean_val, std_val


def _draw_legend(ax, items, stats):
    """Draw aligned custom legend with marker, name, and stats."""
    x_marker, x_name, x_stats = 0.5, 0.55, 0.72
    y_start, y_step = 0.96, 0.08
    fs = utils.LEGEND_SIZE

    for i, (name, display, marker, color) in enumerate(items):
        mean_val, std_val = stats[name]
        y = y_start - i * y_step
        ax.plot(x_marker, y, marker=marker, color=color, markerfacecolor='none', markeredgecolor=color,
                markeredgewidth=0.5, markersize=4, transform=ax.transAxes, clip_on=False)
        ax.text(x_name, y, display, fontsize=fs, va='center', ha='left', transform=ax.transAxes)
        ax.text(x_stats, y, f'({mean_val:.2f} $\\pm$ {std_val:.2f} km)', fontsize=fs, va='center', ha='left', transform=ax.transAxes)


def _plot_different_pop(ax):
    """Supp. Fig. 12a: robustness across population datasets."""
    stats = {}
    for name, (path, marker, color) in POP_DATASETS.items():
        stats[name] = _plot_robustness_series(ax, 4, const.FILE_PATH + path, marker, color)

    ax.set_xlabel('Center number', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_ylabel('Population-weighted average distance (km)', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_xlim(0, 21)
    ax.set_xticks([1, 5, 10, 15, 20], [1, 5, 10, 15, '20+'])
    ax.set_ylim(0, 10)
    ax.minorticks_off()

    items = [(name, name, POP_DATASETS[name][1], POP_DATASETS[name][2]) for name in POP_DATASETS]
    _draw_legend(ax, items, stats)
    utils.FormatAxis(ax)


def _plot_different_thres(ax):
    """Supp. Fig. 12b: robustness across minimum contour area thresholds."""
    stats = {}
    for thres, (marker, color) in THRESHOLDS.items():
        stats[thres] = _plot_robustness_series(ax, int(thres), POP_RASTER, marker, color)

    ax.set_xlabel('Center number', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_xlim(0, 21)
    ax.set_xticks([1, 5, 10, 15, 20], [1, 5, 10, 15, '20+'])
    ax.set_ylim(0, 10)
    ax.minorticks_off()

    items = [(t, f'{t} km$^2$', THRESHOLDS[t][0], THRESHOLDS[t][1]) for t in THRESHOLDS]
    _draw_legend(ax, items, stats)
    utils.FormatAxis(ax)


# ── Figure assembly ───────────────────────────────────────────────────────────

def figure5():
    """Compose and save Figure 5: polycentric vs monocentric average distance."""
    fig, axes = plt.subplots(1, 2, figsize=(18 / 2.54, 6 / 2.54), dpi=600)
    axes[0].axis('off')
    _plot_actuality_vs_mono(axes[1])

    axes[0].text(-0.1, 1.05, 'a', transform=axes[0].transAxes, fontsize=utils.ORDER_SIZE, fontweight='bold')
    axes[1].text(-0.1, 1.05, 'b', transform=axes[1].transAxes, fontsize=utils.ORDER_SIZE, fontweight='bold')

    plt.subplots_adjust(top=0.9, bottom=0.1, left=0.05, right=0.98, wspace=0.1)
    fig.savefig(OUTPUT_FIG5, dpi=600, facecolor='white')
    plt.close(fig)
    print(f"Saved: {OUTPUT_FIG5}")


def supplementary_figure12():
    """Compose and save Supplementary Figure 12: robustness checks."""
    fig, axes = plt.subplots(1, 2, figsize=(18 / 2.54, 6 / 2.54), dpi=600)
    _plot_different_pop(axes[0])
    _plot_different_thres(axes[1])

    axes[0].text(-0.1, 1.05, 'a', transform=axes[0].transAxes, fontsize=utils.ORDER_SIZE, fontweight='bold')
    axes[1].text(-0.1, 1.05, 'b', transform=axes[1].transAxes, fontsize=utils.ORDER_SIZE, fontweight='bold')

    plt.subplots_adjust(top=0.9, bottom=0.1, left=0.05, right=0.98, wspace=0.1)
    fig.savefig(OUTPUT_SUPP, dpi=600, facecolor='white')
    plt.close(fig)
    print(f"Saved: {OUTPUT_SUPP}")


if __name__ == '__main__':
    supplementary_figure12()
