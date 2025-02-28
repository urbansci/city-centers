# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2024-12

import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import Parameters.consts as const
import Parameters.plot_utils as utils

"""
Module for the generation of `` The average coverage area of each center`` plot (Figure 2c).
"""

def PlotAreaPerCenter(ax):
    """
    Plot the average coverage area by each center relative to the center number of the urban area. Urban areas are grouped
    by the center number and those with centers more than 50 are grouped together as one.
    """

    clusters = gpd.read_file(const.FILE_PATH + fileItems[0])
    clusters['PerArea'] = clusters['Area'] / clusters['Center_Num']

    clusters_below_50 = clusters[(clusters['Center_Num'] < 50)]
    avg = clusters_below_50.groupby('Center_Num')['PerArea'].mean()
    std = clusters_below_50.groupby('Center_Num')['PerArea'].std()

    # Group urban areas with centers more than 50
    clusters_above_50 = clusters[clusters['Center_Num'] >= 50]
    clusters_above_50['Center_Num'] = 50
    avg_above_50 = clusters_above_50.groupby('Center_Num')['PerArea'].mean()
    std_above_50 = clusters_above_50.groupby('Center_Num')['PerArea'].std()
    avg = pd.concat([avg, avg_above_50])
    std = pd.concat([std, std_above_50])

    ax.errorbar(avg.index, avg, yerr=std, fmt='o', markersize=5, markerfacecolor='white', markeredgewidth=0.8,
                 capsize=3, capthick=0.8, color='#0D4A70', elinewidth=0.8)

    # Plot mean line
    Average = avg.mean()
    Variance = avg.var()
    print(Average, Variance)
    ax.hlines(Average, ax.get_xlim()[0], ax.get_xlim()[1], color='#C30D23', lw=0.8,  zorder=20)

    ax.set_xticks([0, 10, 20, 30, 40, 50], [0, 10, 20, 30, 40, "50+"])
    ax.set_xlim([0, 51])
    ax.set_ylim([0, 300])
    ax.set_xlabel('Center number', fontsize=utils.LABEL_SIZE+1, labelpad=0)
    ax.set_ylabel(f'Area per center (km$^{2}$)', fontsize=utils.LABEL_SIZE+1, labelpad=utils.LABEL_PAD)

    utils.FormatAxis(ax)
    ax.tick_params(axis='x', which='minor', length=0)
    ax.tick_params(axis='y', which='minor', length=0)

fileItems = ["Result/Shp/Global/GlobalAreas.shp"]

if __name__ == "__main__":
    plt.rcParams['font.family'] = 'Arial'
    f, ax = plt.subplots(1, 1, figsize=(17 / 2.54, 3 / 2.54), dpi=600)

    PlotAreaPerCenter(ax)

    plt.subplots_adjust(top=0.97, bottom=0.15, left=0.05, right=0.98, wspace=0.4, hspace=0.2)
    plt.savefig(r"./Graph/Area/Area Per Center.png")
