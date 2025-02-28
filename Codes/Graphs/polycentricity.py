# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2024-12

import geopandas as gpd
import numpy as np
import pandas as pd
from matplotlib.patches import Arc
import Parameters.consts as const
import Parameters.plot_utils as utils
import matplotlib.pyplot as plt

"""
Module for the generation of the ``Socioeconomic characteristics of monocentric and polycentric cities`` plot (Figure 3b).
"""

def Proportion():
    """
    Get the compositions of the indicators (number, area, population and brightness), broken down into four city categories:
    monocentric (1 center), low polycentric (2-5 centers), moderate polycentric (6-10 centers), and high polycentric (11+ centers) cities.

    Args:
        None.

    Returns:
        pd.Dataframe: The composition for each indicator, broken down by city categories.
        pd.Dataframe: The proportion shared by polycentric cities for each indicator.
    """

    clusters = gpd.read_file(const.FILE_PATH + fileItems[0])
    clusters = clusters[['Center_Num', 'Area', 'Pop', 'Bright']]

    bins = [0, 1, 5, 10, 999]
    labels = ['1', '2', '3', '4']
    clusters['Group'] = pd.cut(clusters['Center_Num'], bins=bins, labels=labels, right=True)
    group = clusters.groupby('Group').sum()
    group['Count'] = clusters.groupby('Group').size()
    group = group[const.COLUMNS]
    # Proportion shared by each category
    group_sum = group.sum(axis=0)
    proportion = group.div(group_sum, axis=1)

    # Proportion shared by polycentric
    prop_polycentric = proportion.drop('1')
    polycentric = prop_polycentric.sum(axis=0)

    print(group_sum, '\n', group, '\n', proportion, '\n', polycentric)
    return proportion, polycentric

def PlotPie(ax):
    """
    Plot pie to show composition of count, area, population and brightness.
    """

    ax1 = ax.inset_axes([0.0, 0.53, 0.5, 0.5])
    ax2 = ax.inset_axes([0.5, 0.53, 0.5, 0.5])
    ax3 = ax.inset_axes([0.0, 0.05, 0.5, 0.5])
    ax4 = ax.inset_axes([0.5, 0.05, 0.5, 0.5])
    axs = [ax1, ax2, ax3, ax4]

    width = 0.4

    for i in range(4):
        data = proportion[const.COLUMNS[i]]
        wedges, _ = axs[i].pie(data, labels=None, colors=utils.AREAS_COLORS, startangle=70, wedgeprops={'width': width})

        d = 0.2  # height of the line
        r_outer = 1

        angle_2 = (wedges[1].theta1 + wedges[1].theta2) / 2
        angle_4 = (wedges[3].theta1 + wedges[3].theta2) / 2
        x2, y2 = np.cos(np.radians(angle_2)) * r_outer, np.sin(np.radians(angle_2)) * r_outer
        x4, y4 = np.cos(np.radians(angle_4)) * r_outer, np.sin(np.radians(angle_4)) * r_outer
        x2_end, y2_end = x2 + np.cos(np.radians(angle_2)) * d, y2 + np.sin(np.radians(angle_2)) * d
        x4_end, y4_end = x4 + np.cos(np.radians(angle_4)) * d, y4 + np.sin(np.radians(angle_4)) * d

        axs[i].plot([x2, x2_end], [y2, y2_end], color='black', lw=utils.LINE_WIDTH)
        axs[i].plot([x4, x4_end], [y4, y4_end], color='black', lw=utils.LINE_WIDTH)

        arc = Arc([0, 0], width=2+2*d, height=2+2*d, angle=0, theta1=angle_2, theta2=angle_4, color='black', lw=utils.LINE_WIDTH)

        axs[i].add_patch(arc)

    ax.axis('off')


fileItems = [r'Result\Shp\Global\GlobalAreas.shp', 'Result\Shp\Global\GlobalAreas_Raster.tif']

if __name__ == "__main__":
    plt.rcParams['font.family'] = 'Arial'
    f, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=(17/2.54, 6/2.54), dpi=600)

    proportion, polycentric = Proportion()

    PlotPie(ax[1])

    ax[0].axis('off')
    ax[0].annotate("a", xy=(0.005, 1.05), xycoords="axes fraction", fontsize=utils.ORDER_SIZE, weight="bold")
    ax[1].annotate("b", xy=(0.01, 1.05), xycoords="axes fraction", fontsize=utils.ORDER_SIZE, weight="bold")

    plt.tight_layout()
    plt.subplots_adjust(left=0, bottom=0.02, right=1, top=0.88, wspace=0, hspace=0)
    plt.savefig("./Graph/Polycentric/Polycentricity.pdf")

