# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2025-11

import geopandas as gpd
import numpy as np
import pandas as pd
from matplotlib.patches import Arc
from Parameters import consts as const
from Parameters import plot_utils as utils
import matplotlib.pyplot as plt

"""
Module for the generation of the ``Characteristics of monocentric and polycentric cities`` plot (Figure 3b and 3c).
"""

def distribution_global():
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
    clusters = clusters[clusters['Center_Num'] > 0]
    clusters = clusters[['Center_Num', 'Area', 'Pop', 'Bright']]

    bins = [0, 1, 5, 10, 999]
    labels = ['1', '2', '3', '4']

    clusters['Group'] = pd.cut(clusters['Center_Num'], bins=bins, labels=labels, right=True)
    group = clusters.groupby('Group').sum()
    group['Count'] = clusters.groupby('Group').size()
    group = group[columns]
    # Proportion shared by each category
    group_sum = group.sum(axis=0)
    proportion = group.div(group_sum, axis=1)

    # Proportion shared by polycentric
    prop_polycentric = proportion.drop('1')
    polycentric = prop_polycentric.sum(axis=0)

    print(group_sum, '\n', group, '\n', proportion, '\n', polycentric)

    group = group.drop('1')
    print(group.sum(axis=0))
    return proportion, polycentric


def plot_global():
    """
    Pie plot to show the urban resources composition by monocentric and polycentric cities (Figure 3b).
    """

    def plot_pie(ax):
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
            data = proportion[columns[i]]
            wedges, _ = axs[i].pie(data, labels=None, colors=utils.AREAS_COLORS, startangle=70,
                                   wedgeprops={'width': width})

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

            arc = Arc([0, 0], width=2 + 2 * d, height=2 + 2 * d, angle=0, theta1=angle_2, theta2=angle_4, color='black',
                      lw=utils.LINE_WIDTH)

            axs[i].add_patch(arc)

        ax.axis('off')

    f, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=(18 / 2.54, 6.3 / 2.54), dpi=600)
    proportion, polycentric = distribution_global()

    plot_pie(ax[1])

    ax[0].axis('off')
    ax[0].annotate("a", xy=(0.005, 1.05), xycoords="axes fraction", fontsize=utils.ORDER_SIZE, weight="bold")
    ax[1].annotate("b", xy=(-0.15, 1.05), xycoords="axes fraction", fontsize=utils.ORDER_SIZE, weight="bold")

    plt.tight_layout()
    plt.subplots_adjust(left=0, bottom=0.02, right=1, top=0.88, wspace=0, hspace=0)
    plt.savefig('D:\\Research\\Graph\\GHSL\\panel\\polycentricity.png', dpi=600)


def distribution_regional():
    """
    Calculate the composition of each indicator for each income group/geographic region, broken down into categories: mono, low poly, moderate poly, and high poly.

    Args:
        None.

    Returns:
        pd.Dataframe: The composition in each group for each indicator, broken down by city categories.
    """

    clusters = gpd.read_file(const.FILE_PATH + fileItems[1])
    clusters['Group'] = clusters['Country'].map(const.GROUPS_MAP)
    clusters['Cont'] = clusters['Country'].map(const.CONTS_MAP)
    clusters['Region'] = clusters['Country'].map(const.REGIONS_MAP)
    clusters['Pop'] = clusters['Pop'] / 1e8

    # Label the clusters by centers numebr
    labels = pd.cut(clusters['Center_Num'], bins=[0, 1, 5, 10, 999], labels=const.LABELS, right=True)
    clusters['Label'] = labels
    clusters['Label'] = clusters['Label'].astype(str)
    clusters['Count'] = 1 / 1e3

    def decomposition(type):
        comp = clusters.groupby(by=[type, 'Label'])[const.COLUMNS].sum()
        comp = comp.reset_index()
        comp = comp.pivot(index=type, columns='Label', values=const.COLUMNS)
        dfs = {value: comp[value].reindex(index=const.GROUPS if type == 'Group'
        else const.CONTS if type == 'Cont'
        else const.REGIONS, columns=const.LABELS) for value in const.COLUMNS}
        return dfs

    dfs_cont = decomposition('Cont')
    dfs_group = decomposition('Group')
    dfs_region = decomposition('Region')
    print("Composition of each indicator: \n", dfs_cont, "\n", dfs_group)

    return dfs_cont, dfs_group, dfs_region


def plot_regional():
    """
    Bar plot to show composition of characteristics by different city categories for different income groups (Figure 3c).
    """

    def panel_barh(ax, df, min_val, fmt, y_ticks=False, y_ticklabel=False):
        """
        Plot horizontal stacked bars to show composition of characteristics by different city categories.
        """
        df = df.div(df.sum(axis=1), axis=0) * 100
        width = 0.7
        df.plot(kind='barh', stacked=True, width=width, color=utils.AREAS_COLORS_L, edgecolor='None', linewidth=0,
                ax=ax,
                zorder=10)

        # Edge
        for i, total in enumerate(df.sum(axis=1)):
            ax.barh(i, width=total, height=width, color='None', edgecolor='black', lw=0.3, zorder=15)

        # Text
        for i, (index, row) in enumerate(df.iterrows()):
            left = 0
            for idx, value in enumerate(row):
                if value > min_val:
                    ax.text(left + value / 2, i, fmt(value), va='center', ha='center',
                            fontsize=utils.TEXT_SIZE, color='black' if idx == 0 else 'white', zorder=15)
                left += value

        utils.FormatAxis(ax, [], spineColor='black', spineWidth=0.3, y_ticks=y_ticks, y_labels=y_ticklabel)

        gridColor = '#E3E2E6'
        ax.grid(axis='x', which='major', color=gridColor, linewidth=0.5, zorder=3)
        ax.invert_yaxis()

        ax.set_xticks(ticks=[0, 20, 40, 60, 80, 100], labels=['0', '20%', '40%', '60%', '80%', '100%'])
        ax.set_xlim(0 - 0.2, 100 + 0.2)
        ax.tick_params(axis='both', length=0, pad=2)

        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.legend().remove()


    fig, ax = plt.subplots(2, 2, figsize=(18 / 2.54, 7 / 2.54), dpi=600)
    dfs_cont, dfs_group, dfs_region = distribution_regional()

    panel_barh(ax[0, 0], dfs_group['Count'], min_val=5, fmt=lambda v: f"{v:.0f}%",
               y_ticks=True, y_ticklabel=True)
    panel_barh(ax[0, 1], dfs_group['Area'], min_val=5, fmt=lambda v: f"{v:.0f}%",)

    panel_barh(ax[1, 0], dfs_group['Pop'], min_val=5, fmt=lambda v: f"{v:.0f}%",
               y_ticks=True, y_ticklabel=True)
    panel_barh(ax[1, 1], dfs_group['Bright'], min_val=5, fmt=lambda v: f"{v:.0f}%",)

    plt.subplots_adjust(top=0.95, right=0.98, left=0.14, bottom=0.10, wspace=0.07, hspace=0.20)
    plt.savefig(r"D:\Research\Graph\GHSL\panel\distribution.pdf")


fileItems = [r'GHSL\Global\globalClusters.shp', r'GHSL\Global\splitClusters.shp']
columns = ['Count', 'Area', 'Pop', 'Bright']

if __name__ == "__main__":
    plot_regional()
