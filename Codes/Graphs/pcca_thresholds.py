# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2024-12

import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from collections import Counter
import Parameters.consts as const
import Parameters.plot_utils as utils

"""
Module for the generation of the ``Parameter estimation in PCCA`` figure (Extended Data Figure 2). 
"""

def PlotPC(ax):
    """
    Plot the percolation graph of specific countries.
    """

    for i, country in enumerate(countries):
        optimalThres = const.OPTIMAL_THRESHOLD[country]
        map = {2:3, 232:246, 150:159, 219:232}
        f = open(const.FILE_PATH + fileItems[0] + f"{map[country]}", 'rb')
        portion = pickle.load(f)

        ax.plot(thresholds, portion[:21], color=utils.COUNTRY_COLORS[i],  linewidth=0.8, label=f'{names[country]}', alpha=1)
        ax.scatter(optimalThres, portion[int(optimalThres/0.5)], color=utils.COUNTRY_COLORS[i], marker='^', s=utils.MARKER_SIZE, edgecolor='grey', linewidth=0.3, zorder=4)

        # index = list(thresholds)
        # index.remove(optimalThres)
        # portion.pop(int(optimalThres/0.5))

    ax.set_xticks(np.arange(0, max(thresholds) + 0.5, 1))
    ax.set_xlabel(r'Threshold (nW$\cdot$cm$^{-2}$$\cdot$sr$^{-1}$)', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_ylabel('Largest Cluster Proportion', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_ylim(-0.03, 0.53)

    ax.legend(loc='best', frameon=False, fontsize=utils.LEGEND_SIZE)
    utils.FormatAxis(ax)


def func(x, a, b):
    return a*x**(-b)

def PlotRankCurve(ax):
    """
    Plot the rank-size curves of area distribution for specific countries.
    """

    for i, country in enumerate(countries):
        cluster = gpd.read_file(const.FILE_PATH + fileItems[1] + f'{country}.shp')
        area = cluster['Area'].copy()
        area.sort_values(ascending=False, inplace=True)

        ax.scatter(np.arange(1, len(area)+1), area, marker='o', facecolor='white', edgecolor=utils.COUNTRY_COLORS[i], linewidth=0.3,
                   s=2, alpha=1, zorder=3)

        # Fit Curve
        x = np.arange(10, len(area)+1, dtype=float)
        y = area.values[9: ] # regress from the 10th point
        popt, pcov = curve_fit(func, x, y)
        a, b = popt
        y_pred = func(x, a, b)
        r2 = r2_score(y, y_pred)

        ax.plot([1, (a/10)**(1/b)], [a, 10], linestyle=(0, (6, 3)), color=utils.COUNTRY_COLORS[i], linewidth=0.8,
                 label=rf'$\beta$={b:.2f}  $R^2$={r2:.2f}', alpha=1, zorder=3)

    ax.set_xlabel('Rank', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_ylabel(f'Area (km$^2$)', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.legend(loc='best', frameon=False, fontsize=utils.LEGEND_SIZE)
    utils.FormatAxis(ax)
    ax.tick_params(axis='x', which='minor', length=0)
    ax.tick_params(axis='y', which='minor', length=0)

def PlotPDF(ax):
    YMin = -7
    value_counts = Counter(const.OPTIMAL_THRESHOLD.values())

    values = np.arange(0.5, 2.5, 0.5)
    counts = [value_counts[val] for val in values]
    residue = len(const.OPTIMAL_THRESHOLD)-sum(counts)-value_counts[None]
    counts.append(residue)
    print(counts)

    ax.bar(np.arange(0.5, 3.0, 0.5), counts, width=0.3, bottom=YMin, align='center', color=utils.BAR_COLOR, edgecolor='black', linewidth=0.3, alpha=None, label=None, zorder=3)

    ax.set_xlabel(r'Threshold (nW$\cdot$cm$^{-2}$$\cdot$sr$^{-1}$)', fontsize=utils.LABEL_SIZE,labelpad=utils.LABEL_PAD)
    ax.set_ylabel('Frequency', fontsize=utils.LABEL_SIZE,labelpad=utils.LABEL_PAD)

    ax.set_xticks(np.arange(0.5, 3.0, 0.5), labels=[0.5, 1.0, 1.5, 2.0, '>2.0'])
    ax.set_ylim(YMin, None)

    utils.FormatAxis(ax)


countries = [2, 232, 150, 219]
names = {2: 'US', 232: 'China', 150: 'Germany', 219: 'Japan'}
thresholds = np.arange(0, 10.5, 0.5)
fileItems = ['Temp/portionR', 'Result/Shp/Cluster/']

if __name__ == "__main__":

    plt.rcParams['font.family'] = 'Arial'
    f, ax = plt.subplots(2, 3, gridspec_kw={'height_ratios': [2, 3]}, figsize=(17/2.54, 13/2.54), dpi=600)

    PlotPC(ax[0, 0])
    PlotRankCurve(ax[0, 1])
    PlotPDF(ax[0, 2])

    ax[1, 0].axis('off')
    ax[1, 1].axis('off')
    ax[1, 2].axis('off')

    ax[0, 0].annotate("a", xy=(-0.2, 1.1), xycoords="axes fraction", fontsize=utils.ORDER_SIZE, weight="bold")
    ax[0, 1].annotate("b", xy=(-0.2, 1.1), xycoords="axes fraction", fontsize=utils.ORDER_SIZE, weight="bold")
    ax[0, 2].annotate("c", xy=(-0.2, 1.1), xycoords="axes fraction", fontsize=utils.ORDER_SIZE, weight="bold")
    ax[1, 0].annotate("d", xy=(-0.2, 1.1), xycoords="axes fraction", fontsize=utils.ORDER_SIZE, weight="bold")

    plt.tight_layout()

    plt.savefig("./Graph/Threshold/Threshold.png")

