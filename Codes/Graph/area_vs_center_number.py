# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2026-01

import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import Parameters.plot_utils as utils
import Parameters.consts as const
import statsmodels.api as sm
import numpy as np

"""
Module for the generation of `` Geographical distribution of global urban centers`` plot (Figure 2b and 2c).
"""

def area_per_center(ax):
    """
    Plot the average coverage area by each center versus the center number of the urban area (Figure 2b).

    Urban areas are grouped by the center number and those with centers more than 20 are grouped together as one.
    """

    clusters = gpd.read_file(const.FILE_PATH + fileItems[0])
    clusters['area_per_center'] = clusters['Area'] / clusters['Center_Num']

    bins = list(range(1, 21, 1))
    bins.append(np.inf)
    labels = list(range(1, 21))
    clusters['bin'] = pd.cut(clusters['Center_Num'], bins=bins, labels=labels, right=False)
    group = clusters.groupby('bin', observed=False)

    data = group['area_per_center']
    avg = data.mean()
    std = data.std()

    ax.errorbar(avg.index, avg.values, yerr=std, fmt='^', markersize=np.sqrt(28), markerfacecolor='white',
                markeredgewidth=0.5, capsize=3, capthick=0.5, elinewidth=0.5, color='#AF58BA')

    # Statistics of means
    Average = avg.mean()
    STD = avg.std()
    print(Average, STD)
    ax.hlines(Average, ax.get_xlim()[0], ax.get_xlim()[1], color='#AF58BA', lw=0.8,  zorder=20)

    ax.set_xticks([0, 5, 10, 15, 20], [0, 5, 10, 15, "20+"])
    ax.set_xlim([0, 21])
    ax.set_ylim([0, 150])
    ax.set_xlabel('Center number', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_ylabel(f'Average area per center (km²)', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD + 1)

    utils.FormatAxis(ax)

def area_versus_center(ax):
    """
    Scaling law analysis of the city area relative to the center number (Figure 2c).

    The area and center number are fitted by OLS under the log-log scale.
    """

    clusters = gpd.read_file(const.FILE_PATH + fileItems[0])

    group = clusters.groupby('Center_Num')
    means = group['Area'].mean()
    std = group['Area'].std()

    x_all = means.index.astype(float)
    y_all = means.values.astype(float)

    mask = means.index > 1
    means_fit = means[mask]

    # fitting for polycentric cities
    x_fit = means_fit.index.astype(float)
    y_fit = means_fit.values.astype(float)
    logx_fit = np.log10(x_fit)
    logy_fit = np.log10(y_fit)

    X_fit = sm.add_constant(logx_fit)
    model = sm.OLS(logy_fit, X_fit).fit()

    intercept = model.params[0]
    slope = model.params[1]
    conf_int = model.conf_int(alpha=0.05)
    ci_low, ci_high = conf_int[1, 0], conf_int[1, 1]

    ax.scatter(clusters['Center_Num'], clusters['Area'], marker='s', s=20, edgecolor='#91d5f2', linewidth=0.25, facecolor='none')
    ax.scatter(x_all, y_all, marker='^', s=28, edgecolor='#AF58BA', linewidth=0.5, facecolor='none')

    ax.plot(x_all, 10 ** model.predict(sm.add_constant(np.log10(x_all))), color='#AF58BA', lw=0.8, zorder=10)

    print(f"Fit Data (n>1): $n={len(x_fit)} \: R^2={model.rsquared:.2f}$ \n"
          + rf"$\beta={slope:.4f} intercept={intercept:.4f} \pm {((ci_high - ci_low) / 2):.4f}$")

    ax.set_xlabel(r'log$_{\mathregular{10}}$ Center number', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_ylabel(r'log$_{\mathregular{10}}$ Area of the city (km²)', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD + 1)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xticks([1, 10 ** 0.5, 10 ** 1, 10 ** 1.5, 10 ** 2], [0, 0.5, 1, 1.5, 2])
    ax.set_yticks([1, 1e1, 1e2, 1e3, 1e4], [0, 1, 2, 3, 4])

    ax.set_xlim(10 ** (-0.08), 10 ** 2.2)
    ax.set_ylim(1, 1e4)

    utils.FormatAxis(ax)

fileItems = [f"GHSL\\Global\\globalClusters.shp"]

if __name__ == "__main__":
    f, ax = plt.subplots(1, 2, figsize=(18 / 2.54, 4.5 / 2.54), dpi=600)

    area_per_center(ax[0])
    area_versus_center(ax[1])

    plt.subplots_adjust(top=0.95, bottom=0.15, left=0.06, right=1.0, wspace=0.15, hspace=0.2)
    plt.savefig(r"D:\Research\Graph\GHSL\panel\area_center.pdf")