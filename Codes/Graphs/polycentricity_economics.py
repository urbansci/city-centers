# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2024-12

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
import Parameters.plot_utils as utils
import Parameters.consts as const
from scipy.stats import pearsonr

"""
Module for the generation of the ``Urban center and economic development`` figure (Figure 4). 
"""

def Distribution():
    """
    Calculate the composition of each indicator for each income group, broken down into categories: mono, low poly, moderate poly, and high poly "

    Args:
        None.

    Returns:
        pd.Dataframe: The composition in each group for each indicator, broken down by city categories.
    """

    clusters = gpd.GeoDataFrame()
    globalClusters = gpd.read_file(const.FILE_PATH + fileItems[3])
    for group in const.GROUPS:
        for i in const.GROUPS_DICT[group]:
            _clusters = globalClusters[globalClusters['Country'] == i]
            _clusters['Group'] = group
            clusters = pd.concat([clusters, _clusters], ignore_index=True)

    # Label the clusters by centers numebr
    labels = pd.cut(clusters['Center_Num'], bins=[0, 1, 5, 10, 999], labels=const.LABELS, right=True)
    clusters['Label'] = labels
    clusters['Count'] = 1
    clusters['Label'] = clusters['Label'].astype(str)
    clusters.to_file(const.FILE_PATH + fileItems[4])

    # Composition of each indicator for each group
    proportion = clusters.groupby(by='Group').apply(lambda x: x.groupby(by='Label')[const.COLUMNS].sum() / x[
        const.COLUMNS].sum())
    proportion = proportion.reset_index()
    proportion = proportion.pivot(index='Group', columns='Label', values=const.COLUMNS)
    dfs = {value: proportion[value].reindex(index=const.GROUPS, columns=const.LABELS) for value in const.COLUMNS}
    print("Composition of each indicator: \n", dfs)

    return dfs

def NationalStats():
    """
    Get the national statistics for countries with more than 10 million population, including GNI per capita, polycentric proportion, urbanized area and belonging continent.

    Args:
        None.

    Returns:
        DataFrame: The dataframe containing above statistics for each country.
    """

    econ = pd.read_excel(const.FILE_PATH + fileItems[5], sheet_name=0, skiprows=3)
    econ = econ[['Country Code', '2020']]
    econ = econ.dropna(subset=['2020'])  # Drop Data Nan
    econ = econ.rename(columns={'2020': 'Econ'})

    pop = pd.read_excel(const.FILE_PATH + r"Data\Class\POP.xls", sheet_name=0, skiprows=3)
    pop = pop[['Country Code', '2020']]
    pop = pop[pop['2020'] > 1e7]
    pop = pop.rename(columns={'2020': 'Pop'})
    econ = pd.merge(econ, pop, left_on='Country Code', right_on='Country Code', suffixes=("", "_"), how='inner')

    countries = gpd.read_file(const.FILE_PATH + fileItems[2])
    contMap = {item: cont for cont, items in const.CONTS_DICT.items() for item in items}
    countries['FID'] = countries.index
    countries['Cont'] = countries.index.map(contMap)

    stats = pd.merge(countries, econ, left_on='SOC', right_on='Country Code', how='inner')
    stats['Poly_Ptn'] = 0.0
    stats['Cluster_Num'] = 0
    stats['Urban_Area'] = 0.0

    globalClusters = gpd.read_file(const.FILE_PATH + fileItems[3])
    for i, country in stats.iterrows():
        index = stats.loc[i, 'FID']
        _clusters = globalClusters[globalClusters['Country'] == index]
        stats.loc[i, 'Poly_Ptn'] = len(_clusters[_clusters['Center_Num'] > 1]) / len(_clusters) if len(_clusters) > 0 else 0
        stats.loc[i, 'Cluster_Num'] = len(_clusters)
        stats.loc[i, 'Urban_Area'] = _clusters['Area'].sum()

    # stats[['geometry', 'SOC', 'Econ', 'Cont', 'Pop', 'Poly_Ptn', 'Cluster_Num', 'Urban_Area']].to_file(
    #     const.FILE_PATH + r'Result\\Shp\\Global\\NationStats_.shp')
    return stats[['SOC', 'Econ', 'Poly_Ptn', 'Cluster_Num', 'Urban_Area', 'Cont']]


def BarH(df, columns):
    """
    Plot horizontal stacked bars to show composition of characteristics by different city categories.

    Args:
        df (dict[str, DataFrame]): Dictionary of dataframes for each column.
        columns (list[str]): Names of columns to be plot in sequence (Only two).

    Returns:
        None.
    """

    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(2, 1, figsize=(11 / 2.54, 5 / 2.54), dpi=600)

    width = 0.7
    df[columns[0]].plot(kind='barh', stacked=True, width=width, color=utils.AREAS_COLORS, edgecolor='None', linewidth=0, ax=ax[0], zorder=10)
    df[columns[1]].plot(kind='barh', stacked=True, width=width, color=utils.AREAS_COLORS, edgecolor='None', linewidth=0, ax=ax[1], zorder=10)

    # Edge
    for i in range(len(const.COLUMNS)):
        ax[0].barh(i, width=1, height=width, color='None', edgecolor='black', lw=0.3, zorder=15)
    for i in range(len(const.COLUMNS)):
        ax[1].barh(i, width=1, height=width, color='None', edgecolor='black', lw=0.3, zorder=15)

    # Text
    for i, (index, row) in enumerate(df[columns[0]].iterrows()):
        left = 0
        for idx, value in enumerate(row):
            if value > 0.05:
                ax[0].text(left + value / 2, i, f"{value:.0%}", va='center', ha='center',
                           fontsize=utils.TEXT_SIZE, color='black' if idx == 0 else 'white', zorder=15)
            left += value

    for i, (index, row) in enumerate(df[columns[1]].iterrows()):
        left = 0
        for idx, value in enumerate(row):
            if value > 0.05:
                ax[1].text(left + value / 2, i, f"{value:.0%}", va='center', ha='center',
                           fontsize=utils.TEXT_SIZE, color='black' if idx == 0 else 'white', zorder=15)
            left += value

    gridColor = '#969696'
    ax[0].grid(axis='x', which='major', color=gridColor, linewidth=0.5, zorder=3)
    ax[1].grid(axis='x', which='major', color=gridColor, linewidth=0.5, zorder=3)
    ax[0].invert_yaxis()
    ax[1].invert_yaxis()

    ax[0].set_yticklabels(const.GROUPS, fontsize=utils.LABEL_SIZE)
    ax[1].set_yticklabels(const.GROUPS, fontsize=utils.LABEL_SIZE)
    ax[0].set_ylabel("")
    ax[1].set_ylabel("")
    ax[1].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax[1].set_xticklabels(['0%', '20%', '40%', '60%', '80%', '100%'], fontsize=utils.LABEL_SIZE)
    ax[0].set_xlim([0, 1])
    ax[1].set_xlim([0, 1])

    utils.FormatAxis(ax[0], ['left'], spineColor=gridColor, x_ticks=False, y_ticks=False, x_labels=False)
    utils.FormatAxis(ax[1], ['left'], spineColor=gridColor, x_ticks=False, y_ticks=False)
    ax[0].legend().remove()
    ax[1].legend().remove()

    plt.subplots_adjust(top=0.95, right=0.975, left=0.20, bottom=0.05, wspace=0, hspace=0.1)

    plt.savefig(r".\Graph\Income\HorizontalBars.png")

def Scatter(df):
    """
    Plot the scatter map of polycentric proportion relative to the national income (GNI per capita), size-coded by number of urban areas.

    Args:
        df (DataFrame): The dataframe containing GNI, number of urban areas as well as polycentric proportion for each country.

    Returns:
         None.
    """

    plt.rcParams['font.family'] = 'Arial'
    fig, ax = plt.subplots(figsize=(11 / 2.54, 5.6 / 2.54), dpi=600)

    Xmin = 10**2.5
    Xmax = 1e5
    bins = [0, 1045, 4095, 12695, 9999999]
    labels = pd.cut(df['Econ'], bins=bins, labels=np.arange(4), right=True)
    df['Label'] = labels

    contColors = {'NA': '#ff7f0e', 'OC': '#2ca02c', 'AF': '#E3E2E6', 'SA': '#8c564b', 'EU': '#9467bd', 'AS': '#1f77b4'}
    contOrder = {'NA': 13, 'OC': 13, 'EU': 13, 'AF': 8, 'SA': 10, 'AS': 9}

    max_size = 80
    min_size = 1

    max_num = df['Urban_Area'].max()
    min_num = df['Urban_Area'].min()
    def transSize(x):
        # Size by log transfomrm
        normalized_size = (np.log10(x) - np.log10(min_num)) / (np.log10(max_num) - np.log10(min_num))
        s = min_size + normalized_size * (max_size - min_size)
        return s

    s = df['Urban_Area'].apply(lambda x: transSize(x))
    ax.scatter(df['Econ'], df['Poly_Ptn'], marker='o', facecolor='None',
                   edgecolor='black', linewidths=0.3, s=s, zorder=15)

    # for i in range(len(df)):
    #     ax.text(
    #         df['Econ'].iloc[i], df['Poly_Ptn'].iloc[i],
    #         str(df['SOC'].iloc[i]), fontsize=3, ha='center', va='center',
    #         zorder=20
    #     )

    for cont in const.CONTS:
        _df = df[df['Cont'] == cont]
        _s = s[_df.index.tolist()]
        ax.scatter(_df['Econ'], _df['Poly_Ptn'], marker='o', c=contColors[cont],
                   edgecolor='black', linewidths=0, s=_s, zorder=contOrder[cont])

    # Dividing line for income group
    for x in bins[1: -1]:
        ax.axvline(x=x, linestyle=(0, (14, 15)), color=utils.GRID_COLOR, linewidth=0.5)

    # Linear Regression
    x = df['Econ'].values
    y = df['Poly_Ptn'].values
    x_log = np.log10(x)

    # Get fitting performance
    x_log_with_const = sm.add_constant(x_log)
    ols_model = sm.OLS(y, x_log_with_const).fit()
    print("Fitting:", ols_model.summary())
    corr_coefficient, p_value = pearsonr(x_log, y)
    print(f"Pearson correlation coefficient: {corr_coefficient}")
    print(f"P-value: {p_value}")

    x_log = x_log.reshape(-1, 1)
    model = LinearRegression()
    model.fit(x_log, y)
    x_fit = [[100], [100000]]
    ax.plot(x_fit, model.predict(np.log10(x_fit)), color='black', linewidth=0.6, zorder=5)

    ax.set_xlabel(r'$\log_{10}\mathrm{GDP \, per \, capita\, } (\$)$', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD-1.5)
    ax.set_ylabel(f'Polycentric proportion', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD+2)

    ax.set_xscale('log')
    ax.set_xlim(Xmin, Xmax)
    ax.set_ylim(None, 0.5)
    ax.set_xticks(ticks=[10**2.5, 1e3, 10**3.5, 1e4, 10**4.5, 1e5])
    ax.set_xticklabels(labels=[2.5, 3, 3.5, 4, 4.5, 5])
    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])

    utils.FormatAxis(ax, ['left', 'bottom'])
    ax.minorticks_off()

    plt.subplots_adjust(top=0.95, right=0.98, left=0.08, bottom=0.10, wspace=0.2, hspace=0.4)
    plt.savefig(r'.\Graph\Income\Economic.png')

fileItems = ["Result\\Shp\\Center\\", "Result\\Shp\\Cluster\\", "Data\\World Country\\CAS Country.shp", "Result\\Shp\\Global\\GlobalAreas.shp", 
            "Result\\Shp\\Global\\Areas_Income.shp", r"Data\Class\GNP.PCAP.xls"]

if __name__ == "__main__":

    dfs = Distribution()
    BarH(dfs, ['Count', 'Pop'])

    proPtn = NationalStats()
    Scatter(proPtn)



