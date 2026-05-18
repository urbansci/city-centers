# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2026-04

import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm

import Parameters.plot_utils as utils
import Parameters.consts as const

"""
Module for the generation of scaling analysis figure (Figure 2).

Panel (a): Average area per center versus center number.
Panel (b): Scaling analysis between total built-up area and center number (log-log).
Panel (c): Robustness analysis — average area per center across minimum contour areas.
Panel (d): Robustness analysis — scaling exponent across minimum contour areas.
"""

# ── Constants ──────────────────────────────────────────────────────────────

CENTERS_PATH = const.FILE_PATH + r'2020\Center\center_a4.shp'
CLUSTERS_PATH = const.FILE_PATH + r'2020\Input\GHSL_UCDB_MTUC_2020_GLOBE_R2024_WGS84.shp'
OUTPUT_PATH = r'D:\Research\Graph\GHSL\panel\scaling.png'

YEAR_ITEMS = {
    '2020': {
        'base': const.FILE_PATH + r'2020\Input\GHSL_UCDB_MTUC_2020_GLOBE_R2024_WGS84.shp',
        'centers': const.FILE_PATH + r'2020\Center\center_a{parm}.shp',
    },
    '2015': {
        'base': const.FILE_PATH + r'2015\Input\GHSL_UCDB_MTUC_2015_GLOBE_R2024_WGS84.shp',
        'centers': const.FILE_PATH + r'2015\Center\center_a{parm}.shp',
    },
}


# ── Helper functions ────────────────────────────────────────────────────

def _build_clusters(centers_path, base_path):
    """Load clusters and assign center count; return only clusters with centers."""
    centers = gpd.read_file(centers_path)
    base_clusters = gpd.read_file(base_path)
    center_num = centers.groupby('ID_UC_G0').size()
    base_clusters['Center_Num'] = base_clusters['ID_UC_G0'].map(center_num).fillna(0).astype(int)
    return base_clusters[base_clusters['Center_Num'] > 0]


# ── Panel plot functions ────────────────────────────────────────────────

def _plot_area_per_center(ax, clusters, min_centers=4):
    """Panel (a): average area per center vs center number with error bars."""
    clusters = clusters.copy()
    clusters['area_per_center'] = clusters['Area'] / clusters['Center_Num']
    bins = list(range(1, 21)) + [np.inf]
    labels = list(range(1, 21))
    clusters['bin'] = pd.cut(clusters['Center_Num'], bins=bins, labels=labels, right=False)
    data = clusters.groupby('bin', observed=False)['area_per_center']
    avg, std = data.mean(), data.std()

    ax.errorbar(avg.index, avg.values, yerr=std, fmt='^', markersize=np.sqrt(28),
                markerfacecolor='white', markeredgewidth=0.5, capsize=3,
                capthick=0.5, elinewidth=0.5, color='#AF58BA')

    mean_val = avg[avg.index >= min_centers].mean()
    ax.hlines(mean_val, ax.get_xlim()[0], ax.get_xlim()[1], color='#AF58BA', lw=0.8, zorder=20)
    ax.text(0.90, 0.85, f'Mean = {mean_val:.1f} km²\n',
            transform=ax.transAxes, fontsize=utils.TEXT_SIZE, va='top', ha='right')

    ax.set_xlabel('Center number', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_ylabel('Average area per center (km²)', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD + 1)
    ax.set_ylim(0, 150)
    ax.set_yticks([0, 50, 100, 150])
    ax.set_xlim(0, 21)
    ax.set_xticks([1, 5, 10, 15, 20], ['1', '5', '10', '15', '20+'])
    utils.FormatAxis(ax)


def _plot_area_versus_center(ax, clusters, min_centers=4):
    """Panel (b): log-log scaling of total area vs center number with OLS fit."""
    clusters = clusters.copy()
    means = clusters.groupby('Center_Num')['Area'].mean()

    x_all = means.index.astype(float)
    y_all = means.values.astype(float)
    fit_means = means[means.index >= min_centers]

    logx_fit = np.log10(fit_means.index.astype(float))
    logy_fit = np.log10(fit_means.values.astype(float))
    model = sm.OLS(logy_fit, sm.add_constant(logx_fit)).fit()
    slope = model.params[1]
    ci_low, ci_high = model.conf_int(alpha=0.05)[1]

    ax.scatter(clusters['Center_Num'], clusters['Area'], marker='s', s=20,
               edgecolor='#91d5f2', linewidth=0.25, facecolor='none')
    ax.scatter(x_all, y_all, marker='^', s=28,
               facecolors='white', edgecolors='#AF58BA', linewidths=0.5, zorder=5)
    ax.plot(x_all, 10 ** model.predict(sm.add_constant(np.log10(x_all))),
            color='#AF58BA', lw=0.8, zorder=10)

    stat_text = (rf'$\beta$ = {slope:.2f} $\pm$ {((ci_high - ci_low) / 2):.2f}' + '\n'
                 + f'$R$²= {model.rsquared:.2f}')
    ax.text(0.7, 0.4, stat_text, transform=ax.transAxes, fontsize=utils.TEXT_SIZE,
            va='top', ha='left', linespacing=1.8)

    ax.set_xlabel(r'log$_{\mathregular{10}}$ Center number', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_ylabel(r'log$_{10}$ Area of the city' + ' (km²)', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD + 1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xticks([1, 10 ** 0.5, 10 ** 1, 10 ** 1.5, 10 ** 2], [0, 0.5, 1, 1.5, 2])
    ax.set_yticks([1, 1e1, 1e2, 1e3, 1e4], [0, 1, 2, 3, 4])
    ax.set_xlim(10 ** (-0.08), 10 ** 2.3)
    ax.set_ylim(1, 1e4)
    utils.FormatAxis(ax)


def _plot_robustness(ax_c, ax_d, min_centers=4):
    """Panels (c/d): robustness of mean area and scaling exponent across min contour areas."""
    parms = list(range(1, 11))
    series = [
        {'year': '2015', 'color': '#009ADE', 'marker': 'o', 'label': '2015'},
        {'year': '2020', 'color': '#AF58BA', 'marker': 's', 'label': '2020'},
    ]

    for sr in series:
        yi = YEAR_ITEMS[sr['year']]
        mc_means, mc_stds = [], []
        mc_betas, mc_ci_lows, mc_ci_highs = [], [], []

        for parm in parms:
            clusters = _build_clusters(yi['centers'].format(parm=parm), yi['base'])
            clusters['area_per_center'] = clusters['Area'] / clusters['Center_Num']
            avg = clusters.groupby('Center_Num')['area_per_center'].mean()
            valid_avg = avg[avg.index >= min_centers]
            mc_means.append(valid_avg.mean())
            mc_stds.append(valid_avg.std())

            group_means = clusters.groupby('Center_Num')['Area'].mean()
            fit_means = group_means[group_means.index >= min_centers]
            x_fit = fit_means.index.astype(float)
            y_fit = fit_means.values.astype(float)
            if len(x_fit) >= 3:
                logx = np.log10(x_fit)
                logy = np.log10(y_fit)
                model = sm.OLS(logy, sm.add_constant(logx)).fit()
                ci = model.conf_int(alpha=0.05)
                mc_betas.append(model.params[1])
                mc_ci_lows.append(ci[1, 0])
                mc_ci_highs.append(ci[1, 1])
            else:
                mc_betas.append(np.nan)
                mc_ci_lows.append(np.nan)
                mc_ci_highs.append(np.nan)

        c, m = sr['color'], sr['marker']

        def _plot_lm(ax, x, y, color=c, marker=m, label=None):
            ax.plot(x, y, color=color, linewidth=0.8, zorder=2)
            ax.scatter(x, y, marker=marker, s=16, facecolor='white',
                       edgecolor='white', linewidth=0, zorder=8)
            ax.scatter(x, y, marker=marker, s=16, facecolor='none',
                       edgecolor=color, linewidth=0.8, zorder=10, label=label)

        ax_c.errorbar(parms, mc_means, yerr=mc_stds, fmt='none', ecolor=c,
                      capsize=2, capthick=0.5, elinewidth=0.5, zorder=2)
        _plot_lm(ax_c, parms, mc_means, color=c, marker=m, label=sr['label'])

        ax_d.fill_between(parms, mc_ci_lows, mc_ci_highs, facecolor=c,
                          edgecolor='none', alpha=0.1, zorder=1)
        _plot_lm(ax_d, parms, mc_betas, color=c, marker=m, label=sr['label'])

    ax_c.set_ylabel('Average area per center (km²)', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax_c.set_xlabel('Minimum contour area (km²)', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax_c.set_ylim(0, 150)
    ax_c.set_xticks(ticks=range(1, 11))
    ax_c.legend(fontsize=utils.LEGEND_SIZE, markerscale=1, loc='lower right', handletextpad=0.3)
    utils.FormatAxis(ax_c)

    ax_d.axhline(1, color='gray', lw=0.8, ls='--', zorder=0)
    ax_d.set_ylabel(r'$\beta$', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD - 1)
    ax_d.set_xlabel('Minimum contour area (km²)', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax_d.set_ylim(0.8, 1.2)
    ax_d.set_xticks(ticks=range(1, 11))
    utils.FormatAxis(ax_d)


# ── Main ────────────────────────────────────────────────────────────────

def figure2():
    """Compose and save the complete Figure 2 (2x2 panels)."""
    WIDTH_CM, HEIGHT_CM = 13.0, 8.0

    print("Loading data...")
    clusters = _build_clusters(CENTERS_PATH, CLUSTERS_PATH)

    fig, axes = plt.subplots(2, 2, figsize=(WIDTH_CM / 2.54, HEIGHT_CM / 2.54))

    print("Panel (a/b): Area analysis...")
    _plot_area_per_center(axes[0, 0], clusters)
    _plot_area_versus_center(axes[0, 1], clusters)

    print("Panel (c/d): Robustness analysis...")
    _plot_robustness(axes[1, 0], axes[1, 1])

    fig.subplots_adjust(left=0.10, right=0.97, top=0.96, bottom=0.09, wspace=0.25, hspace=0.35)

    print(f"Saving to {OUTPUT_PATH}...")
    fig.savefig(OUTPUT_PATH, dpi=600, facecolor='white')
    plt.close(fig)
    print("Done.")


if __name__ == "__main__":
    figure2()
