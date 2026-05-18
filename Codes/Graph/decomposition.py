# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2026-04

"""
Decomposition of scaling laws (Figure 3 & Supplementary Figure 11).

Figure 3: Five-panel decomposition showing that the linear center–area scaling
    arises from two independent sublinear relationships (A~P and N~P).
    (a) Per-country area vs center number with highlighted series.
    (b) Center number vs population (US cities, log-log).
    (c) Area vs population (US cities, log-log).
    (d) Histogram of scaling exponents across countries.
    (e) Scatter of gamma vs alpha colored by continent.

Supplementary Figure 11: Robustness of per-country scaling across minimum
    contour area thresholds (2, 3, 4 km²).
"""

import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import Parameters.plot_utils as utils
import Parameters.consts as const
import statsmodels.api as sm
import numpy as np
from rasterstats import zonal_stats
import rasterio
from adjustText import adjust_text

# ── Constants ─────────────────────────────────────────────────────────────────

COUNTRY_SERIES = {
    0:   ('China',   '#d62728', 'o'),
    142: ('US',      '#AF58BA', '^'),
    1:   ('Japan',   '#8c564b', 's'),
    93:  ('UK',      '#ff7f0e', 'D'),
    4:   ('Russia',  '#2ca02c', 'v'),
}

CONTINENT_STYLE = {
    'Asia':     ('#d62728', 'o'),
    'Europe':   ('#8c564b', 's'),
    'Americas': ('#AF58BA', '^'),
    'Africa':   ('#ff7f0e', 'D'),
    'Oceania':  ('#2ca02c', 'v'),
}

ALL_COUNTRY_CONTINENT = {
    0:   'Asia',      # China
    142: 'Americas',  # US
    5:   'Asia',      # India
    4:   'Europe',    # Russia
    152: 'Americas',  # Brazil
    1:   'Asia',      # Japan
    170: 'Asia',      # Indonesia
    93:  'Europe',    # UK
    192: 'Americas',  # Canada     
    177: 'Oceania',   # Australia
}

FINAL_COUNTRIES = {192, 1, 93, 5, 0, 142, 170, 4, 152, 177}
# Canada, Japan, UK, India, China, US, Indonesia, Russia, Brazil, Australia

POP_RASTER = const.FILE_PATH + r"Data\Pop\WorldPop\ppp_2020_1km_Aggregated.tif"
WB_BOUNDARIES = const.FILE_PATH + r"Data\World Country\World Bank Boundaries.gpkg"

CENTERS_PATH = const.FILE_PATH + r'2020\Center\center_a{parm}.shp'
CLUSTERS_PATH = const.FILE_PATH + r'2020\Input\GHSL_UCDB_MTUC_2020_GLOBE_R2024_WGS84.shp'

MIN_CENTERS = 4

_pop_series = None
_wb_names = None


# ── Data ──────────────────────────────────────────────────────────────────────

def _get_wb_names():
    """Load World Bank country name mapping (cached)."""
    global _wb_names
    if _wb_names is None:
        wb = gpd.read_file(WB_BOUNDARIES, ignore_geometry=True)
        _wb_names = dict(zip(wb['id'], wb['NAM_0']))
    return _wb_names


def _build_clusters(parm=4):
    """Load clusters with center count and population for a given contour area threshold."""
    global _pop_series
    centers_path = CENTERS_PATH.format(parm=parm)
    centers = gpd.read_file(centers_path)
    base = gpd.read_file(CLUSTERS_PATH)
    center_num = centers.groupby('ID_UC_G0').size()
    base['Center_Num'] = base['ID_UC_G0'].map(center_num).fillna(0).astype(int)
    clusters = base[base['Center_Num'] > 0].copy()

    if _pop_series is None:
        with rasterio.open(POP_RASTER) as src:
            raster_crs = src.crs
        gdf_proj = base.to_crs(raster_crs) if base.crs != raster_crs else base
        stats = zonal_stats(gdf_proj, POP_RASTER, stats=['sum'])
        _pop_series = pd.Series([s['sum'] if s['sum'] is not None else 0 for s in stats], index=base['ID_UC_G0'])

    clusters['Pop'] = clusters['ID_UC_G0'].map(_pop_series).fillna(0)
    return clusters


def _country_betas(clusters):
    """Compute scaling exponents (A~N, N~P, A~P) per country via OLS."""
    clusters = clusters[(clusters['Area'] > 0) & (clusters['Pop'] > 0)].copy()
    clusters = clusters[clusters['Country'].isin(FINAL_COUNTRIES)]
    results = []
    wb_names = _get_wb_names()
    for country, group in clusters.groupby('Country'):
        sub = group[group['Center_Num'] >= MIN_CENTERS]
        
        means_a = sub.groupby('Center_Num')['Area'].mean()
        log_a_gm = np.log10(means_a.values)
        log_n_gm = np.log10(means_a.index.values.astype(float))    
        mdl_an = sm.OLS(log_a_gm, sm.add_constant(log_n_gm)).fit()

        means_p = sub.groupby('Center_Num')['Pop'].mean()    
        log_p_gm = np.log10(means_p.values)
        mdl_np = sm.OLS(log_n_gm, sm.add_constant(log_p_gm)).fit()

        log_p = np.log10(sub['Pop'].values.astype(float))
        log_a = np.log10(sub['Area'].values.astype(float))
        mdl_ap = sm.OLS(log_a, sm.add_constant(log_p)).fit()
        
        ci_an = mdl_an.conf_int(alpha=0.05)
        ci_np = mdl_np.conf_int(alpha=0.05)
        ci_ap = mdl_ap.conf_int(alpha=0.05)
        results.append({
            'Country': int(country),
            'Name': wb_names.get(int(country), str(int(country))),
            'n': len(sub),
            'beta': mdl_an.params[1],
            'ci_AN': (ci_an[1, 1] - ci_an[1, 0]) / 2,
            'r2_AN': mdl_an.rsquared,
            'gamma': mdl_np.params[1],
            'ci_NP': (ci_np[1, 1] - ci_np[1, 0]) / 2,
            'alpha': mdl_ap.params[1],
            'ci_AP': (ci_ap[1, 1] - ci_ap[1, 0]) / 2,
            'ratio': mdl_ap.params[1] / mdl_np.params[1],
        })
    return pd.DataFrame(results)


# ── Panel plots ────────────────────────────────────────────────────────────────
def _plot_multi_country_area_vs_center(ax, clusters, bg_country_ids=None,
                                       xlim=(10**0.4, 10**2.2), ylim=(10**2, 10**5),
                                       xticks=None, yticks=None):
    """Fig. 3a: per-country area vs center number with highlighted series."""
    if bg_country_ids is not None:
        for cid in bg_country_ids:
            if cid in COUNTRY_SERIES:
                continue
            sub = clusters[(clusters['Country'] == cid) & (clusters['Center_Num'] >= MIN_CENTERS)]
            if len(sub) == 0:
                continue
            means = sub.groupby('Center_Num')['Area'].mean()
            ax.scatter(means.index.astype(float), means.values.astype(float),
                       marker='o', s=25, edgecolor='#c8c8c8', linewidth=0.6,
                       facecolor='none', zorder=3)

    entries = []
    for cid, (name, color, marker) in COUNTRY_SERIES.items():
        sub = clusters[(clusters['Country'] == cid) & (clusters['Center_Num'] >= MIN_CENTERS)].copy()
        if len(sub) == 0:
            continue
        means = sub.groupby('Center_Num')['Area'].mean()
        if len(means) < 2:
            continue
        x = means.index.astype(float)
        y = means.values.astype(float)
        logx = np.log10(x)
        logy = np.log10(y)
        model = sm.OLS(logy, sm.add_constant(logx)).fit()
        beta = model.params[1]
        ci = model.conf_int(alpha=0.05)
        half_ci = (ci[1, 1] - ci[1, 0]) / 2
        entries.append((beta, name, color, marker, x, y, half_ci, model.rsquared))

    entries.sort(key=lambda e: e[0])
    for beta, name, color, marker, x, y, half_ci, r2 in entries:
        print(f"  {name:<10s}  beta = {beta:.3f} +/- {half_ci:.3f},  R2 = {r2:.3f}")
        ax.scatter(x, y, marker=marker, s=25, edgecolor=color,
                   linewidth=0.6, facecolor='none', zorder=10)

    x_ref = np.array([10**0.4, 10**2.2])
    ax.plot(x_ref, 10**2.5 * (x_ref / x_ref[0])**1, color='gray', lw=0.8,
            ls=(0, (8, 4)), zorder=0)

    ax.set_xlabel(r'log$_{10}$ Center number', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_ylabel(r'log$_{10}$ Area (km²)', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD + 1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if xticks is None:
        xticks = [10**0.5, 10**1.0, 10**1.5, 10**2.0]
        xticklabels = [0.5, 1.0, 1.5, 2.0]
    else:
        xticklabels = [round(np.log10(v), 2) for v in xticks]
    if yticks is None:
        yticks = [10**2, 10**3, 10**4, 10**5]
        yticklabels = [2, 3, 4, 5]
    else:
        yticklabels = [round(np.log10(v), 2) for v in yticks]
    ax.set_xticks(xticks, xticklabels)
    ax.set_yticks(yticks, yticklabels)
    utils.FormatAxis(ax)


def _plot_pop_vs_center(ax, clusters):
    """Fig. 3b: log-log scaling of center number vs population (N ~ P^gamma)."""
    clusters = clusters[clusters['Center_Num'] >= MIN_CENTERS].copy()
    means = clusters.groupby('Center_Num')['Pop'].mean()
    x_fit = means.values.astype(float)
    y_fit = means.index.values.astype(float)
    logx = np.log10(x_fit)
    logy = np.log10(y_fit)
    model = sm.OLS(logy, sm.add_constant(logx)).fit()
    slope = model.params[1]
    ci = model.conf_int(alpha=0.05)

    ax.scatter(clusters['Pop'], clusters['Center_Num'], marker='s', s=20,
               edgecolor='#91d5f2', linewidth=0.25, facecolor='none')
    ax.scatter(x_fit, y_fit, marker='^', s=28, edgecolor='#AF58BA', linewidth=0.5, facecolor='none')
    x_line = np.linspace(4.8, 7.3, 100)
    ax.plot(10**x_line, 10**model.predict(sm.add_constant(x_line)),
            color='#AF58BA', lw=0.8, zorder=10)

    half_ci = (ci[1, 1] - ci[1, 0]) / 2
    print(f"Panel (b) N~P:  beta = {slope:.3f} +/- {half_ci:.3f},  R2 = {model.rsquared:.3f}")
    
    ax.set_xlabel(r'log$_{10}$ Population', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_ylabel(r'log$_{10}$ Center number', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD + 1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(10**4.5, 10**7.5)
    ax.set_xticks([10**5, 10**6, 10**7], [5, 6, 7])
    ax.set_ylim(10**0, 10**2.5)
    ax.set_yticks([10**0, 10**0.5, 10**1.0, 10**1.5, 10**2.0, 10**2.5], [0, 0.5, 1.0, 1.5, 2.0, 2.5])
    utils.FormatAxis(ax)


def _plot_area_vs_pop(ax, clusters):
    """Fig. 3c: log-log scaling of area vs population (A ~ P^alpha)."""
    clusters = clusters.copy()
    clusters = clusters[(clusters['Area'] > 0) & (clusters['Pop'] > 0) & (clusters['Center_Num'] >= MIN_CENTERS)]
    log_pop = np.log10(clusters['Pop'].values.astype(float))
    log_area = np.log10(clusters['Area'].values.astype(float))
    model = sm.OLS(log_area, sm.add_constant(log_pop)).fit()
    slope = model.params[1]
    ci = model.conf_int(alpha=0.05)

    ax.scatter(clusters['Pop'], clusters['Area'], marker='s', s=20,
               edgecolor='#91d5f2', linewidth=0.25, facecolor='none')
    x_line = np.linspace(4.8, 7.3, 100)
    ax.plot(10 ** x_line, 10 ** model.predict(sm.add_constant(x_line)),
            color='#AF58BA', lw=0.8, zorder=10)

    half_ci = (ci[1, 1] - ci[1, 0]) / 2
    print(f"Panel (c) A~P:  beta = {slope:.3f} +/- {half_ci:.3f},  R2 = {model.rsquared:.3f}")
    
    ax.set_xlabel(r'log$_{10}$ Population', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_ylabel(r'log$_{10}$ Area (km²)', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD + 1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(10**4.5, 10**7.5)
    ax.set_xticks([10**5, 10**6, 10**7], [5, 6, 7])
    ax.set_ylim(10**1.5, 10**4)
    ax.set_yticks([10**1.5, 10**2, 10**2.5, 10**3, 10**3.5, 10**4], [1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
    utils.FormatAxis(ax)


def _plot_beta_histogram(ax, df, legend_kw=None):
    """Fig. 3d: histogram of beta and alpha/gamma across countries."""
    all_vals = np.concatenate([df['beta'].values, df['ratio'].values])
    bw = 0.125
    lo = 1.0 - bw / 2 - np.ceil((1.0 - bw / 2 - all_vals.min()) / bw) * bw
    hi = 1.0 + bw / 2 + np.ceil((all_vals.max() - 1.0 - bw / 2) / bw) * bw
    bins = np.arange(lo, hi + 0.001, bw)
    ax.hist(df['ratio'], bins=bins, alpha=0.6,
            label=r'$\alpha/ \gamma$', color='#AF58BA', edgecolor='white', linewidth=0.5)
    ax.hist(df['beta'], bins=bins, alpha=0.6,
            label=r'$\beta$', color='#91d5f2', edgecolor='white', linewidth=0.5)
    ax.axvline(1, color='gray', lw=0.8, ls='--', zorder=0)
    
    ax.set_xlim(0.68, 1.32)
    ax.set_xticks([0.75, 1.0, 1.25])
    ax.set_ylim(0, 10)
    ax.set_yticks([0, 5, 10])
    ax.set_xlabel(r'Scaling exponent', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_ylabel('Number of countries', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    if legend_kw is not None:
        ax.legend(**legend_kw)
    utils.FormatAxis(ax)


def _plot_continent_scatter(ax, df, annotate=False):
    """Fig. 3e: scatter of gamma vs alpha colored by continent."""
    df = df.copy()
    df['continent'] = df['Country'].map(ALL_COUNTRY_CONTINENT)
    texts = []
    for continent, (color, marker) in sorted(CONTINENT_STYLE.items()):
        mask = df['continent'] == continent
        if not mask.any():
            continue
        ax.scatter(df.loc[mask, 'gamma'], df.loc[mask, 'alpha'], marker=marker, s=18,
                   edgecolor=color, linewidth=0.5, facecolor='none', label=continent, zorder=10)
        if annotate:
            for _, row in df[mask].iterrows():
                short = row['Name'].split()[0]
                texts.append(ax.text(row['gamma'], row['alpha'], short,
                                     fontsize=utils.LEGEND_SIZE - 1, va='center',
                                     ha='left', color=color))
    r_lims = [0.4, 1.0]
    ax.plot(r_lims, r_lims, color='gray', lw=0.8, ls='--', zorder=0)
    ax.set_xlim(r_lims)
    ax.set_ylim(r_lims)
    ax.set_xticks([0.4, 0.6, 0.8, 1.0])
    ax.set_yticks([0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel(r'$\gamma$', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD)
    ax.set_ylabel(r'$\alpha$', fontsize=utils.LABEL_SIZE, labelpad=utils.LABEL_PAD + 1)
    if annotate and texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle='-', color='gray', lw=0.4),
                    expand=(1.3, 1.5), force_text=(0.5, 0.8))
    utils.FormatAxis(ax)


# ── Figure assembly ────────────────────────────────────────────────────────────

def final_figure(parm=4, annotate_e=False):
    """Compose Figure 3: decomposition of scaling laws (5 panels)."""
    clusters = _build_clusters(parm)

    df = _country_betas(clusters)
    df = df.sort_values('beta').reset_index(drop=True)
    print(f"{'Name':<30} {'n':>5} {'beta':>7} {'ci':>6} {'gamma':>7} {'ci':>6} {'alpha':>7} {'ci':>6} {'a/g':>7}")
    for _, row in df.iterrows():
        print(f"  {row['Name']:<28} {row['n']:>5} "
              f"{row['beta']:>7.3f} {row['ci_AN']:>6.3f} "
              f"{row['gamma']:>7.3f} {row['ci_NP']:>6.3f} "
              f"{row['alpha']:>7.3f} {row['ci_AP']:>6.3f} "
              f"{row['ratio']:>7.3f}")

    us_clusters = clusters[clusters['Country'] == 142]

    fig = plt.figure(figsize=(16 / 2.54, 10 / 2.54), dpi=600)
    gs = GridSpec(2, 4, figure=fig, wspace=0.45, hspace=0.45, width_ratios=[1, 0.6, 1, 1])

    ax_a = fig.add_subplot(gs[:, :2])
    ax_b = fig.add_subplot(gs[0, 2])
    ax_c = fig.add_subplot(gs[0, 3])
    ax_d = fig.add_subplot(gs[1, 2])
    ax_e = fig.add_subplot(gs[1, 3])

    a_lims = dict(xlim=(10**0.4, 10**2.2), xticks=[10**0.5, 10**1.0, 10**1.5, 10**2.0])
    _plot_multi_country_area_vs_center(ax_a, clusters,
                                       bg_country_ids=df['Country'].tolist(),
                                       **a_lims)
    _plot_pop_vs_center(ax_b, us_clusters)
    _plot_area_vs_pop(ax_c, us_clusters)
    _plot_beta_histogram(ax_d, df, legend_kw=dict(
        fontsize=utils.LEGEND_SIZE, loc='upper left',
        bbox_to_anchor=(0.58, 1.1), handlelength=1.3, handletextpad=0.4))
    _plot_continent_scatter(ax_e, df, annotate=annotate_e)

    ax_a.text(-0.1, 1.0, 'a', transform=ax_a.transAxes, fontsize=utils.ORDER_SIZE, fontweight='bold', va='bottom', ha='right')
    ax_b.text(-0.2, 1.0, 'b', transform=ax_b.transAxes, fontsize=utils.ORDER_SIZE, fontweight='bold', va='bottom', ha='right')
    ax_c.text(-0.2, 1.0, 'c', transform=ax_c.transAxes, fontsize=utils.ORDER_SIZE, fontweight='bold', va='bottom', ha='right')
    ax_d.text(-0.2, 1.0, 'd', transform=ax_d.transAxes, fontsize=utils.ORDER_SIZE, fontweight='bold', va='bottom', ha='right')
    ax_e.text(-0.2, 1.0, 'e', transform=ax_e.transAxes, fontsize=utils.ORDER_SIZE, fontweight='bold', va='bottom', ha='right')

    fig.subplots_adjust(top=0.96, bottom=0.1, left=0.08, right=0.98)
    suffix = '_annotated' if annotate_e else ''
    out = rf"D:\Research\Graph\GHSL\panel\decomposition.png"
    fig.savefig(out, dpi=600)
    print(f"Saved: {out}")


def robustness_figure():
    """Generate Supplementary Figure 11: robustness across minimum contour areas."""
    fig = plt.figure(figsize=(18 / 2.54, 6 / 2.54), dpi=600)
    gs = GridSpec(1, 3, figure=fig, wspace=0.35)

    for col, parm in enumerate([2, 3, 4]):
        clusters = _build_clusters(parm)
        df = _country_betas(clusters)
        df = df.sort_values('beta').reset_index(drop=True)

        ax = fig.add_subplot(gs[0, col])
        lims = dict(xlim=(10**0.4, 10**2.5),
                    xticks=[10**0.5, 10**1.0, 10**1.5, 10**2.0, 10**2.5],
                    ylim=(10**1, 10**4),
                    yticks=[10**1, 10**2, 10**3, 10**4])
        _plot_multi_country_area_vs_center(ax, clusters,
                                           bg_country_ids=df['Country'].tolist(),
                                           **lims)
        label = chr(ord('a') + col)
        ax.set_title(rf'Minimum contour area = {parm} km²', fontsize=utils.LABEL_SIZE + 1)
        ax.text(-0.15, 1.05, label, transform=ax.transAxes,
                fontsize=utils.ORDER_SIZE, fontweight='bold', va='bottom', ha='right')

    fig.subplots_adjust(top=0.88, bottom=0.15, left=0.07, right=0.98)
    out = r"D:\Research\Graph\GHSL\panel\Chain rule robustness.png"
    fig.savefig(out, dpi=600)
    print(f"Saved: {out}")


if __name__ == "__main__":
    robustness_figure()
