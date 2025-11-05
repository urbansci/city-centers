# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2024-12

"""
Module storing customer plot constants and functions.
"""

import matplotlib.pyplot as plt

# COLOR
EDGE_COLOR = 'black'
BACK_COLOR = '#F8F8F8'
BAR_COLOR = '#1887CC'
GRID_COLOR = '#969696'
AREAS_COLORS = ['#E3E2E6', '#009ADE', '#F28522', '#AF58BA']  # Colors for urban area categories with different center number
AREAS_COLORS_L = ['#EFEEF2', '#009FE6', '#FF8D26', '#CB61D9']
GROUP_COLORS = ['#6741D9', '#298FCC', '#E6842E', '#24B24C', 'black']  # Colors for income groups.
COUNTRY_COLORS = ['#0000cc', '#a0527a',  '#c78f14'] # Colors for different countries

# PAD
TITLE_PAD = 8
LABEL_PAD = 3
TICK_PAD = 0.2

# FONT SIZE
TICK_SIZE = 7
LEGEND_SIZE = 7
LABEL_SIZE = 7
TEXT_SIZE = 7
ORDER_SIZE = 10
TITLE_SIZE = 7

# MARKER SIZE
MARKER_SIZE = 12
BAR_WIDTH = 0.75
LINE_WIDTH = 0.25

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = TEXT_SIZE
plt.rcParams['axes.titlesize'] = TITLE_SIZE
plt.rcParams['axes.labelsize'] = LABEL_SIZE
plt.rcParams['xtick.labelsize'] = TICK_SIZE
plt.rcParams['ytick.labelsize'] = TICK_SIZE
plt.rcParams['legend.fontsize'] = LEGEND_SIZE
plt.rcParams['axes.labelpad'] = LABEL_PAD

plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600

plt.rcParams['lines.linewidth'] = 0.5

plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['axes.grid'] = False

plt.rcParams['legend.frameon'] = False


def FormatAxis(ax, positions=['left', 'bottom'], spineColor='black', spineWidth=0.5, x_ticks=True, y_ticks=True, x_labels=True, y_labels=True, minor=False):
    """
    Format display style of axis.

    Args:
        ax(ax) : The axes to be formatted.
        positions(list) : A list of spines whose ticks and visibility are retained.
        spineColor(str) : The color of the visible spines.
        x_ticks(bool) : Whether to show x ticks.
        y_ticks(bool) : Whether to show y ticks.
        x_labels(bool) : Whether to show x ticklabels.
        y_labels(bool) : Whether to show y ticklabels.
    Returns:
        None
    """

    for pos in ['left', 'top', 'bottom', 'right']:
        if pos in positions:
            ax.spines[pos].set_linewidth(spineWidth)
            ax.spines[pos].set_color(spineColor)
        else:
            ax.spines[pos].set_linewidth(0)

    ax.tick_params(axis='x', length=2 if x_ticks else 0, width=spineWidth, colors='black', rotation=0, pad=TICK_PAD if x_ticks else 2 + TICK_PAD, labelbottom=x_labels, labelsize=TICK_SIZE)
    ax.tick_params(axis='y', length=2 if y_ticks else 0, width=spineWidth, colors='black', rotation=0, pad=TICK_PAD if x_ticks else 2 + TICK_PAD, labelleft = y_labels, labelsize=TICK_SIZE)
    if minor == False:
        ax.minorticks_off()