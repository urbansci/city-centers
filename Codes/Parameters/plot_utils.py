# -*- coding:utf-8 -*-
# @author  : Shuai Pang
# @time    : 2024-12

"""
Module storing customer plot constants and functions.
"""

# COLOR
EDGE_COLOR = 'black'
BACK_COLOR = '#F8F8F8'
BAR_COLOR = '#1887CC'
GRID_COLOR = '#969696'
AREAS_COLORS = ['#E3E2E6', '#009ADE', '#F28522', '#AF58BA']  # Colors for urban area categories with different center number
GROUP_COLORS = ['#6741D9', '#298FCC', '#E6842E', '#24B24C', 'black']  # Colors for income groups.
COUNTRY_COLORS = ['#0000cc', '#a0527a', '#5ab1b2', '#c78f14'] # Colors for different countries

# PAD
TITLE_PAD = 8
LABEL_PAD = 3
TICK_PAD = 0.2

# FONT SIZE
TICK_SIZE = 6
LEGEND_SIZE = 6
LABEL_SIZE = 6
TEXT_SIZE = 6
ORDER_SIZE = 10
TITLE_SIZE = 7

# MARKER SIZE
MARKER_SIZE = 12
BAR_WIDTH = 0.75
LINE_WIDTH = 0.25


def FormatAxis(ax, positions=['left', 'bottom'], spineColor='black', x_ticks=True, y_ticks=True, x_labels=True, y_labels=True):
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
            ax.spines[pos].set_linewidth(0.5)
            ax.spines[pos].set_color(spineColor)
        else:
            ax.spines[pos].set_linewidth(0)

    ax.tick_params(axis='x', length=2 if x_ticks else 0, width=0.5, colors='black', pad=TICK_PAD if x_ticks else 2 + TICK_PAD, labelbottom=x_labels, labelsize=TICK_SIZE)
    ax.tick_params(axis='y', length=2 if y_ticks else 0, width=0.5, colors='black', pad=TICK_PAD if x_ticks else 2 + TICK_PAD, labelleft = y_labels, labelsize=TICK_SIZE)