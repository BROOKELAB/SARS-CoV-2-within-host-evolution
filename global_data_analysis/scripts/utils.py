import numpy as np


def plot_style(grey='#333333'):
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = 'Arial'
    mpl.rcParams['font.weight'] = 'light'
    mpl.rcParams['text.color'] = grey
    mpl.rcParams['axes.labelcolor'] = grey
    mpl.rcParams['xtick.color'] = grey
    mpl.rcParams['ytick.color'] = grey
    # Font sizes
    mpl.rcParams['figure.titlesize'] = 18
    mpl.rcParams['axes.titlesize'] = 18
    mpl.rcParams['axes.labelsize'] = 18
    mpl.rcParams['xtick.labelsize'] = 18
    mpl.rcParams['ytick.labelsize'] = 18
    # Border colors
    mpl.rcParams['axes.edgecolor'] = grey

