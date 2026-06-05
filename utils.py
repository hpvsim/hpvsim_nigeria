"""Utilities for hpvsim_nigeria: fonts and calibration datafiles."""
import sciris as sc


def set_font(size=None, font='Libertinus Sans'):
    """Set the plotting font."""
    sc.fonts(add=sc.thisdir(aspath=True) / 'assets' / 'LibertinusSans-Regular.otf')
    sc.options(font=font, fontsize=size)


def make_datafiles():
    """Calibration target files for Nigeria."""
    return [
        'data/nigeria_cancer_cases.csv',
        'data/nigeria_cin_types.csv',
        'data/nigeria_cancer_types.csv',
    ]
