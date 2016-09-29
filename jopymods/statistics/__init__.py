#coding:utf8

"""
Collections of functions for plotting stuff mostly using matplotlib.pyplot
"""

from ._chi_square import *
from ._Kolmogorow import *

__all__ = [_s for _s in dir() if not _s.startswith('_')]
