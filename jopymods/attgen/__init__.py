#coding:utf8

"""
Collections of functions for plotting stuff mostly using matplotlib.pyplot
"""

from ._Johannes import *
from ._Patrick import *
from ._Tomasz import*

__all__ = [_s for _s in dir() if not _s.startswith('_')]
