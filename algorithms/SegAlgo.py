import copy
import os
import numpy as np
import vispy
from utilities import *
from utilities.utilities import *

class SegAlgo(object):
    def __init__(self,mesh,volume):
        self._mesh = mesh
        self._volume = volume

    def Optimize(self,*args,**kwargs):
        pass