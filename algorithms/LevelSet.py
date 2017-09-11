from SegAlgo import *

class LevelSet(SegAlgo):
    def __init__(self,mesh, volume):
        super(LevelSet,self).__init__(mesh,volume)

    @fn_timer
    def Optimize(self, *args, **kwargs):
        print 'Optimize of LevelSet'