
import numpy as np 
import pandas as pd 
import pickle

from matplotlib import pyplot as pp 
pp.ion()

class Module :
    def __init__( self, geo_tag ):
        self.geo_tag = geo_tag

        with open(f'{self.geo_tag}.pkl','rb') as ff :
            self.data = pickle.load(ff)

    def do_stuff( self ):


        expressions = self.data['expressions_gs']

        e = expressions.to_numpy()

        print( e.shape )

        m = np.mean( e, axis=1 )
        s = np.std( e, axis=1 )
        
        pp.close('all')
        pp.figure()
        pp.plot( np.sort( m ) )
        pp.grid()
        pp.xlabel('Sorted genes based on expression average')
        pp.ylabel('Gene Expression - Mean')
        pp.title(f'Dataset {self.geo_tag}')
        pp.savefig(f'{self.geo_tag}_mean.pdf')

        pp.figure()
        pp.plot( np.sort( s ) )
        pp.grid()
        pp.xlabel('Sorted genes based on expression std')
        pp.ylabel('Gene Expression - STD')
        pp.title(f'Dataset {self.geo_tag}')
        pp.savefig(f'{self.geo_tag}_std.pdf')



