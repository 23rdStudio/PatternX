import pickle
import numpy as np 
import pandas as pd 

from matplotlib import pyplot as pp 

dname = "GSE32887"

class Module :
    def _load_data( self ):
        
        with open(f'{dname}.pkl','rb') as ff :
            data = pickle.load(ff)

        self.data = data

    def __init__( self ):
        self._load_data()

    def do_stuff( self ):

        annots = self.data['annots']
        expressions = self.data['expressions']
        
        print( expressions.keys() )

        mat = expressions['probe_set_id'].to_numpy()

        m = np.mean( mat, axis=1 )

        pp.plot( np.sort(m) )

