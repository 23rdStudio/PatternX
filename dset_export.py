import os
import yaml
import pickle
import pandas as pd 
import numpy as np 

from matplotlib import pyplot as pp 
pp.ion()


class Module :

    def load_data( self, geo_tag ):

        data_path = f'{geo_tag}.pkl'

        if not os.path.exists( data_path ) :
            return None 

        with open(data_path,'rb') as ff :
            data = pickle.load(ff)

        return data 

    def __init__( self ):

        with open('config.yaml','r') as ff :
            cfg = yaml.safe_load(ff)

        dataset = cfg["data"]["dataset"]

        self.data = []

        for d in dataset :
            geo_tag = d['geo_tag']
            pattern = d['pattern']

            data = self.load_data( d["geo_tag"] )

            if data is None :
                continue 

            d = {}
            d['geo_tag'] = geo_tag 
            d['pattern'] = pattern 
            d['data'] = data 

            self.data.append( d )


        print( len(self.data ))

    def do_stuff( self ):

        patterns = ['1','2a','3']

        genes_set = None 

        dd = []

        for d in self.data :
            pattern = d['pattern']

            
            if not pattern in patterns :
                continue

            data = d['data']['expressions_gs']
            annots = d['data']['annots']['annots']

            columns = data.columns.to_numpy()
            gset = set(data.index.to_list())

            if genes_set is None :
                genes_set = gset 
            else :
                genes_set = genes_set.intersection( gset )

            annots_arr = annots.to_numpy().ravel().astype( int )

            keep = np.where( annots_arr == 0 )[0]

            columns = columns[keep]
            data = data[columns]

            dd.append( data )

        selected_genes = list(genes_set)

        dd_cat = []

        for d in dd : 
            dd_cat.append( d.loc[selected_genes,:] )

        data = pd.concat( dd_cat, axis=1 )

        pp.matshow( data )

        print( data.shape )





        
