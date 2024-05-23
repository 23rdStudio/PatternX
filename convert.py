import numpy as np 
import pandas as pd 
import pickle


class Module :
    def __init__( self ):
        with open('gene_symbol2probeid.pkl','rb') as ff :
            self.mapping = pickle.load(ff)

    def convert( self, data_name ):

        with open(f'{data_name}.pkl','rb') as ff :
            data = pickle.load( ff )

        expressions = data['expressions'] 

        index = expressions.index 

        c = 0
        count = 0 

        indices = []
        columns = expressions.columns 

        dd = []

        for gs, psids in self.mapping.items() :

            rows = []
            for p in psids :
                if p in index :
                    rows.append( expressions.loc[p].to_numpy().ravel() )

            if len( rows ) == 0 :
                continue 

            rows = np.array( rows )

            if len(rows) == 0 :
                continue

            d = np.mean( rows, axis=0 )


            dd.append( d )
            indices.append( gs )

        data_gs = pd.DataFrame( dd, index=indices, columns=columns )
        

        data['expressions_gs'] =  data_gs 

        with open(f'{data_name}.pkl','wb') as ff :
            pickle.dump(data, ff)

    def do_stuff( self ):
        geo_tags = ['GSE148346', 'GSE6281', 'GSE11792', 'GSE10433', 'GSE6475', 'GSE92566', 'GSE32887' ]
        
        for geo_tag in geo_tags :
            print( geo_tag )
            self.convert( geo_tag )
