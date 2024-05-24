import numpy as np 
import pandas as pd 
import anndata 
import pickle 
import omicverse as ov

def load_dataset( geotag, batch ):
    with open( f'{geotag}.pkl','rb') as ff : 
        dataset = pickle.load(ff)

    data = dataset["expressions_gs"]
    adata = anndata.AnnData( data.T )
    adata.obs['batch'] = batch

    return adata 

if __name__=="__main__" :

    geotags = ['GSE6281', 'GSE6475', 'GSE10433', 'GSE11792', 'GSE32887','GSE92566','GSE148346']

    adata_list = []

    for idx, geotag in enumerate(geotags) :
        adata = load_dataset( geotag, str(idx) )
        adata_list.append( adata )

    adata = anndata.concat( adata_list, merge="same" )

    ov.bulk.batch_correction( adata, batch_key="batch" ) 

    data = adata.layers['batch_correction']

    print( data )

    print( adata.X )