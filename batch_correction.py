import os
import numpy as np 
import datetime
import pandas as pd 
import anndata 
import pickle 
#import omicverse as ov
#import pydeseq2.preprocessing
import scanpy
#import pydeseq2
import platform 

from matplotlib import pyplot as pp

if platform.system() == "Darwin" :
    data_root = "/Users/heydar/Work/void/data/bio/brain"
elif platform.system() == "Linux" :
    data_root = ""

brain_dataset = "Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715.h5"

def plot_gene_count( gene_counts ):

    X = []
    Y = []

    for g in gene_counts :
        X.append( g[1] )
        Y.append( g[0] )

    pp.plot( X,Y, '.-' )
    pp.grid()
    pp.xlabel('Datasets')
    pp.ylabel('Number of genes')
    pp.xticks(rotation=45)
    pp.tight_layout()
    pp.savefig('plot.pdf')
    

def load_dataset( geotag, batch ):
    with open( f'{geotag}.pkl','rb') as ff : 
        dataset = pickle.load(ff)

    annots = dataset['annots'].to_numpy()
    pattern = np.full_like( annots, dataset['pattern'] )
    dset = np.full_like( annots, geotag )

    info = np.concatenate( [ dset, pattern, annots ], axis=1 )
    info = pd.DataFrame( info, columns=['dataset','pattern','annots'], 
                         index=dataset['annots'].index )

    data = dataset["expressions_gs"]
    adata = anndata.AnnData( data.T )
    adata.obs['batch'] = batch

    return adata, info  

def load_brain_dataset( path, batch ):
    adata = scanpy.read( path )
    data = adata.layers['counts'].astype( np.float64 )

    sample_ids = list(adata.obs.index)
    gene_symbols = list(adata.var.external_gene_name)

    pattern = list(adata.obs.Pattern)
    pattern = np.array( pattern ).reshape((-1,1))

    dset = np.full_like( pattern, 'brain' )

    annots = np.concatenate( [ dset, pattern, np.zeros(pattern.shape) ], axis=1 )
    annots = pd.DataFrame( annots, columns=['dataset','pattern','annots'], index=sample_ids )

    dframe = pd.DataFrame( data, index=sample_ids, columns=gene_symbols )

    out_adata = anndata.AnnData( dframe )
    out_adata.obs['batch'] = batch 

    return out_adata, annots

if __name__=="__main__" :

    geo_names = ['brain', 'GSE6281', 'GSE6475', 'GSE10433', 'GSE11792', 'GSE32887','GSE92566','GSE148346']

    geotags = ['GSE6281', 'GSE6475', 'GSE10433', 'GSE11792', 'GSE32887','GSE92566','GSE148346']

    adata_list = []
    annots_list = []

    batch_idx = 0

    data_path = os.path.join( data_root, brain_dataset )

    adata, annots = load_brain_dataset( data_path, str(batch_idx))
    adata_list.append( adata )
    annots_list.append( annots )

    idy = 0
    gene_counts = []
    gene_counts.append( (adata.shape[1], geo_names[idy]) )

    for idx, geotag in enumerate(geotags) :
        idy += 1

        batch_idx += 1
        adata, annots = load_dataset( geotag, str(batch_idx) )
        adata_list.append( adata )
        annots_list.append( annots )

        aa = anndata.concat( adata_list, merge="same" )
        gene_counts.append( (aa.shape[1], geo_names[idy]) )

    plot_gene_count( gene_counts )

if __name__=="__main__2" :

    adata = anndata.concat( adata_list, merge="same" )

    annots = pd.concat( annots_list, axis=0 )

    ov.bulk.batch_correction( adata, batch_key="batch" ) 

    data = adata.to_df( layer="batch_correction" )

    # Removing Negative Values
    data[ data < 0 ] = 0

    data_norm = pydeseq2.preprocessing.deseq2_norm( data )[0]

    out = {}
    out['data'] = data 
    out['data_norm'] = data_norm 
    out['annots'] = annots 

    tag = datetime.datetime.now().strftime('%Y%m%d')

    with open(f'dataset_{tag}.pkl','wb') as ff :
        pickle.dump( out, ff )
