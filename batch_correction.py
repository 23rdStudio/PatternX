import numpy as np 
import pandas as pd 
import anndata 
import pickle 
import omicverse as ov
import pydeseq2.preprocessing
import scanpy
import pydeseq2

def load_dataset( geotag, batch ):
    with open( f'{geotag}.pkl','rb') as ff : 
        dataset = pickle.load(ff)

    data = dataset["expressions_gs"]
    adata = anndata.AnnData( data.T )
    adata.obs['batch'] = batch

    return adata 

def load_brain_dataset( path, batch ):
    adata = scanpy.read( path )
    data = adata.layers['counts'].astype( np.float64 )

    sample_ids = list(adata.obs.index)
    gene_symbols = list(adata.var.external_gene_name)

    dframe = pd.DataFrame( data, index=sample_ids, columns=gene_symbols )

    out_adata = anndata.AnnData( dframe )
    out_adata.obs['batch'] = batch 

    return out_adata

if __name__=="__main__" :

    geotags = ['GSE6281', 'GSE6475', 'GSE10433', 'GSE11792', 'GSE32887','GSE92566','GSE148346']

    adata_list = []

    batch_idx = 0

    adata = load_brain_dataset("Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715.h5", str(batch_idx))
    adata_list.append( adata )

    for idx, geotag in enumerate(geotags) :
        batch_idx += 1
        adata = load_dataset( geotag, str(batch_idx) )
        adata_list.append( adata )

    adata = anndata.concat( adata_list, merge="same" )

    ov.bulk.batch_correction( adata, batch_key="batch" ) 

    data = adata.to_df( layer="batch_correction" )

    # Removing Negative Values
    data[ data < 0 ] = 0

    data_norm = pydeseq2.preprocessing.deseq2_norm( data )[0]

    print( data_norm )