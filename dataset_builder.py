import os 
import platform
import pickle
import scanpy
import anndata 
import numpy as np
import pandas as pd
import omicverse as ov


if platform.system() == "Darwin" :
    data_root = "/Users/heydar/Work/void/data/bio/brain"
elif platform.system() == "Linux" :
    data_root = ""

brain_dataset = "Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715.h5"
geotags = [ ('GSE6281','2a'), ('GSE6475','3'), ('GSE10433','3'), 
            ('GSE11792','3'), ('GSE32887','4b'),('GSE92566','4a'), 
            ('GSE148346','1') ]

class DatasetBuilder :
    def _load_brain( self ):
        db_path = os.path.join( data_root, brain_dataset ) 
        adata   = scanpy.read( db_path )

        index = adata.to_df().index
        pattern = adata.obs["Pattern"].to_numpy()
        annots = np.zeros( len(index), dtype=int )

        aa = pd.DataFrame( index=index )
        aa['annots'] = annots 
        aa['pattern'] = pattern 
        aa['dataset'] = 'brain'

        out = {}
        out['norm'] = adata.to_df()
        out['raw'] = adata.to_df(layer="counts")
        out['annots'] = aa 

        return out

    def _load_dataset( self, geotag, pattern ):
        with open( f'{geotag}.pkl','rb') as ff : 
            dataset = pickle.load(ff)

        annots = dataset['annots']
        annots['pattern'] = pattern
        annots['annots'] = annots['annots'].to_numpy().astype(int)
        annots['dataset'] = geotag

        out = {}
        out['raw'] = dataset['expressions_gs'].T
        out['annots'] = dataset['annots']

        return out

    def _load_data( self ):
        datasets = {}
        datasets['brain'] = self._load_brain()

        for geotag, pattern in geotags :
            datasets[geotag] = self._load_dataset( geotag, pattern )

        return datasets

    def _get_common_genes( self ):

        adata_list = []
        for name, dataset in self.datasets.items() :
            adata = anndata.AnnData( dataset['raw'] )
            adata.obs['batch'] = name
            adata_list.append(adata)

        adata = anndata.concat( adata_list, merge='same' )

        data = adata.to_df()

        return data.columns

    def get_adata( self, name ):
        # brain-raw : raw 
        # brain-norm : norm 
        # patterx : brain + datasets

        if name == "brain-raw" :
            return anndata.AnnData( self.datasets['brain']['raw'] ), self.datasets['brain']['annots']
        if name == "brain-norm" :
            return anndata.AnnData( self.datasets['brain']['norm'] ), self.datasets['brain']['annots']
        if name == "patternx" :
            adata_list = []
            annots_list = []
            for name, dataset in self.datasets.items() :
                adata = anndata.AnnData( dataset['raw'] )
                adata.obs['batch'] = name
                adata_list.append(adata)
                annots_list.append( dataset['annots'] )

            adata = anndata.concat( adata_list, merge='same' )
            annots = pd.concat( annots_list )

            return adata, annots

    def __init__( self ):
        self.datasets = self._load_data()
        self.common_genes = self._get_common_genes()

    def do_stuff( self ):

        adata, annots = self.get_adata( 'patternx' )
        adata, annots = self.get_adata( 'brain-norm' )

        print( annots )

