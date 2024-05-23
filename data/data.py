import os
import copy
import numpy as np
import GEOparse
import pandas as pd

class Data :
    def _get_type( self ):

        types = set()
            
        for gsm_name, gsm in self.gse.gsms.items() :
            types.add( gsm.metadata['type'][0] )

        types = list( types )

        assert len(types) == 1, "Multiple types found"

        return types[0]

    def _load_data_type_rna( self ):
        genes = None 
        sample_ids = []
        expressions = []

        for gsm_name, gsm in self.gse.gsms.items():
            g = gsm.table['ID_REF'].to_numpy()

            if genes is not None :
                assert (g == genes).all(), f"Dataset {self.geo_tag} : Gene names for different samples does not match."
            else :
                genes = g

            sample_ids.append(gsm_name)
            expressions.append( gsm.table['VALUE'].to_numpy() )

        expressions = np.array( expressions )

        return genes, sample_ids, expressions

    def _load_data_type_sra( self ):
        print("Loading SRA")

        #self.gse.download_supplementary_files(download_sra=True,directory=self.destdir)
        pass

    def _parse_data( self ):

        genes = None 
        sample_ids = []
        expressions = {}

        data_type = self.type

        if data_type == "RNA" :
            genes, sample_ids, expressions = self._load_data_type_rna()
        elif data_type == "SRA" :
            pass
            #self._load_data_type_sra()
        else :
            assert False, f"Data type {data_type} is not supported"

        if genes is not None and len(expressions) > 0 and len(sample_ids) > 0 :
            self._valid = True

            self._genes = genes
            self._valid_genes = np.arange( len(self._genes) )
            self._sample_ids = sample_ids
            self._expressions = expressions

    def __init__( self, geo_tag, pattern, cache_root ):
        self.geo_tag = geo_tag 
        self.pattern = pattern
        
        self._valid = False

        destdir = os.path.join( cache_root, pattern )

        if not os.path.isdir( destdir ):
            os.mkdir( destdir )

        self.destdir = destdir

        self.gse = GEOparse.get_GEO(geo=self.geo_tag, destdir=self.destdir, silent=True)

        self._type = None
        self._parse_data()

    @property 
    def valid( self ):
        return self._valid

    @property 
    def genes( self ):
        return self._genes[self._valid_genes]

    @property 
    def sample_ids( self ):
        return self._sample_ids

    @property 
    def type( self ):
        if self._type is None :
            self._type = self._get_type()
        return self._type

    @property 
    def expressions( self ):
        sample_ids = self.sample_ids 
        genes = self.genes 
        expressions = np.array( self._expressions ).T 
        expressions = pd.DataFrame( expressions, index=genes, columns=sample_ids )

        return expressions

    def shape( self ):
        print( self.geo_tag )
        pass

    def samples( self ):
        print( self.geo_tag )
        for key, item in self.gse.gpls.items() : 
            print( item.columns )

    def do_stuff( self ):
        print( self.pattern, self.geo_tag )
        self._parse_data()

    def reload_data( self ):
        self._parse_data()

    def print_sample_metadata( self ):

        sample_ids = self.sample_ids

        for s in sample_ids :

            print()
            print(f'## {s}')
            self.gse.gsms[s].show_metadata()



