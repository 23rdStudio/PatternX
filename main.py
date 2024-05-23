import os
import sys
import yaml
import numpy as np
from data import load_dataset
from data import Data
import datetime
import pickle
import pandas as pd

class Module :
    def load_config( self ):
        with open('config.yaml', 'r') as ff :
            cfg = yaml.safe_load(ff)
        self.cfg = cfg

        with open('gene_symbol2probeid.pkl','rb') as ff :
            self.gene_mapping = pickle.load(ff)

    def __init__( self ):
        self.load_config()

    def examine_type( self, ignore_sra=False ):
        cache_root = self.cfg['data']['cache']

        count = 0

        for idx, item in enumerate(self.cfg['data']['dataset']) :
            d = Data( item['geo_tag'], item['pattern'], cache_root=cache_root )

            if ignore_sra and d.type == "SRA" :
                continue 

            print( count+1, item['geo_tag'], item['pattern'], d.type )
            print( "Data shape [ num samples, num genes ] :", d.expressions.shape )
            count += 1

    def examine_genes( self, ignore_sra=False ):
        cache_root = self.cfg['data']['cache']

        count = 0

        for idx, item in enumerate(self.cfg['data']['dataset']) :
            d = Data( item['geo_tag'], item['pattern'], cache_root=cache_root )

            if ignore_sra and d.type == "SRA" :
                continue 

            genes = d.genes
            print( count+1, item['geo_tag'], item['pattern'], d.type )
            print( f"Number of genes : { len(genes) }")

            print( np.sort(genes) )

            print( " ****************** ") 

            count += 1

    def map2gene_symbol( self, expressions ):
        sample_ids = expressions.columns.to_numpy()
        gene_symbols = []
        gene_values = []

        for key, items in self.gene_mapping.items() :
    
            count = 0
            values = np.zeros( expressions.shape[1] )

            for i in items :
                if i in expressions.index :
                    count += 1
                    vals = expressions.loc[[i]].to_numpy().ravel()
                    values += vals 

            if count > 0 :
                values /= vals 

                gene_symbols.append( key )
                gene_values.append( values )

        data = pd.DataFrame( gene_values, index=gene_symbols, columns=sample_ids )

        return data



    def examine_gene_pool( self ):

        gene_pool = np.sort( self.dataset[1].genes )

        for idx, d in enumerate(self.dataset) : 
            count = 0

            for g in d.genes :
                if g in gene_pool :
                    count = count + 1

            print( idx+1, d.pattern, d.geo_tag, 
                    " - Number of genes in the pool : ", count , 
                    " - Fraction of the genes in the pool : ", count / len(d.genes) )

    def load_dataset( self ):
        geo_tags = ['GSE148346', 'GSE6281', 'GSE11792', 'GSE10433', 'GSE6475', 'GSE92566', 'GSE32887' ]
    


        self.dataset = load_dataset( self.cfg['data'], only_valid=True, geo_tags=geo_tags )

    def idetify_patterns( self ):

        pattern_map = {}

        for d in self.cfg['data']['dataset'] :
            pattern_map[ d['geo_tag'] ]= d['pattern']

        geo_tags = ['GSE148346', 'GSE6281', 'GSE11792', 'GSE10433', 'GSE6475', 'GSE92566', 'GSE32887' ]

        return pattern_map

    def load_annots( self, geo_tag ):
    
        out = []

        with open( f'{geo_tag}_annots.md', 'r' ) as ff :
            lines = ff.readlines()

            gsms = []
            annots = []

            values = []

            for line in lines :
                values.append( line.strip()[3:] )

            for gsm, annot in zip( values[::2], values[1::2] ):
                gsms.append( gsm )
                annots.append( annot )
    
        return pd.DataFrame(data=annots, index=gsms, columns=['annots'])
        
        
    def export( self ):

        pattern_map = self.idetify_patterns()

        tag = datetime.datetime.now().strftime('%Y%m%d')


        for d in self.dataset :
            geo_tag = d.geo_tag

            print( geo_tag )

            pattern = pattern_map[geo_tag]
            annots = self.load_annots( geo_tag )
            expressions = d.expressions

            data = {}
            data['pattern'] = pattern 
            data['annots'] = annots 
            data['expressions'] = expressions 

            with open(f'{geo_tag}.pkl','wb') as ff :
                pickle.dump( data, ff )

    def do_stuff( self ):
        gene_pool = []

        for d in self.dataset :
            if len(d.genes) > len(gene_pool) :
                gene_pool = d.genes

        gene_pool = np.sort( gene_pool )

        print( gene_pool )

        for d in self.dataset :
            
            count = 0

            for g in d.genes :
                if g in gene_pool :
                    count = count + 1

            print( d.pattern, d.geo_tag, count, count / len(d.genes) )

    def do_stuff2( self, idx ):
        d = self.dataset[idx]

        return d 

if __name__=="__main__" :
    m = Module()
    m.load_dataset()

    print( f'# {m.dataset[int(sys.argv[1])].geo_tag}')

    m.dataset[int(sys.argv[1])].print_sample_metadata()
