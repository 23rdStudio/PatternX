import os
import yaml
import numpy as np
from data import load_dataset
from data import Data

class Module :
    def load_config( self ):
        with open('config.yaml', 'r') as ff :
            cfg = yaml.safe_load(ff)
        self.cfg = cfg

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

    def examine_gene_pool( self ):

        gene_pool = np.sort( self.dataset[1].genes )

        for idx, d in enumerate(self.dataset) : 
            count = 0

            for g in d.genes :
                if g in gene_pool :
                    count = count + 1

            print( idx+1, d.pattern, d.geo_tag, " - Number of genes in the pool : ", count , " - Fraction of the genes in the pool : ", count / len(d.genes) )

    def load_dataset( self ):
        self.dataset = load_dataset( self.cfg['data'], only_valid=True )

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

if __name__=="__main__" :

    m = Module()

    m.load_dataset()
