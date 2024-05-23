import pandas as pd 
import numpy as np
import pickle


if __name__=="__main__" :

    with open('gene_symbol2probeid.pkl','rb') as ff :
        mapping = pickle.load(ff)

    with open('GSE10433_20240401.pkl','rb') as ff :
        data = pickle.load(ff)

    expressions = data['expressions']

    sample_ids = expressions.columns.to_numpy()
    gene_symbols = []
    gene_values = []

    for key, items in mapping.items() :
    
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

    print( data )


