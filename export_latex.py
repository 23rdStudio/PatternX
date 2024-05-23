import pickle 

if __name__=="__main__" :
    with open('GSE6281.pkl','rb') as ff :
        data = pickle.load(ff)

    data = data['expressions_gs']

    print( data )

if __name__=="__main__2" :


    columns = data.columns.to_numpy()
    index = data.index.to_numpy()

    data = data[ columns[[0,1,2,3,4,5]] ].iloc[ [0,1,2,3,4,5,6,7,8,9] ]

    res = data.style.to_latex()
    print( res )


