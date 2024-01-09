from .data import Data 
from tqdm import tqdm

def load_dataset( cfg, only_valid=False ):
    cache_path = cfg['cache']

    dataset = []
 
    for item in tqdm(cfg['dataset']) :
        should_load = item.get('should_load',True)

        if not should_load :
            continue

        d = Data( item['geo_tag'], item['pattern'], cache_path )

        if d.valid or not only_valid :
            dataset.append(d)

    return dataset
