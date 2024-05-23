from .data import Data 
from tqdm import tqdm

def load_dataset( cfg, only_valid=False, geo_tags=None ):
    cache_path = cfg['cache']

    dataset = []
 
    for item in cfg['dataset'] :
        should_load = item.get('should_load',True)

        if not should_load :
            continue

        if geo_tags is not None :
            if not item['geo_tag'] in geo_tags :
                continue

        d = Data( item['geo_tag'], item['pattern'], cache_path )

        if d.valid or not only_valid :
            dataset.append(d)

    return dataset
