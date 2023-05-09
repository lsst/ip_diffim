#!/usr/bin/env python3
import sys
sys.path.append('~/transinet4lsst')
from mask_tools import should_ignore, get_enabled_flags_names

def vis(diff1file, diff2file, truth):
    import sys
    import pandas as pd
    from vis_tools import df_to_fk5_tempfile

    #--- load truth list and save as a temp text file
    truth = truth[abs(truth['delta_flux'])>100]

    # remove any row of truth for which the should_ignore(mask) is True
    print(len(truth))
    truth = truth[~truth.apply(lambda row: should_ignore(row['mask']), axis=1)]
    print(len(truth))

    truthTempFile = df_to_fk5_tempfile(truth,
                                   color='green')

    #--- Run DS9
    cmd = f'ds9 -mode pan {diff1file} -scale limits -1 1 -linear {diff2file} -single -lock frame wcs -wcs degrees -regions load all {truthTempFile.name}'
    import os
    os.system(cmd)


if __name__ == '__main__':

    collection = sys.argv[1] if len(sys.argv)>1 else 'latest' #pick the latest collection in the butler registry by default

    # Create a lsst butler object and fetch the two difference images.
    from lsst.daf.butler import Butler
    butler = Butler('/repo/apv')

    #--- Get a glimpse of the butler repository

    #- get the latest collection
    if collection == 'latest':
        collection = butler.registry.queryCollections()[-1]
        print(f'Using the latest collection: {collection}')

    import sys
    sys.path.append('/home/nima/lsst-repos/ap_pipe-notebooks/tools')
    from reg_tools import list_datasets
    rr = list_datasets(butler,
                      collections=collection,
                      dataID = '',
                      datasetType='*')

    ids = []
    for r in rr:
        ids.append(r.dataId.byName())
        ids[-1]['type']=r.datasetType.name

    import pandas as pd
    df = pd.DataFrame(ids)
    print(df)


    #--- Get the difference image file URIs
    visit = 982985
    detector = 164
    diff1file = butler.getURI('goodSeeingDiff_differenceTempExp', collections=collection,
                              dataId={'instrument': 'LSSTCam-imSim', 'visit': visit, 'detector': detector}).ospath
    diff2file = butler.getURI('goodSeeingTNDiff_differenceTempExp', collections='test_transinet/20230502T193020Z', # Let's use a pre-computed AL diff image.
                              dataId={'instrument': 'LSSTCam-imSim', 'visit': visit, 'detector': detector}).ospath

    #--- Get the truth file URI (which is right next to the calexp file)
    sciFile = butler.getURI('calexp',
                            collections='ap_verify-output',
                            dataId={'instrument': 'LSSTCam-imSim', 'visit': visit, 'detector': detector}).ospath
    truthfile = sciFile + '.SNVAR.parq'
    print(truthfile)



    #--- Obtain mask info for truth rows
    truth = pd.read_parquet(truthfile, engine='pyarrow')
    sci = butler.get('calexp',
                     collections='ap_verify-output',
                     dataId={'instrument': 'LSSTCam-imSim', 'visit': visit, 'detector': detector})

    X,Y = sci.getWcs().skyToPixelArray(ra=truth['ra'],dec=truth['dec'],degrees=True)
    truth = truth.assign(x=X,y=Y) # add x,y columns to truth
    mask = sci.getMaskedImage().getMask()
    truth['mask'] = mask.getArray()[truth['y'].astype(int),truth['x'].astype(int)]


    vis(diff1file, diff2file, truth)
