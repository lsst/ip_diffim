#!/usr/bin/env python3

def vis(diff1file, diff2file, truthfile):
    import sys
    import pandas as pd
    from vis_tools import df_to_fk5_tempfile

    #--- load truth list and save as a temp text file
    truth = pd.read_parquet(truthfile, engine='pyarrow')
    truth = truth[abs(truth['delta_flux'])>100]

    truthTempFile = df_to_fk5_tempfile(truth,
                                   color='green')

    #--- Run DS9                                                                                                                                                                                                                                              
    cmd = f'ds9 -mode pan {diff1file} {diff2file} -single -lock scale -lock frame wcs -wcs degrees -scale limits -1 1 -linear -regions load all {truthTempFile.name}'
    import os
    os.system(cmd)


if __name__ == '__main__':

    # Create a lsst butler object and fetch the two difference images.
    from lsst.daf.butler import Butler
    butler = Butler('/repo/apv', collections = 'test_transinet/20230313T083548Z')

    #--- Get a glimpse of the butler repository
    import sys
    sys.path.append('/home/nima/lsst-repos/ap_pipe-notebooks/tools')
    from reg_tools import list_datasets
    rr = list_datasets(butler,
                      collections=butler.collections,
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
    diff1file = butler.getURI('goodSeeingDiff_differenceTempExp', collections='test_transinet/20230309T212156Z',
                              dataId={'instrument': 'LSSTCam-imSim', 'visit': 943296, 'detector': 168}).ospath
    diff2file = butler.getURI('goodSeeingTNDiff_differenceTempExp',
                              dataId={'instrument': 'LSSTCam-imSim', 'visit': 943296, 'detector': 168}).ospath

    #--- Get the truth file URI (which is right next to the calexp file)
    sciFile = butler.getURI('calexp',
                            collections='ap_verify-output',
                            dataId={'instrument': 'LSSTCam-imSim', 'visit': 943296, 'detector': 168})
    truthfile = sciFile.ospath + '.SNVAR.parq'
    print(truthfile)

    vis(diff1file, diff2file, truthfile)

