import numpy as np
import pandas as pd
import os

import eazy.hdf5
from eazy import photoz

import argparse

parser = argparse.ArgumentParser(
                    prog='Eazy photoz',
                    description='Gets redshifts',
                    epilog='Hello darkness my old friend')

parser.add_argument('--fname', default='smiles', type=str)  
parser.add_argument('--mock', default=False, type=bool)
parser.add_argument('--translate_file', default="ALL", type=str)
parser.add_argument('--templates_file', default="templates/LarsonSEDTemplates/tweak_fsps_QSF_12_v3_newtemplates.param", type=str)
parser.add_argument('--fname_append', default = "", type = str)

# sfhz templates templates/sfhz/agn_blue_sfhz_13.param

def eazy_params(main_output_directory, main_output_file, catalogue_fname, translate_file, zPDF_fname, 
                templates_file = 'templates/LarsonSEDTemplates/tweak_fsps_QSF_12_v3_newtemplates.param'):
    params = {}

    params['OUTPUT_DIRECTORY'] = main_output_directory
    params['MAIN_OUTPUT_FILE'] = main_output_file
    params['CATALOG_FILE'] = catalogue_fname
    params['TEMPLATES_COMBOS'] = -1
    params['CAT_HAS_EXTCORR'] = True
    params['Z_STEP'] = 0.01
    params['Z_MIN'] = 0.01
    params['Z_MAX'] = 40
    params['APPLY_PRIOR'] = 'n'
    params['APPLY_IGM'] = 'y'
    params['PRIOR_ABZP'] = 8.9
    params['TEMPLATES_FILE'] = templates_file
    # params['TEMPLATES_FILE'] = 'templates/fsps_full/tweak_fsps_QSF_12_v3.param'
    params['TEMP_ERR_FILE'] = 'templates/TEMPLATE_ERROR.eazy_v1.0'
    params['TEMP_ERR_A2'] = 0.0
    params['SYS_ERR'] = 0.0
    params['NOT_OBS_THRESHOLD'] = -90
    #params['N_MIN_COLORS'] = 3
    params['FILTERS_RES'] = 'FILTER.RES.latest'

    return [params, translate_file, zPDF_fname]

def run_eazy(params):
    translate_file = params[1]
    self = photoz.PhotoZ(
        translate_file = translate_file, 
        params = params[0], 
        load_prior = False, 
        n_proc=8
    )
    self.fit_catalog(self.idx, n_proc=8, prior = False, beta_prior = False)
    zout, hdu = self.standard_output()

    ## Save the zPDF
    zPDF = {}
    zPDF['zgrid'] = self.zgrid
    zPDF_fname = params[2]
    for i in range(self.lnp.shape[0]):
        zPDF['gal'+str(i+1)] = np.exp(self.lnp[i,:])
    zPDF = pd.DataFrame(zPDF)
    zPDF.to_csv(zPDF_fname)

    eazy.hdf5.write_hdf5(self, h5file=self.param['MAIN_OUTPUT_FILE'] + '.h5')

    return self

def init_from_h5(params):
    new = eazy.hdf5.initialize_from_hdf5(h5file=params['MAIN_OUTPUT_FILE'] + '.h5')

if __name__ == '__main__':

    args = vars(parser.parse_args())
    fname = args['fname']
    fname_append = args['fname_append']
    templates_file = args['templates_file']

    print(args)

    main_dir = "/Users/jordan/Documents/Epochs/"
    if args['mock']:
        main_output_directory = os.path.join(main_dir, "source_inject", f"{fname}", "eazy")
        os.makedirs(main_output_directory, exist_ok=True)

        catalogue_fname = os.path.join(main_dir, "source_inject", f"{fname}", f"{fname}_mock_eazy_photometry.csv")

        main_output_file = f"{main_output_directory}/{fname}_zphot{fname_append}.out"
        zPDF_fname = f"{main_output_directory}/{fname}_mock_photometry_eazy_zPDF{fname_append}.csv"
        translate_file = os.path.join(main_dir, "source_inject", f"eazy_zphot_mock.translate")
    else:
        main_output_directory = os.path.join(main_dir, "data/eazy", f"{fname}{fname_append}")
        os.makedirs(main_output_directory, exist_ok=True)

        catalogue_fname = os.path.join(main_dir, "data/eazy", f"{fname}_phot_eazy.csv")

        if args['translate_file'] == "ALL":
            translate_file = os.path.join(main_dir, "data/eazy", f"eazy_{fname}_zphot.translate")
        elif args['translate_file'] == "NIRCAM":
            translate_file = os.path.join(main_dir, "data/eazy",f"eazy_{fname}_zphot_nircam.translate")

        main_output_file = f"{main_output_directory}/{fname}_zphot{fname_append}.out"
        zPDF_fname = f"{main_output_directory}/{fname}_eazy_zPDF{fname_append}.csv"

    print(main_output_file)
    print(templates_file)
    input('Continue')

    params = eazy_params(
        main_output_directory,
        main_output_file,
        catalogue_fname,
        translate_file,
        zPDF_fname,
        templates_file
    )

    run_eazy(
        params
    )
