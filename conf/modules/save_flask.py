from numpy import log, pi
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
import twopoint
from astropy.io import fits
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline

def interp_loglog(x, xp, fp, logx=True, logy=True, **kwargs):
    """
    Performs linear interpolation where x and/or y are log-scaled.
    
    Parameters
    ----------
    x : array_like
        The x-coordinates at which to evaluate the interpolated values.
    xp : 1-D sequence of floats
        The x-coordinates of the data points, must be increasing if argument period is not specified. Otherwise, xp is internally sorted after normalizing the periodic boundaries with xp = xp % period.
    fp : 1-D sequence of float or complex
        The y-coordinates of the data points, same length as xp.
    logx : bool, optional
        Wheather to log-scale x, by default True
    logy : bool, optional
        Wheather to log-scale y, by default True
    
    Returns
    -------
    y : float or complex (corresponding to fp) or ndarray
        The interpolated values, same shape as x.
    """
    if logx:
        funcx = np.log
    else:
        funcx = lambda x: x

    if logy:
        funcy_in = np.log
        funcy_out = np.exp
    else:
        funcy_in = lambda x: x
        funcy_out = lambda x: x
        
    return funcy_out(np.interp(funcx(x), funcx(xp), funcy_in(fp), **kwargs))

def setup(options):
    # Load twopoint file
    filename = options.get_string(option_section, "data_file")
    covmat_name = options.get_string(option_section, "covmat_name", default="COVMAT")
    ell_min_log = options.get_int(option_section, "ell_min_log")
    two_point_data = twopoint.TwoPointFile.from_fits(filename, covmat_name)
    hdus = fits.open(filename)

    # Get data sets
    data_sets = options.get_string(option_section, "data_sets", default="").split()

    # Load band-power window matrices
    bpws_prefix = options.get_string(option_section, 'bpws_prefix', default='bpws')
    bpws = {}
    data_file = fits.open(filename)
    print("data_file opened")
    for name in data_sets:
        s = two_point_data.get_spectrum(name)       
        bpws[name] = {}
        for b1, b2 in s.bin_pairs:
            bpws[name][b1,b2] = np.array(data_file['{}_{}_{}_{}'.format(bpws_prefix, name,b1,b2).capitalize()].data)[s.angular_bin[s.get_pair_mask(b1,b2)]]
            lmaxp1 = bpws[name][b1,b2].shape[1]
            # Add extra bpws for galaxy_cl
            if name == 'galaxy_cl':
                assert b1 == b2
                for b in range(b1, 7): # TODO: HARDCODED; LJF - right now use 6 for redmagic and 7 for maglim
                    bpws[name][b, b1] = bpws[name][b1, b1]

    ell_lims_min = hdus['ELL_LIMS'].data['ell_lims_min']
    ell_lims_max = hdus['ELL_LIMS'].data['ell_lims_max']
    ells = [np.arange(ell_low, ell_up+1).mean() for ell_low, ell_up in zip(ell_lims_min, ell_lims_max)]
    ell_edges = [ell_lims_min[0]-0.5]+list(ell_lims_max+0.5)

    return {'bpws':bpws, 'data_sets':data_sets, 'ells':ells, 'ell_edges':ell_edges, 'lmaxp1':lmaxp1, 'ell_min_log':ell_min_log}


def execute(block, config):
    bpws = config['bpws']
    data_sets = config['data_sets']
    ells = config['ells']
    ell_edges = config["ell_edges"]
    lmaxp1 = config["lmaxp1"]
    ell_min_log = config["ell_min_log"]

    # ell_interp = np.arange(ell_min_log,lmaxp1)
    # TODO: HARDCODED
    ell_interp = np.arange(ell_min_log,20000)

    for name in data_sets:
        # Test that pre-binning Cl are computed at all integer ells
        # np.testing.assert_array_almost_equal(block[name, 'ell'], np.arange(len(block[name, 'ell'])))
        # Test that pre-binning Cl are computed at all integer ells up to ell_min_log
        np.testing.assert_array_almost_equal(block[name, 'ell'][:ell_min_log], np.arange(ell_min_log))

        # Rename section
        block._copy_section(name, name+'_pre_binned_bpws')
        block[name, 'bin_avg'] = True
        for b1,b2 in bpws[name].keys():
            a = 'bin_{}_{}'.format(b1,b2)
            cl_in = block[name+'_pre_binned_bpws', a]
            clout = np.concatenate([cl_in[:ell_min_log], interp_loglog(ell_interp, block[name, 'ell'][ell_min_log:], cl_in[ell_min_log:], logx=True, logy=np.all(cl_in[ell_min_log:]>0.))])
            tmp = np.array([np.arange(len(clout))[1:], clout[1:]]).T
            if name == 'galaxy_cl':
                ofn = f'Y3_5x2pt_Nsource4-Cl_f{b2}z{b2}f{b1}z{b1}.dat'
            elif name == 'galaxy_shear_cl':
                ofn = f'Y3_5x2pt_Nsource4-Cl_f10z{b2}f{b1}z{b1}.dat'
            elif name == 'shear_cl':
                ofn = f'Y3_5x2pt_Nsource4-Cl_f10z{b2}f10z{b1}.dat'
            np.savetxt('Cl_flaskv2p0_nolimber_emu_Nsource4/' + ofn, tmp, fmt=('%d', '%e'))
            # block[name, a] = np.dot(bpws[name][b1,b2], clout)

        # Save lens n(z)'s for FLASK
        # TODO: HARDCODED normalizations here based on example3x2
        #       If based on number densities, that should be an input parameter
        # norm = [2.409e-02, 4.235e-02, 6.855e-02, 3.505e-02, 3.469e-02, 3.4e-2]
        # Those normalizations below are based on 2x2 config space papers (https://arxiv.org/pdf/2105.13546.pdf)
        # REDMAGIC:
        # norm = [0.022, 0.038, 0.058, 0.029, 0.025]
        # MAGLIM:
        norm = [0.150, 0.107, 0.109, 0.146, 0.106, 0.100]
        for b in range(1, 7): # LJF - 6 for redmagic, 7 for maglim
            norm_tmp = spline(block['nz_lens', 'z'], block['nz_lens', f'bin_{b}'], ext='zeros').integral(0, 5) 
            tmp = np.array([block['nz_lens', 'z'], block['nz_lens', f'bin_{b}'] * norm[b-1] / norm_tmp]).T
            np.savetxt(f'nl_3x2v2p0_f{b}.dat', tmp, fmt='%e')

        block[name, 'ell'] = ells
        block[name, 'ell_edges'] = ell_edges

    # signal that everything went fine
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
