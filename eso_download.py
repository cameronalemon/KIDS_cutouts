import os 
import sys

import numpy as np

from astropy import table
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
import astropy.io.fits as fits

from pyvo.dal import tap
import matplotlib.pyplot as plt


ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
tapobs = tap.TAPService(ESO_TAP_OBS)

def ESO_query(ra, dec, table='KIDS'):
    """
    Function to query the ObsCore table at ESO

    parameter obs_collection: see here for a list: 
    parameter ra: right ascension in degrees
    parameter dec: declination in degrees

    return: table of images that contain the given ra, dec in the specified table
    """

    query = """SELECT *
    FROM ivoa.ObsCore
    WHERE CONTAINS(point('', """+str(ra)+""", """+str(dec)+"""), s_region)=1
    AND dataproduct_type in ('image')
    AND obs_collection = '"""+table+"""'
    """

    res = tapobs.search(query=query, maxrec=1000)
    return res.to_table()


def construct_eso_download_html(dp_id, ra, dec, radius):
    """
    Returns the html link from the ESO dataportal
    """
    dp_id_html = dp_id.replace(':', '%3A')
    radius = "{:.5f}".format(radius)

    html_link = 'https://dataportal.eso.org/dataportal_new/soda/sync?ID='
    html_link += dp_id_html
    html_link += '&PREFIX='
    html_link += str(ra)
    html_link += '_'
    html_link += str(dec)
    html_link += '&CIRCLE='
    html_link += str(ra)
    html_link += '+'
    html_link += str(dec)
    html_link += '+'
    html_link += str(radius)

    return html_link


def download_cutouts(ra, dec, radius, table='KIDS', bands=['g_SDSS', 'r_SDSS', 'i_SDSS'], filename_prefix=None, savedir='./'):
    
    res = ESO_query(ra=ra, dec=dec, table=table)
    bands_exist = [band in res['filter'] for band in bands]
    
    #if not all the bands requested exist, then raise an Exception
    if False in bands_exist:
        print('The following bands did not exist:', list(np.array(bands)[np.array([not x for x in bands_exist])]))
        return Exception

    if not filename_prefix:
        filename_prefix = table+'_RA_'+str(ra)+'_DEC_'+str(dec)

    #download the image from each band
    filenames = []
    for band in bands:
        #pass the correct filename to the link constructor
        index = np.where(res['filter']==band)[0]
        dp_id = res['dp_id'][index][0]

        #construct the html
        html_link = construct_eso_download_html(dp_id, ra, dec, radius)

        #download the cutouts
        filename = savedir+filename_prefix+'_'+band+'.fits'
        os.system('wget -O '+filename+' "'+html_link+'"')
        filenames.append(filename)

    return filenames


def make_gri(savename, g_file, r_file, i_file):
    g = fits.open(g_file)[0].data
    r = fits.open(r_file)[0].data
    i = fits.open(i_file)[0].data


    fig, ax = plt.subplots(figsize=(5, 5))
    ax.imshow(g, interpolation='nearest', origin='lower')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.margins(x=0)

    plt.savefig(savename, bbox_inches='tight')
    plt.close()



if __name__ == '__main__':
    ra, dec = 184.29040, -2.93950
    radius = 3./3600.

    files = download_cutouts(ra=ra, dec=dec, radius=radius, filename_prefix='test')

    gri_file = make_gri(savename='test_KIDS.png', g_file=files[0], r_file=files[1], i_file=files[2])
