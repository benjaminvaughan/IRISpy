################################################################################
# Name : nh_iris_map
# Purpose : Copy of the IDL code that Mike Zemcov wrote to talk to the Mosaique
# scripts, but ported over to python. It continues multiple potential fields
# depending on what you want the location of the origin for your new header
# information to be. This creates a mosaic with all of the IRIS maps that
# intersect with the given field
#Author : Benjamin Vaughan
#Start Date : Oct 11, 2019
#Additional Info
#
################################################################################
from mosaic import mosaic
import numpy as np
import config
from astropy.io import fits
from utility import *
from astropy.wcs import WCS as world
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import FK4

#choose the field that you want
#fieldnum is the label of the field
fieldnum = 1
#num corresponds to the field number in the fields.txt list and not related to the matlab label, used to index fields.txt
field = config.fields[:,fieldnum] # : selects both RA and Dec

#note that the fields in fields.txt are given in the ICRS coordinate frame, but the IRAS data has B1950 Astrometry.
ra_c = field[0] *u.deg
dec_c = field[1] * u.deg
coord = SkyCoord(ra=ra_c, dec=dec_c)
b1950_coord = coord.transform_to(FK4(equinox='B1950'))
ra_c_b = float(b1950_coord.ra.to_string(decimal=True))
dec_c_b = float(b1950_coord.dec.to_string(decimal=True))

pixsize = 0.0250000000000 #pixels / arcsecond of observation to interpolate IRIS
# maps to.
naxis = 22
name = 'rxj1347'
#square patch of sky, this can be changed if desired.
naxis1 = naxis2 = naxis

y = np.arange(0, naxis1)
x = np.arange(0, naxis2)
yvec = np.repeat(x[:, np.newaxis], naxis2, axis=1)
xvec = np.repeat(y[np.newaxis, :], naxis1, axis=0)
head = make_header(pixsize, naxis, ra_c_b, dec_c_b) #this function assumes a square.

w = world(head)
c = pixel_to_skycoord(xvec, yvec, w, origin=0)
c = c.transform_to('icrs')

#this is converting the pixel coords to right ascension and declination in fk4
ra = np.asarray(c.ra.to_string(decimal=True), dtype=float)
dec = np.asarray(c.dec.to_string(decimal=True), dtype=float)

#get information from the IRAS astrometry positions
map = mosaic(head, band=4)

hdu = fits.PrimaryHDU(map, head)
hdul = fits.HDUList([hdu])
hdul.writeto(config.DataDir + 'iris_0' + name + '_fx.fits')
#ra
hdu = fits.PrimaryHDU(ra, head)
hdul = fits.HDUList([hdu])
hdul.writeto(config.DataDir + 'iris_0' + name + '_ra.fits')
#dec
hdu = fits.PrimaryHDU(dec, head)
hdul = fits.HDUList([hdu])
hdul.writeto(config.DataDir + 'iris_0' + name + '_dc.fits')
