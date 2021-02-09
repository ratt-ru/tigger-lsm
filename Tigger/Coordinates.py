# -*- coding: utf-8 -*-
#
# % $Id$
#
#
# Copyright (C) 2002-2011
# The MeqTree Foundation &
# ASTRON (Netherlands Foundation for Research in Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>,
# or write to the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
from __future__ import print_function, division, absolute_import
import Tigger
from Tigger import startup_dprint

startup_dprint(1, "start of Coordinates")

import math
import numpy
import traceback
import warnings
import numpy as np
from numpy import sin, cos, arcsin

startup_dprint(1, "imported numpy")

from astropy.io import fits as pyfits

startup_dprint(1, "imported pyfits")

DEG = math.pi / 180

startup_dprint(1, "importing WCS")

# If we're being imported outside the main app (e.g. a script is trying to read a Tigger model,
# whether TDL or otherwise), then pylab may be needed by that script for decent God-fearing
# purposes. Since WCS is going to pull it in anyway, we try to import it here, and if that
# fails, replace it by dummies.
if not Tigger.matplotlib_nuked:
    try:
        import pylab
    except:
        Tigger.nuke_matplotlib()

# some locales cause WCS to complain that "." is not the decimal separator, so reset it to "C"
import locale

locale.setlocale(locale.LC_NUMERIC, 'C')

from astropy.wcs import WCS, FITSFixedWarning
import PyWCSTools.wcs

startup_dprint(1, "imported WCS")
warnings.simplefilter('ignore', category=FITSFixedWarning)

def angular_dist_pos_angle(ra1, dec1, ra2, dec2):
    """Computes the angular distance between the two points on a sphere, and
    the position angle (North through East) of the direction from 1 to 2."""
    # I lifted this somewhere
    sind1, sind2 = sin(dec1), sin(dec2)
    cosd1, cosd2 = cos(dec1), cos(dec2)
    cosra, sinra = cos(ra1 - ra2), sin(ra1 - ra2)

    adist = numpy.arccos(min(sind1 * sind2 + cosd1 * cosd2 * cosra, 1))
    pa = numpy.arctan2(-cosd2 * sinra, -cosd2 * sind1 * cosra + sind2 * cosd1)
    return adist, pa


def angular_dist_pos_angle2(ra1, dec1, ra2, dec2):
    """Computes the angular distance between the two points on a sphere, and
    the position angle (North through East) of the direction from 1 to 2."""
    # I re-derived this from Euler angles, but it seems to be identical to the above
    ra = ra2 - ra1
    sind0, sind, cosd0, cosd = sin(dec1), sin(dec2), cos(dec1), cos(dec2)
    sina, cosa = sin(ra) * cosd, cos(ra) * cosd
    x = cosa * sind0 - sind * cosd0
    y = sina
    z = cosa * cosd0 + sind * sind0
    print(x, y, z)
    PA = numpy.arctan2(y, -x)
    R = numpy.arccos(z)

    return R, PA


def angular_dist_pos_angle2(ra1, dec1, ra2, dec2):
    """Computes the angular distance between the two points on a sphere, and
    the position angle (North through East) of the direction from 1 to 2."""
    # I re-derived this from Euler angles, but it seems to be identical to the above
    ra = ra2 - ra1
    sind0, sind, cosd0, cosd = sin(dec1), sin(dec2), cos(dec1), cos(dec2)
    sina, cosa = sin(ra) * cosd, cos(ra) * cosd
    x = cosa * sind0 - sind * cosd0
    y = sina
    z = cosa * cosd0 + sind * sind0
    print(x, y, z)
    PA = numpy.arctan2(y, -x)
    R = numpy.arccos(z)
    return R, PA


def _deg_to_dms(x, prec=0.01):
    """Converts x (in degrees) into d,m,s tuple, where d and m are ints.
    prec gives the precision, in arcseconds."""
    mins, secs = divmod(round(x * 3600 / prec) * prec, 60)
    mins = int(mins)
    degs, mins = divmod(mins, 60)
    return degs, mins, secs


def ra_hms(rad, scale=12, prec=0.01):
    """Returns RA as tuple of (h,m,s)"""
    # convert negative values
    while rad < 0:
        rad += 2 * math.pi
    # convert to hours
    rad *= scale / math.pi
    return _deg_to_dms(rad, prec)


def dec_dms(rad, prec=0.01):
    return dec_sdms(rad, prec)[1:]


def dec_sdms(rad, prec=0.01):
    """Returns Dec as tuple of (sign,d,m,s). Sign is "+" or "-"."""
    sign = "-" if rad < 0 else "+"
    d, m, s = _deg_to_dms(abs(rad) / DEG, prec)
    return (sign, d, m, s)


def ra_hms_string(rad):
    return "%dh%02dm%05.2fs" % ra_hms(rad)


def dec_sdms_string(rad):
    return "%s%dd%02dm%05.2fs" % dec_sdms(rad)


def radec_string(ra, dec):
    return "%s %s" % (ra_hms_string(ra), dec_sdms_string(dec))


class _Projector(object):
    """This is an abstract base class for all projection classes below. A projection class can be used to create projector objects for
    conversion between world (ra,dec) and projected (l,m) coordinates.

    * A projector is instantiated as proj = Proj(ra0,dec0)      # ra0,dec0 is projection centre
    * converts ra,dec->l,m as
          l,m = proj.lm(ra,dec)
    * converts l,m->ra,dec as
          ra,dec = proj.radec(l,m)
    * converts angular offsets (from 0,0 point) into l,m:
          l,m = proj.offset(dra,ddec)

    Alternativelty, there are class methods which do not require one to instantiate a projector object:

    * Proj.radec_lm(ra,dec,ra0,dec0)
    * Proj.lm_radec(l,m,ra0,dec0)
    * Proj.offset_lm(dra,ddec,ra0,dec0)
    """

    def __init__(self, ra0, dec0, has_projection=False):
        self.ra0, self.dec0, self.sin_dec0, self.cos_dec0 = ra0, dec0, sin(dec0), cos(dec0)
        self._has_projection = has_projection

    def has_projection(self):
        return bool(self._has_projection)

    def __eq__(self, other):
        """By default, two projections are the same if their classes match, and their ra0/dec0 match."""
        return type(self) is type(other) and self.ra0 == other.ra0 and self.dec0 == other.dec0

    def __ne__(self, other):
        return not self == other

    @classmethod
    def radec_lm(cls, ra, dec, ra0, dec0):
        return cls(ra0, dec0).lm(ra, dec)

    @classmethod
    def lm_radec(cls, l, m, ra0, dec0):
        return cls(ra0, dec0).radec(l, m)

    @classmethod
    def offset_lm(cls, dra, ddec, ra0, dec0):
        return cls(ra0, dec0).offset(dra, ddec)

    def lm(self, ra, dec):
        raise TypeError("lm() not yet implemented in projection %s" % type(self).__name__)

    def offset(self, dra, ddec):
        raise TypeError("offset() not yet implemented in projection %s" % type(self).__name__)

    def radec(self, l, m):
        raise TypeError("radec() not yet implemented in projection %s" % type(self).__name__)


def get_wcs_info(hdr):
    naxis = hdr['NAXIS']
    ra_axis = dec_axis = None
    refpix = [hdr["CRPIX%d" % (iaxis+1)]-1 for iaxis in range(naxis)]
    for iaxis in range(naxis):
        name = hdr.get("CTYPE%d" % (iaxis+1), '').upper()
        if name.startswith("RA"):
            ra_axis = iaxis
        elif name.startswith("DEC"):
            dec_axis = iaxis
    wcs = WCS(hdr)  
    refsky = wcs.wcs_pix2world([refpix], 0)[0,:]    
    return wcs, refpix, refsky, ra_axis, dec_axis


class Projection(object):
    """Projection is a container for the different projection classes.
    Each Projection class can be used to create a projection object: proj = Proj(ra0,dec0), with lm(ra,dec) and radec(l,m) methods.
    """

    class FITSWCSpix(_Projector):
        """FITS WCS projection, as determined by a FITS header. lm is in pixels (0-based)."""

        def __init__(self, header):
            """Constructor. Create from filename (treated as FITS file), or a FITS header object"""
            # attach to FITS file or header
            if isinstance(header, str):
                header = pyfits.open(header)[0].header

            try:
                self.wcs, self.refpix, self.refsky, self.ra_axis, self.dec_axis = get_wcs_info(header)
                if self.ra_axis is None or self.dec_axis is None:
                    raise RuntimeError("Missing RA or DEC axis")
                ra0, dec0 = self.refsky[self.ra_axis], self.refsky[self.dec_axis]
                self.xpix0, self.ypix0 = self.refpix[self.ra_axis], self.refpix[self.dec_axis]
                refpix1 = self.refpix.copy()
                refpix1[self.ra_axis] += 1
                refpix1[self.dec_axis] += 1
                delta = self.wcs.wcs_pix2world([refpix1], 0)[0] - self.refsky
                self.xscale = delta[self.ra_axis] * DEG
                self.yscale = delta[self.dec_axis] * DEG
                has_projection = True
            except Exception as exc:
                traceback.print_exc()
                print("No WCS in FITS file, falling back to pixel coordinates.")
                ra0 = dec0 = self.xpix0 = self.ypix0 = 0
                self.xscale = self.yscale = DEG / 3600
                has_projection = False
            _Projector.__init__(self, ra0 * DEG, dec0 * DEG, has_projection=has_projection)

        def lm(self, ra, dec):
            if not self.has_projection():
                return numpy.sin(ra) / self.xscale, numpy.sin(dec) / self.yscale
            if numpy.isscalar(ra) and numpy.isscalar(dec):
                if ra - self.ra0 > math.pi:
                    ra -= 2 * math.pi
                if ra - self.ra0 < -math.pi:
                    ra += 2 * math.pi
                skyvec = self.refsky.copy()
                skyvec[self.ra_axis] = ra / DEG
                skyvec[self.dec_axis] = dec / DEG
                pixvec = self.wcs.wcs_world2pix([skyvec], 0)[0] 
                return pixvec[self.ra_axis], pixvec[self.dec_axis]
            else:
                if numpy.isscalar(ra):
                    ra = numpy.array(ra)
                if numpy.isscalar(dec):
                    dec = numpy.array(dec)
                n = max(len(ra), len(dec))
                skymat = numpy.array([self.refsky for _ in range(n)])
                skymat[:, self.ra_axis] = ra / DEG
                skymat[:, self.dec_axis] = dec / DEG
                ra = skymat[:, self.ra_axis]
                ra[ra - self.ra0 > 180] -= 360
                ra[ra - self.ra0 < -180] += 360
                ## when fed in arrays of ra/dec, wcs.wcs2pix will return a nested list of
                ## [[l1,m1],[l2,m2],,...]. Convert this to an array and extract columns.
                lm = self.wcs.wcs_world2pix(skymat, 0)
                return lm[:, self.ra_axis], lm[:, self.dec_axis]

        def radec(self, l, m):
            if not self.has_projection():
                return numpy.arcsin(l * self.xscale), numpy.arcsin(m * self.yscale)
            if numpy.isscalar(l) and numpy.isscalar(m):
                pixvec = self.refpix.copy()
                pixvec[self.ra_axis] = l
                pixvec[self.dec_axis] = m
                skyvec = self.wcs.wcs_pix2world([pixvec], 0)[0]
                ra, dec = skyvec[self.ra_axis], skyvec[self.dec_axis]
            else:
                ## this is slow as molasses because of the way astLib.WCS implements the loop. ~120 seconds for 4M pixels
                ## when fed in arrays of ra/dec, wcs.wcs2pix will return a nested list of
                ## [[l1,m1],[l2,m2],,...]. Convert this to an array and extract columns.
                #        radec = numpy.array(self.wcs.pix2wcs(l,m))
                #        ra = radec[...,0]
                #        dec = radec[...,1]
                ### try a faster implementation -- oh well, only a bit faster, ~95 seconds for the same
                ### can also replace list comprehension with map(), but that doesn't improve things.
                ### Note also that the final array constructor takes ~10 secs!
                radec = numpy.array(
                    [PyWCSTools.wcs.pix2wcs(self.wcs.WCSStructure, x, y) for x, y in zip(l + 1, m + 1)])
                ra = radec[..., 0]
                dec = radec[..., 1]
            return ra * DEG, dec * DEG

        def offset(self, dra, ddec):
            return self.xpix0 + dra / self.xscale, self.ypix0 + ddec / self.xscale

        def __eq__(self, other):
            """By default, two projections are the same if their classes match, and their ra0/dec0 match."""
            return type(self) is type(other) and (
            self.ra0, self.dec0, self.xpix0, self.ypix0, self.xscale, self.yscale) == (
                   other.ra0, other.dec0, other.xpix0, other.ypix0, other.xscale, other.yscale)


    ## OMS 9/2/2021: Retiring FITSWCS, as it was only being used as a base for SinWCS() before, so it's cleaner to do SinWCS directly.
    ## There is one place that Tigger uses FITSWCS, but only to get the header info, not for coordinate conversions.

    class SinWCS(FITSWCSpix):
        """
        A sin WCS projection centred on the given ra0/dec0 coordinates,
        with 0,0 being the reference pixel, 
        """

        def __init__(self, ra0, dec0):
            hdu = pyfits.PrimaryHDU()
            hdu.header.set('NAXIS', 2)
            hdu.header.set('NAXIS1', 3)
            hdu.header.set('NAXIS2', 3)
            hdu.header.set('CTYPE1', 'RA---SIN')
            hdu.header.set('CDELT1', -1. / 60)
            hdu.header.set('CRPIX1', 2)
            hdu.header.set('CRVAL1', ra0 / DEG)
            hdu.header.set('CUNIT1', 'deg     ')
            hdu.header.set('CTYPE2', 'DEC--SIN')
            hdu.header.set('CDELT2', 1. / 60)
            hdu.header.set('CRPIX2', 2)
            hdu.header.set('CRVAL2', dec0 / DEG)
            hdu.header.set('CUNIT2', 'deg     ')
            Projection.FITSWCSpix.__init__(self, hdu.header)
            self._l0 = self.refpix[self.ra_axis]
            self._m0 = self.refpix[self.dec_axis]

        def lm(self, ra, dec):
            l, m = super().lm(ra, dec)
            return sin((l - self._l0) * self.xscale), sin((m - self._m0)*self.yscale)

        def radec(self, l, m):
            return super().radec(arcsin(l / self.xscale + self._l0), arcsin(m / self.yscale + self._m0))

        def offset(self, dra, ddec):
            return sin(dra), sin(ddec)

