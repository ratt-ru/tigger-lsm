#!/usr/bin/env python3
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

import os
import sys

from astropy.io import fits as pyfits
from past.builtins import cmp

def main():

    import Tigger.Models.Formats
    from Tigger.Models.Formats import ASCII

    AUTO = "auto"
    full_formats = Tigger.Models.Formats.listFormatsFull()
    input_formats = [name for name, (load, save, doc, extensions) in full_formats if load] + [AUTO]

    # setup some standard command-line option parsing
    #
    from optparse import OptionParser

    parser = OptionParser(usage="""%prog: [options] input_image sky_model [output_image]""",
                          description="""Restores sources from sky model into an input image, writes result to output image. If
an output image is not specified, makes a name for it automatically.""")
    parser.add_option("-t", "--type", choices=input_formats,
                      help="Input model type (%s). Default is %%default." % (", ".join(input_formats)))
    parser.add_option("--format", type="string",
                      help="""Input format, for ASCII or BBS tables. For ASCII tables, default is "%s". For BBS tables, the default format is specified in the file header.""" % ASCII.DefaultDMSFormatString)
    parser.add_option("-n", "--num-sources", dest="nsrc", type="int", action="store",
                      help="Only restore the NSRC brightest sources")
    parser.add_option("-s", "--scale", dest="fluxscale", metavar="FLUXSCALE[,N]", action="store",
                      help="rescale model fluxes by given factor. If N is given, rescale N brightest only.")
    parser.add_option("-b", "--restoring-beam", type="string", metavar="BMAJ[,BMIN,PA]",
                      help="specify restoring beam size, overriding BMAJ/BMIN/BPA keywords in input image. " +
                           "Use a single value (arcsec) for circular beam, or else " +
                           "supply major/minor size and position angle (deg).")
    parser.add_option("-p", "--psf-file", dest="psf", action="store",
                      help="determine restoring beam size by fitting PSF file, overriding BMAJ/BMIN/BPA keywords in input image.")
    parser.add_option("--clear", action="store_true",
                      help="clear contents of FITS file before adding in sources")
    parser.add_option("--pb", action="store_true",
                      help="apply model primary beam function during restoration, if it's defined, and source is not tagged 'nobeam'")
    parser.add_option("--beamgain", action="store_true",
                      help="apply beamgain atribute during restoration, if it's defined, and source is not tagged 'nobeam'")
    parser.add_option("--ignore-nobeam", action="store_true",
                      help="apply PB or beamgain even if source is tagged 'nobeam'")
    parser.add_option("-F", "--freq", type="float", metavar="MHz", default=0,
                      help="use this frequency (for spectral indices and primary beams)")
    parser.add_option("-f", dest="force", action="store_true",
                      help="overwrite output image even if it already exists")
    parser.add_option("-v", "--verbose", dest="verbose", type="int", action="store",
                      help="set verbosity level (0 is silent, higher numbers mean more messages)")
    parser.add_option("-T", "--timestamps", action="store_true",
                      help="enable timestamps in debug messages (useful for timing)")
    parser.set_defaults(n=0, fluxscale='1')

    (options, rem_args) = parser.parse_args()

    # get filenames
    if len(rem_args) == 2:
        input_image, skymodel = rem_args
        name, ext = os.path.splitext(input_image)
        output_image = name + ".restored" + ext
    elif len(rem_args) == 3:
        input_image, skymodel, output_image = rem_args
    else:
        parser.error("Insufficient number of arguments. Use -h for help.")

    # check for overwritten output
    if os.path.exists(output_image) and not options.force:
        parser.error("File %s already exists, use the -f option to overwrite." % output_image)

    # find Tigger
    try:
        import Tigger
    except ImportError:
        sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
        try:
            import Tigger
        except:
            print("Unable to import the Tigger package. Please check your installation and PYTHONPATH.")
            sys.exit(1)

    #Tigger.nuke_matplotlib();  # don't let the door hit you in the ass, sucka

    from Tigger.Tools import Imaging
    from Tigger.Tools.Imaging import FWHM, DEG, ARCSEC

    Imaging._verbosity.set_verbose(options.verbose)
    Imaging._verbosity.enable_timestamps(options.timestamps)

    # read model and sort by apparent brightness
    # figure out input type
    try:
        input_type, import_func, dum, input_doc = Tigger.Models.Formats.resolveFormat(skymodel,
                                                                                      options.type if options.type != AUTO else None)
    except:
        print("Unable to determine model type for %s, please specify one explicitly with the -t/--type option." % skymodel)
        sys.exit(1)

    print("Reading %s (%s)" % (skymodel, input_doc))
    model = import_func(skymodel, format=options.format)

    Imaging.dprintf(1, "Read %d sources from %s\n", len(model.sources), skymodel)

    sources = sorted(model.sources, key=lambda a: a.brightness()) #, lambda a, b: cmp(b.brightness(), a.brightness()))

    # apply counts and flux scales
    if options.nsrc:
        sources = sources[:options.nsrc]
        Imaging.dprintf(1, "Using %d brightest sources\n", len(sources))

    if options.fluxscale != '1':
        if "," in options.fluxscale:
            scale, n = options.fluxscale.split(",", 1)
            scale = float(scale)
            n = int(n)
            Imaging.dprintf(1, "Flux of %d brightest sources will be scaled by %f\n", n, scale)
        else:
            scale = float(options.fluxscale)
            n = len(sources)
            Imaging.dprintf(1, "Flux of all model sources will be scaled by %f\n", n, scale)
        for src in sources[:n]:
            src.flux.rescale(0.01)

    # open input image
    input_hdu = pyfits.open(input_image)[0]

    # get restoring beam size
    if options.restoring_beam:
        ff = options.restoring_beam.split(",")
        try:
            if len(ff) == 1:
                gx = gy = float(ff[0])
                grot = 0
                print("User-specified restoring beam of %.2f\"" % gx)
            else:
                gx, gy, grot = list(map(float, ff))
                print("User-specified restoring beam of %.2f\" by %.2f\" at PA %.2f deg" % (gx, gy, grot))
        except:
            print("Invalid -b/--restoring-beam setting.")
            sys.exit(1)
        gx /= FWHM * ARCSEC
        gy /= FWHM * ARCSEC
        grot /= DEG
    elif options.psf:
        # fit the PSF
        gx, gy, grot = Imaging.fitPsf(options.psf)
        print("Fitted restoring beam to PSF file %s: %.2f\" by %.2f\" at PA %.2f deg" % (
        options.psf, gx * FWHM * ARCSEC, gy * FWHM * ARCSEC, grot * DEG))
    else:
        # else look in input header
        gx, gy, grot = [input_hdu.header.get(x, None) for x in ('BMAJ', 'BMIN', 'BPA')]
        if any([x is None for x in (gx, gy, grot)]):
            print("Unable to determine restoring beam size, no BMAJ/BMIN/BPA keywords in input image.")
            print("Try using the -b/-p options to specify an explicit restoring beam.")
            sys.exit(1)
        print("Restoring beam (as per input header) is %.2f\" by %.2f\" at PA %.2f deg" % (gx * 3600, gy * 3600, grot))
        gx /= DEG * FWHM
        gy /= DEG * FWHM
        grot /= DEG

    pbexp = None
    freq = options.freq * 1e+6 or model.refFreq() or 1400 * 1e+6

    if options.pb and model.primaryBeam():
        try:
            pbexp = eval('lambda r,fq:' + model.primaryBeam())
            dum = pbexp(0, 1e+9)  # evaluate at r=0 and 1 GHz as a test
            if not isinstance(dum, float):
                raise TypeError("Primary beam expression does not evaluate to a float")
        except Exception as exc:
            print("Bad primary beam expression '%s': %s" % (options.pb, str(exc)))
            sys.exit(1)
        if not freq:
            print("Model must contain a reference requency, or else specify one with --freq.")
            sys.exit(1)

    # read, restore, write
    print("Restoring model into input image %s" % input_image)
    if options.clear:
        input_hdu.data[...] = 0
    Imaging.restoreSources(input_hdu, sources, gx, gy, grot, primary_beam=pbexp, freq=freq,
                           apply_beamgain=options.beamgain, ignore_nobeam=options.ignore_nobeam)

    print("Writing output image %s" % output_image)
    if os.path.exists(output_image):
        os.remove(output_image)
    input_hdu.writeto(output_image)
