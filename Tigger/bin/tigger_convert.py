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

import glob
import math
import re
import sys
import traceback

import numpy as np
import os.path
import Cattery.Siamese.OMS.Utils as Utils
from Tigger.Coordinates import get_wcs_info

# add brick
from astropy.io import fits
from astropy.wcs import WCS

DEG = math.pi / 180

NATIVE = "Tigger"


def Jones2Mueller_circular(J):
    S = np.matrix([[1, 0, 0, 1], [0, 1, 1j, 0], [0, 1, -1j, 0], [1, 0, 0, -1]])
    # Compute the Mueller matrix
    MM = (S.I) * np.kron(J, J.H) * S
    return np.real(MM)


def Jones2Mueller_linear(J):
    S = np.matrix([[1, 1, 0, 0], [0, 0, 1, 1j], [0, 0, 1, -1j], [1, -1, 0, 0]])
    # Compute the Mueller matrix
    MM = (S.I) * np.kron(J, J.H) * S
    return np.real(MM)


## Griffin's old version, for linear. Possibly the order is wrong
#  A=np.matrix([[1,0,0,1],[1,0,0,-1],[0,1,1,0],[0,1j,-1j,0]])
#  M=A*np.kron(J,J.conj())*np.linalg.inv(A)
#  return np.real(M)

def arc2lm(l0, m0, arclen=2. * np.pi, nsteps=360):
    """Return cartesian positions that sample an arc of a circle (similar to np.linspace)
    l0,m0: initial cartesian position to determine radius and starting point
    arclen: angle, in radians, to sample, value should be between 0 and 2pi
    nsteps: number of samples"""
    r = np.sqrt(float(l0) ** 2. + float(m0) ** 2.)
    angle = np.arctan2(m0, l0)
    da = np.linspace(0., arclen, num=nsteps)
    l = r * np.cos(angle + da)
    m = r * np.sin(angle + da)
    return l, m


def rotatelm(l0, m0, rotangle):
    """Rotate (l0,m0) to a new (l,m) based on angle"""
    r = np.sqrt(float(l0) ** 2. + float(m0) ** 2.)
    angle = np.arctan2(m0, l0)
    l = r * np.cos(angle + rotangle)
    m = r * np.sin(angle + rotangle)
    return l, m


def main():
    import Kittens.utils

    _verbosity = Kittens.utils.verbosity(name="convert-model")
    dprint = _verbosity.dprint
    dprintf = _verbosity.dprintf

    # find Tigger
    try:
        import Tigger
    except ImportError:
        dirname = os.path.dirname(os.path.realpath(__file__))
        # go up the directory tree looking for directory "Tigger"
        while len(dirname) > 1:
            if os.path.basename(dirname) == "Tigger":
                break
            dirname = os.path.dirname(dirname)
        else:
            print("Unable to locate the Tigger directory, it is not a parent of %s. Please check your installation"
                    "and/or PYTHONPATH." % os.path.realpath( __file__))
            sys.exit(1)
        sys.path.append(os.path.dirname(dirname))
        try:
            import Tigger
        except:
            print("Unable to import the Tigger package from %s. Please check your installation and PYTHONPATH." %
                    dirname)
            sys.exit(1)

    # some things can implicitly invoke matplotlib, which can cry when no X11 is around
    # so to make sure thingfs work in pipelines, we explicitly disable this here, unless we're asked for plots
    if not "--enable-plots" in sys.argv:
        Tigger.nuke_matplotlib();  # don't let the door hit you in the ass, sucka

    from Tigger import Coordinates
    import Tigger.Models.Formats
    import Tigger.Models.ModelClasses

    AUTO = "auto"
    full_formats = Tigger.Models.Formats.listFormatsFull()
    input_formats = [name for name, (load, save, doc, extensions) in full_formats if load] + [AUTO]
    output_formats = [name for name, (load, save, doc, extensions) in full_formats if save] + [AUTO]

    from Tigger.Models.Formats import ASCII

    # setup some standard command-line option parsing
    #
    from optparse import OptionParser, OptionGroup

    parser = OptionParser(usage="""%prog: sky_model [output_model]""",
                          description="""Converts sky models into Tigger format and/or applies various processing options.
Input 'sky_model' may be any model format importable by Tigger, recognized by extension, or explicitly specified via an option switch.
'output_model' is always a native Tigger model. If an output model is not specfied, the conversion is done in-place if the input model
is a Tigger model (-f switch must be specified to allow overwriting), or else a new filename is generated.""")

    group = OptionGroup(parser, "Input/output and conversion options")
    parser.add_option_group(group)
    group.add_option("-f", "--force", action="store_true",
                     help="Forces overwrite of output model.")
    group.add_option("-t", "--type", choices=input_formats,
                     help="Input model type (%s). Default is %%default." % (", ".join(input_formats)))
    group.add_option("-o", "--output-type", choices=output_formats, metavar="TYPE",
                     help="Output model type (%s). Default is %%default." % (", ".join(output_formats)))
    group.add_option("-a", "--append", metavar="FILENAME", action="append",
                     help="Append another model to input model. May be given multiple times.")
    group.add_option("--append-type", choices=input_formats, metavar="TYPE",
                     help="Appended model type (%s). Default is %%default." % (", ".join(input_formats)))
    group.add_option("--format", type="string",
                     help="""Input format, for ASCII or BBS tables. For ASCII tables, default is "%s". For BBS tables, the default format is specified in the file header.""" % ASCII.DefaultDMSFormatString)
    group.add_option("--append-format", type="string", default="",
                     help="""Format of appended file, for ASCII or BBS tables. Default is to use --format.""")
    group.add_option("--output-format", type="string", metavar="FORMAT",
                     help="""Output format, for ASCII or BBS tables. If the model was originally imported from an ASCII or BBS table, the default output format will be the same as the original format.""")
    group.add_option("--help-format", action="store_true",
                     help="Prints help on format strings.")
    group.add_option("--min-extent", type="float", metavar="ARCSEC",
                     help="Minimal source extent, when importing NEWSTAR or ASCII files. Sources with a smaller extent will be treated as point sources. Default is %default.")

    group = OptionGroup(parser, "Options to select a subset of the input")
    parser.add_option_group(group)
    group.add_option("-T", "--tags", type="string", action="append", metavar="TAG",
                     help="Extract sources with the specified tags.")
    group.add_option("--select", type="string", metavar='TAG<>VALUE', action="append",
                     help="Selects a subset of sources by comparing the named TAG to a float VALUE. '<>' " +
                          "represents the comparison operator, and can be one of == (or =),!=,<=,<,>,>=. Alternatively, " +
                          "you may use the FORTRAN-style operators .eq.,.ne.,.le.,.lt.,.gt.,.ge. Multiple " +
                          "select options may be given, in which case the effect is a logical-AND. Note that VALUE may be "
                          "followed by one of the characters d, m or s, in which case it will be converted from degrees, "
                          "minutes or seconds into radians. This is useful for selections such as \"r<5d\".")
    group.add_option("--remove-nans", action="store_true",
                     help="Removes the named source(s) from the model. NAME may contain * and ? wildcards.")

    group = OptionGroup(parser, "Options to manipulate fluxes etc.")
    parser.add_option_group(group)
    group.add_option("--app-to-int", action="store_true",
                     help="Treat fluxes as apparent, and rescale them into intrinsic using the " +
                          "supplied primary beam model (see --primary-beam option).")
    group.add_option("--int-to-app", action="store_true",
                     help="Treat fluxes as intrinsic, and rescale them into apparent using the " +
                          "supplied primary beam model (see --primary-beam option).")
    group.add_option("--newstar-app-to-int", action="store_true",
                     help="Convert NEWSTAR apparent fluxes in input model to intrinsic. Only works for NEWSTAR or NEWSTAR-derived input models.")
    group.add_option("--newstar-int-to-app", action="store_true",
                     help="Convert NEWSTAR intrinsic fluxes in input model to apparent. Only works for NEWSTAR or NEWSTAR-derived input models.")
    group.add_option("--center", type="string", metavar='COORDINATES',
                     help="Override coordinates of the nominal field center specified in the input model. Use the form " +
                          "\"Xdeg,Ydeg\" or \"Xdeg,Yrad\" to specify RA,Dec in degrees or radians, or else a " +
                          "a pyrap.measures direction string of the form " + \
                          "REF,C1,C2, for example \"j2000,1h5m0.2s,+30d14m15s\". See the pyrap.measures documentation for more details.")
    group.add_option("--refresh-r", action="store_true",
                     help="Recompute the 'r' (radial distance from center) attribute of each source based on the current field center.")
    group.add_option("--ref-freq", type="float", metavar="MHz",
                     help="Set or change the reference frequency of the model.")

    group = OptionGroup(parser, "Primary beam-related options")
    parser.add_option_group(group)
    group.add_option("--primary-beam", type="string", metavar="EXPR",
                     help="""Apply a primary beam expression to estimate apparent fluxes. Any valid Python expression using the variables 'r' and 'fq' is accepted. Use "refresh" to re-estimate fluxes using the current expression.
                    Example (for the WSRT-like 25m dish PB): "cos(min(65*fq*1e-9*r,1.0881))**6".
                    OR: give a set of FITS primary beam patterns of the form e.g. FILENAME_$(xy)_$(reim).fits, these are the same FITS files used in MeqTrees pybeams_fits.""")
    group.add_option("--linear-pol", action="store_true",
                     help="Use XY basis correlations for beam filenames and Mueller matrices. Default is RL.")
    group.add_option("--fits-l-axis", type="string", default="-X",
                     help="CTYPE for L axis in the FITS PB file. Note that our internal L points East (increasing RA), if the "
                          "FITS beam axis points the opposite way, prefix the CTYPE with a '-'' character.")
    group.add_option("--fits-m-axis", type="string", default="Y",
                     help="CTYPE for M axis in the FITS PB file. Note that our internal M points North (increasing Dec), if the "
                          "FITS beam axis points the opposite way, prefix the CTYPE with a '-'' character.")
    group.add_option("--beam-freq", type="float", metavar="MHz",
                     help="use given frequency for primary beam model, rather than the model reference frequency")
    group.add_option("--beam-clip", type="float", metavar="GAIN", default=0.001,
                     help="when using a FITS beam, clip (power) beam gains at this level to keep intrinsic source fluxes from blowing up. Sources below this beamgain will be tagged 'nobeam'. Default: %default")
    group.add_option("--beam-spi", type="float", metavar="MHz",
                     help="perform a spectral index fit to each source based on a frequency dependent FITS beam, requires --primary-beam option to be used with a FITS file. " +
                          "Apply this spectral index to LSM sources. " +
                          "Must supply a band width (centred on --beam-freq) over which the beam spi is estimated")
    group.add_option("--force-beam-spi-wo-spectrum", action="store_true",
                     help="apply beam-derived spectral indices even to sources without an intrinsic spectrum. Default " +
                          "is to only apply to sources that already have a spectrum."
                     )
    group.add_option("--beam-nopol", action="store_true",
                     help="apply intensity beam model only, ignoring polarization. Default is to use polarization."
                     )
    group.add_option("--beam-diag", action="store_true",
                     help="use diagonal Jones terms only for beam model. Default is to use all four terms if available."
                     )
    group.add_option("--pa", type="float", default=None,
                     help="Rotate the primary beam pattern through a parallactic angle (in degrees).")
    group.add_option("--pa-range", type="str", default=None, metavar="FROM,TO",
                     help="Rotate the primary beam pattern through a range of parallactic angles (in degrees) and use the average value over PA.")
    group.add_option("--pa-from-ms", type="str", default=None, metavar="MS1[:FIELD1],MS2:[FIELD2],...",
                     help="Rotate the primary beam pattern through a range of parallactic angles as given by the MS and field ID (default 0), " +
                          "and take the average over time. This is more accurate than --pa-range.")
    group.add_option("--beam-average-jones", action="store_true",
                     help="Correct approach to rotational averaging is to convert Jones(PA) to Mueller(PA), then average " +
                          "over PA. Tigger versions<=1.3.3 took the incorrect approach of averaging Jones over PA, then converting " +
                          "to Mueller. Use this option to mimic the old approach.")

    group = OptionGroup(parser, "Options to cluster and rename sources")
    parser.add_option_group(group)
    group.add_option("--cluster-dist", type="float", metavar="ARCSEC",
                     help="Distance parameter for source clustering, 0 to disable. Default is %default.")
    group.add_option("--rename", action="store_true",
                     help="Rename sources according to the COPART (cluster ordering, P.A., radial distance, type) scheme")
    group.add_option("--radial-step", type="float", metavar="ARCMIN",
                     help="Size of one step in radial distance for the COPART scheme. Default is %default'.")
    group.add_option("--merge-clusters", type="string", metavar="TAG(S)",
                     help="Merge source clusters bearing the specified tags, replacing them with a " + "single point source. Multiple tags may be given separated by commas. " +
                          "Use 'ALL' to merge all clusters.")
    group.add_option("--prefix", type="string",
                     help="Prefix all source names with the given string")

    group = OptionGroup(parser, "Other model manipulation options")
    parser.add_option_group(group)
    group.add_option("--remove-source", type="string", action="append",
                     metavar="NAME",
                     help="Removes the named source(s) from the model. NAME may contain * and ? wildcards.")
    group.add_option("--add-brick", type="string", action="append",
                     metavar="NAME:FILE[:PAD_FACTOR:[TAGS:...]]",
                     help="Adds a uv-brick to the model. NAME is a source name, FILE is a " +
                          "FITS file, PAD_FACTOR is set to 1 if not specified. TAGS is a list of boolean flags.")
    group.add_option("--recenter", type="string", metavar='COORDINATES',
                     help="Shift the sky model from the nominal center to a different field center. COORDINATES specified as per the --center option.")

    group = OptionGroup(parser, "Debugging and verbosity options")
    parser.add_option_group(group)
    group.add_option("-v", "--verbose", action="count",
                     help="increases verbosity.")
    group.add_option("-d", "--debug", dest="debug", type="string", action="append", metavar="Context=Level",
                     help="(for debugging Python code) sets verbosity level of the named Python context. May be used multiple times.")
    group.add_option("--enable-plots", action="store_true",
                     help="enables various diagnostic plots")

    parser.set_defaults(cluster_dist=60, min_extent=0, format=None, type='auto', output_type='auto', radial_step=10,
                        ref_freq=-1)

    (options, rem_args) = parser.parse_args()
    min_extent = (options.min_extent / 3600) * DEG

    if options.help_format:
        print(ASCII.FormatHelp)
        sys.exit(0)

    # get filenames
    if len(rem_args) == 1:
        skymodel = rem_args[0]
        output = None
    elif len(rem_args) == 2:
        skymodel, output = rem_args
    else:
        parser.error("Incorrect number of arguments. Use -h for help.")

    if options.app_to_int and options.int_to_app:
        parser.error("Can't use --app-to-int and --int-to-app together.")
    if options.newstar_app_to_int and options.newstar_int_to_app:
        parser.error("Can't use --newstar-app-to-int and --newstar-int-to-app together.")

    global measures_dmdq
    measures_dmdq = None


    def pyrap_dmdq():
        """Helper function: imports pyrap.measures, and returns dm,dq objects"""
        global measures_dmdq
        if measures_dmdq is None:
            try:
                import pyrap.measures
                import pyrap.quanta
            except:
                traceback.print_exc()
                print("Failed to import pyrap.measures, which is required by one of the options you specified.")
                print("You probably need to install the 'pyrap' package for this to work.")
                sys.exit(1)
            measures_dmdq = pyrap.measures.measures(), pyrap.quanta
        return measures_dmdq


    def convert_coordinates(coords):
        """Converts a measures coordinate string into a ra,dec pair (radians at J2000)"""
        match = re.match("^([\d.]+)(rad|deg|),([-]?[\d.]+)(rad|deg|)$", coords)
        if match:
            ra = float(match.group(1))
            dec = float(match.group(3))
            return ra * (DEG if match.group(2) == "deg" else 1), dec * (DEG if match.group(4) == "deg" else 1)
        dm, dq = pyrap_dmdq()
        try:
            coord_dir = dm.direction(*(coords.split(',')))
            coord_dir = dm.measure(coord_dir, 'j2000')
            qq = dm.get_value(coord_dir)
            return [q.get_value('rad') for q in qq]
        except:
            print("Error parsing or converting coordinate string '%s', see traceback:" % coords)
            traceback.print_exc()
            sys.exit(1)


    # figure out center and recenter option
    if options.recenter:
        recenter_radec = convert_coordinates(options.recenter)
    if options.center:
        center_radec = convert_coordinates(options.center)
        options.refresh_r = True
    else:
        center_radec = None

    # check the 'select' option
    select_predicates = {
        '=': lambda x, y: x == y,
        '==': lambda x, y: x == y,
        '!=': lambda x, y: x != y,
        '>=': lambda x, y: x >= y,
        '<=': lambda x, y: x <= y,
        '>': lambda x, y: x > y,
        '<': lambda x, y: x < y,
        '.eq.': lambda x, y: x == y,
        '.ne.': lambda x, y: x != y,
        '.ge.': lambda x, y: x >= y,
        '.le.': lambda x, y: x <= y,
        '.gt.': lambda x, y: x > y,
        '.lt.': lambda x, y: x < y
    }
    select_units = dict(d=DEG, m=DEG / 60, s=DEG / 3600)

    selections = []
    for selstr in (options.select or []):
        match = re.match("^(?i)([^=<>!.]+)(%s)([^dms]+)([dms])?" % "|".join(
            [key.replace('.', '\.') for key in list(select_predicates.keys())]), selstr)
        if not match:
            parser.error("Malformed --select string '%s'." % selstr)
        try:
            value = float(match.group(3))
        except:
            parser.error("Malformed --select string '%s': right-hand side is not a number." % selstr)
        scale = select_units.get(match.group(4), 1.)
        selections.append((selstr, match.group(1), select_predicates[match.group(2).lower()], value * scale))

    # figure out input type
    try:
        input_type, import_func, dum, input_doc = Tigger.Models.Formats.resolveFormat(skymodel,
                                                                                      options.type if options.type != AUTO else None, io="input")
    except:
        print("Unable to determine model type for %s, please specify one explicitly with the -t/--type option." % skymodel)
        sys.exit(1)

    # figure out output type, if explicitly specified
    output_type = None
    if output is None and options.output_type == AUTO:
        options.output_type = "Tigger"

    if options.output_type != AUTO:
        output_type, dum, export_func, output_doc = Tigger.Models.Formats.getFormat(options.output_type)
        output_extensions = Tigger.Models.Formats.getFormatExtensions(options.output_type)
        if not export_func:# or not extensions: #not defined @oms
            print("Output model type '%s' is not supported." % options.output_type)
            sys.exit(1)

    # figure out output name, if not specified
    if output is None:
        if not output_type:
            print("An output filename and/or an explicit output model type (-o/--output-type) must be specfified.")
            sys.exit(1)
        # get base input name
        # if input extension is "lsm.html", then split off two extensions, not just one
        basename = os.path.splitext(skymodel)[0]
        if skymodel.endswith(".lsm.html"):
            basename = os.path.splitext(basename)[0]
        output = basename + output_extensions[0]
    # else output name is specified, use this to determine format unless it is explicitly set
    elif not output_type:
        try:
            output_type, dum, export_func, output_doc = Tigger.Models.Formats.resolveFormat(output, None, io="output")
        except:
            export_func = None
        if not export_func:
            print("Unable to determine model type for %s, please specify one explicitly with the -o/--output-type option." % output)
            sys.exit(1)

    # check if we need to overwrite
    if os.path.exists(output) and not options.force:
        print("Output file %s already exists. Use the -f switch to overwrite." % output)
        sys.exit(1)

    print("Reading %s (%s)" % (skymodel, input_doc))

    # load the model
    try:
        model = import_func(skymodel, min_extent=min_extent, format=options.format, center=center_radec,
                            verbose=options.verbose)
    except Exception as exc:
        if options.verbose:
            traceback.print_exc()
        print("Error loading model:", str(exc))
        sys.exit(1)
    sources = model.sources
    if not sources:
        print("Input model %s contains no sources" % skymodel)
    else:
        print("Model contains %d sources" % len(sources))

    # append, if specified
    if options.append:
        for modelnum, filename in enumerate(options.append):
            # figure out input type
            try:
                append_type, append_func, dum, append_doc = Tigger.Models.Formats.resolveFormat(filename,
                                                                                                options.append_type if options.append_type != AUTO else None, io="input")
            except:
                print("Unable to determine model type for %s, please specify one explicitly with the --append-type option." % filename)
                sys.exit(1)
            print("Reading %s (%s)" % (filename, append_doc))
            # read model to be appended
            model2 = append_func(filename, min_extent=min_extent, format=options.append_format or options.format)
            if model2.sources:
                sources += model2.sources
                for src in model2.sources:
                    src.name = "M%d:%s" % (modelnum, src.name)
                # recompute 'r' attribute (unless --center is in effect, in which case it's going to be done anyway below)
                if options.refresh_r:
                    for src in model2.sources:
                        src.setAttribute('r', Coordinates.angular_dist_pos_angle(src.pos.ra, src.pos.dec, *model.fieldCenter())[0])
            print("Appended %d sources from %s (%s)" % (len(model2.sources), filename, append_doc))

    # apply center, if specified
    if options.center:
        print("Center of field set to %s" % options.center)
        model.setFieldCenter(*center_radec)

    # apply selection by tag
    if options.tags:
        tags = []
        for ot in options.tags:
            tags += ot.split(",")
        for tag in tags:
            sources = [src for src in sources if getattr(src, tag, False)]
        if not sources:
            print("No sources left after selection by tag (-T/--tag) has been applied.")
            sys.exit(0)
        print("Selection by tag (%s) reduces this to %d sources" % (", ".join(options.tags), len(sources)))

    # apply selection by NaN
    if options.remove_nans:
        sources = [src for src in sources if not any([math.isnan(x)
                                                      for x in (src.pos.ra, src.pos.dec, src.flux.I)])]
        if not sources:
            print("No sources left after applying --remove-nans.")
            sys.exit(0)
        print("Removing NaN positions and fluxes reduces this to %d sources" % len(sources))

    # remove sources
    if options.remove_source:
        import fnmatch

        remove_names = set()
        for patt in options.remove_source:
            if patt[0] == "'" and patt[-1] == "'":
                patt = patt[1:-1]
            match = fnmatch.filter([src.name for src in sources], patt.replace("\\", ""))
            remove_names.update(match)
            print("Removing sources: %s matches %s" % (patt, ",".join(sorted(match))))
        sources = [src for src in sources if src.name not in remove_names]

    # add brick
    if options.add_brick:
        for brickspec in options.add_brick:
            # get names, check for uniqueness
            try:
                ff = brickspec.split(':')
                srcname = ff[0]
                fitsfile = ff[1]
                pad = float(ff[2] or '1') if len(ff) > 2 else 1
                tags = ff[3:] if len(ff) > 3 else []
            except:
                parser.error("Invalid --add-brick setting %s" % brickspec)
            if [src.name for src in sources if src.name == srcname]:
                print("Error: model already contains a source named '%s'" % srcname)

            input_hdu = fits.open(fitsfile)[0]
            hdr = input_hdu.header
            max_flux = float(input_hdu.data.max())

            wcs, refpix, refsky, ra_axis, dec_axis = get_wcs_info(hdr)
            naxis = hdr["NAXIS"]
            if ra_axis is None or dec_axis is None:
                print("Can't find RA and/or DEC axis in this FITS image")
                sys.exit(1)

            ra0, dec0 = refsky[ra_axis], refsky[dec_axis]

            # convert reference pixel from degrees to rad
            ra0 *= DEG
            dec0 *= DEG
            # get NX, NY
            nx, ny = input_hdu.data.shape[-1:-3:-1]

            # get half-size of image
            sx = nx/2 * abs(hdr[f"CDELT{ra_axis+1}"])
            sy = nx/2 * abs(hdr[f"CDELT{dec_axis+1}"])
            print("Image half-size is %f,%f deg" % (sx, sy))
            sx *= DEG
            sy *= DEG

            from Tigger.Models import ModelClasses, SkyModel

            pos = ModelClasses.Position(ra0, dec0)
            flux = ModelClasses.Flux(max_flux)
            shape = ModelClasses.FITSImage(sx, sy, 0, fitsfile, nx, ny, pad=pad)
            source = SkyModel.Source(srcname, pos, flux, shape=shape)
            for tag in tags:
                source.setAttribute(tag, True)
            if not options.refresh_r:
                source.setAttribute('r', Coordinates.angular_dist_pos_angle(ra0, dec0, *model.fieldCenter())[0])
            sources.append(source)
            print("Adding FITS source %s (%s,pad=%f) with tags %s" % (srcname, fitsfile, pad, tags))

    # convert apparent flux to intrinsic using the NEWSTAR beam gain
    if options.newstar_app_to_int:
        nsrc = 0
        for src in sources:
            bg = getattr(src, 'newstar_beamgain', None)
            if getattr(src, 'flux_apparent', None) and bg is not None:
                src.setAttribute('Iapp', src.flux.I)
                for pol in 'IQUV':
                    if hasattr(src.flux, pol):
                        setattr(src.flux, pol, getattr(src.flux, pol) / bg)
                src.removeAttribute('flux_apparent')
                src.setAttribute('flux_intrinsic', True)
                nsrc += 1
        print("Converted NEWSTAR apparent to intrinsic flux for %d model sources" % nsrc)
        if len(sources) != nsrc:
            print("  (%d sources were skipped for whatever reason.)" % (len(model.sources) - nsrc))
    elif options.newstar_int_to_app:
        nsrc = 0
        for src in sources:
            bg = getattr(src, 'newstar_beamgain', None)
            if getattr(src, 'flux_intrinsic', None) and bg is not None:
                src.setAttribute('Iapp', src.flux.I * bg)
                for pol in 'IQUV':
                    if hasattr(src.flux, pol):
                        setattr(src.flux, pol, getattr(src.flux, pol) * bg)
                src.removeAttribute('flux_intrinsic')
                src.setAttribute('flux_apparent', True)
                nsrc += 1
        print("Converted NEWSTAR apparent to intrinsic flux for %d model sources" % nsrc)
        if len(sources) != nsrc:
            print("  (%d sources were skipped for whatever reason.)" % (len(model.sources) - nsrc))

    # set refrence frequency
    if options.ref_freq >= 0:
        model.setRefFreq(options.ref_freq * 1e+6)
        print("Setting reference frequency to %f MHz" % options.ref_freq)

    # recenter
    if options.recenter:
        print("Shifting model to new center %s" % options.recenter)
        ra0, dec0 = model.fieldCenter()
        field_center = ra1, dec1 = recenter_radec
        ddec = dec1 - dec0
        cosd0, sind0 = math.cos(ddec), math.sin(ddec)
        for src in sources:
            ra, dec = src.pos.ra, src.pos.dec
            x, y, z = math.cos(ra - ra0) * math.cos(dec), math.sin(ra - ra0) * math.cos(dec), math.sin(dec)
            x1 = cosd0 * x - sind0 * z
            y1 = y
            z1 = sind0 * x + cosd0 * z
            src.pos.ra = ra1 + (math.atan2(y1, x1) if (x1 or y1) else 0)
            src.pos.dec = math.asin(z1)
        # reset model center
        model.setFieldCenter(ra1, dec1)

    # recompute radial distance
    if options.refresh_r:
        print("Recomputing the 'r' attribute based on the field center")
        model.recomputeRadialDistance()


    # select
    def getTagValue(src, tag):
        """Helper function: looks for the given tag in the source, or in its sub-objects"""
        for obj in src, src.pos, src.flux, getattr(src, 'shape', None), getattr(src, 'spectrum', None):
            if obj is not None and hasattr(obj, tag):
                return getattr(obj, tag)
        return None


    for selstr, tag, predicate, value in selections:
        # get tag value
        srctag = [(src, getTagValue(src, tag)) for src in model.sources]
        sources = [src for src, tag in srctag if tag is not None and predicate(tag, value)]
        print("Selection '%s' leaves %d out of %d sources" % (selstr, len(sources), len(model.sources)))
        if len(sources) != len(model.sources):
            model.setSources(sources)

    # set PB expression and estimate apparent fluxes
    pb = options.primary_beam
    if pb == "refresh":
        pb = model.primaryBeam()
        if pb:
            print("Recalculating apparent fluxes")
        else:
            print("No primary beam expression in model, ignoring '--primary-beam refresh' option")
    if options.app_to_int or options.int_to_app:
        pb = pb or model.primaryBeam()
        if pb:
            print("Converting apparent fluxes to intrinsic" if options.app_to_int else "Converting intrinsic fluxes to apparent")
        else:
            print("No primary beam expression in model and no --primary-beam option given, cannot convert between apparent and intrinsic.")
            sys.exit(1)
    if pb:
        fitsBeam = False
        if pb.lower().endswith('.fits'):  # if pb is a FITS file, load interpolator
            fitsBeam = True

            # Following code is nicked from Cattery/Siamese/OMS/pybeams_fits.py
            CORRS_XY = "xx", "xy", "yx", "yy"
            CORRS_RL = "rr", "rl", "lr", "ll"
            REIM = "re", "im"
            REALIMAG = dict(re="real", im="imag")



            def make_beam_filename(filename_pattern, corr, reim):
                """Makes beam filename for the given correlation and real/imaginary component (one of "re" or "im")"""
                return Utils.substitute_pattern(filename_pattern,
                                                corr=corr.lower(), xy=corr.lower(), CORR=corr.upper(), XY=corr.upper(),
                                                reim=reim.lower(), REIM=reim.upper(), ReIm=reim.title(),
                                                realimag=REALIMAG[reim].lower(), REALIMAG=REALIMAG[reim].upper(),
                                                RealImag=REALIMAG[reim].title())


            """Makes beam interpolator node for the given filename pattern."""
            filename_real = []
            filename_imag = []
            # load beam interpolator
            import Cattery.Siamese.OMS.InterpolatedBeams as InterpolatedBeams

            vbs = []
            for icorr, corr in enumerate(CORRS_XY if options.linear_pol else CORRS_RL):
                if icorr in (1, 2):
                    print('  omitting %s beam due to --beam-diag' % corr)
                    vbs.append(0)
                else:
                    # make FITS images or nulls for real and imaginary part
                    filenames = [make_beam_filename(pb, corr, 're'), make_beam_filename(pb, corr, 'im')]
                    print('Loading FITS Beams', filenames[0], filenames[1])
                    vb = InterpolatedBeams.LMVoltageBeam(verbose=(options.verbose or 0) - 2, l_axis=options.fits_l_axis,
                                                         m_axis=options.fits_m_axis)
                    vb.read(*filenames)
                    vbs.append(vb)

            model.setPrimaryBeam(vbs)
            # get frequency
            # fq = model.refFreq() or 1.4e+9
            beamRefFreq = (options.beam_freq or 0) * 1e+6 or model.refFreq() or 1424500000.12
            print("Using FITS beams with reference frequency %f MHz" % (beamRefFreq * 1e-6))

        else:  # else, assume pb is an expession
            try:
                pbexp = eval('lambda r,fq:' + pb)
                dum = pbexp(0, 1e+9);  # evaluate at r=0 and 1 GHz as a test
                if not isinstance(dum, float):
                    raise TypeError("does not evaluate to a float")
            except Exception as exc:
                print("Bad primary beam expression '%s': %s" % (pb, str(exc)))
                sys.exit(1)
            model.setPrimaryBeam(pb)
            # get frequency
            # fq = model.refFreq() or 1.4e+9
            fq = (options.beam_freq or 0) * 1e+6 or model.refFreq() or 1424500000.12
            print("Using beam expression '%s' with reference frequency %f MHz" % (pb, fq * 1e-6))

        nsrc = 0
        # ensure that every source has an 'r' attribute
        if not options.refresh_r:
            for src in sources:
                if not hasattr(src, 'r'):
                    src.setAttribute('r',
                                     Coordinates.angular_dist_pos_angle(src.pos.ra, src.pos.dec, *model.fieldCenter())[
                                         0])
        # evaluate sources
        if not (options.app_to_int or options.int_to_app):
            for src in sources:
                r = getattr(src, 'r', None)
                if r is not None:
                    bg = pbexp(r, fq)
                    src.setAttribute('beamgain', bg)
                    src.setAttribute('Iapp', src.flux.I * bg)
                    nsrc += 1
            print("Applied primary beam expression to %d model sources" % nsrc)
        else:
            # precompute PAs if fitsBeams are used
            if fitsBeam:
                if options.pa_from_ms is not None:
                    ms_strings = options.pa_from_ms.split(",")
                    ms_field = []
                    if len(ms_strings) > 1:
                        for ms_string in ms_strings:
                            match = re.match("^(.*?)(:[0-9]+)?$", ms_string)
                            if match:
                                msname, field = match.group(1), int(match.group(2)[1:]) if match.group(2) else 0
                            else:
                                msname, field = options.pa_from_ms, 0
                            ms_field.append((msname, field))
                    else:
                        ms_string = ms_strings[0]
                        match = re.match("^(.*?)(:[0-9]+)?$", ms_string)
                        if match:
                            msname, field = match.group(1), int(match.group(2)[1:]) if match.group(2) else 0
                            if os.path.exists(msname + "/SUBMSS"):
                                ms_field = [(ms, field) for ms in glob.glob(msname + "/SUBMSS/*") if os.path.isdir(ms)]
                            else:
                                ms_field = [[msname, 0]]
                    from pyrap.tables import table

                    dm, dq = pyrap_dmdq()
                    pas = []
                    zenith = dm.direction('AZEL', '0deg', '90deg')
                    for ms, field in ms_field:
                        print("Getting PA range from MS %s, field %d" % (ms, field))
                        tab = table(ms)
                        antpos = table(tab.getkeyword("ANTENNA")).getcol("POSITION")
                        ra, dec = table(tab.getkeyword("FIELD")).getcol("PHASE_DIR", field, 1)[0][0]
                        # make position measure from antenna 0
                        pos0 = dm.position('itrf', *[dq.quantity(x, 'm') for x in antpos[0]])
                        dm.do_frame(pos0)
                        # make direction measure from field centre
                        fld = dm.direction('J2000', dq.quantity(ra, "rad"), dq.quantity(dec, "rad"))
                        tab = tab.query("FIELD_ID==%d" % field)
                        # get unique times
                        times = np.array(sorted(set(tab.getcol("TIME")[~tab.getcol("FLAG_ROW")])))
                        pa1 = [(dm.do_frame(dm.epoch("UTC", dq.quantity(t, "s"))) and dm.posangle(fld,
                                                                                                  zenith).get_value(
                            "rad")) for t in times]
                        pas += pa1
                        pa1 = np.array(pa1) / DEG
                        if options.enable_plots:
                            import pylab

                            pylab.plot((times - times[0]) / 3600, pa1)
                            pylab.xlabel("Time since beginning of observation, hours")
                            pylab.ylabel("PA, degrees")
                            pylab.savefig(os.path.basename(ms) + ".parangle.png")
                            print("Saved plot " + os.path.basename(ms) + ".parangle.png")
                        print("MS %s, PA range is %fdeg to %fdeg" % (ms, pa1[0], pa1[-1]))
                    # get lm's rotated through those ranges
                    pa_range = np.array(pas)
                elif options.pa_range is not None:
                    try:
                        ang0, ang1 = list(map(float, options.pa_range.split(",", 1)))
                    except:
                        parser.error("Incorrect --pa-range option. FROM,TO values expected.")
                    pa_range = np.arange(ang0, ang1 + 1, 1) * DEG
                elif options.pa is not None:
                    pa_range = options.pa * DEG
                else:
                    pa_range = None
                if options.verbose:
                    print("PA (deg):", " ".join(["%f" % (x / DEG) for x in pa_range]) if np.iterable(
                        pa_range) else pa_range)
            if options.enable_plots:
                import pylab

                pylab.figure()
            for src in sources:
                r = getattr(src, 'r', None)
                if r is not None:
                    if fitsBeam:
                        # this is where the interpolator is called to determine the beam gain
                        # AIPS Memo 27 Sin Projection
                        ra0, dec0 = model.fieldCenter()
                        # ra0 = sources[0].pos.ra
                        # dec0 = sources[0].pos.dec
                        l = math.cos(src.pos.dec) * math.sin(src.pos.ra - ra0)
                        m = math.sin(src.pos.dec) * math.cos(dec0) - math.cos(src.pos.dec) * math.sin(dec0) * math.cos(
                            src.pos.ra - ra0)

                        # rotate through (range of) PA value(s), if such option is supplied above
                        if pa_range is not None:
                            l, m = rotatelm(l, m, pa_range)

                        Jones2Mueller = Jones2Mueller_linear if options.linear_pol else Jones2Mueller_circular

                        jones = [vb.interpolate(l, m, freq=beamRefFreq) if vb else np.array(0) for vb in vbs]
                        # incorrect old-style Jones averaging
                        if options.beam_average_jones:
                            a, b, c, d = [j.mean() for j in jones]
                            mueller = Jones2Mueller(np.matrix([[a, b], [c, d]]))
                            if options.verbose > 1:
                                print("%s: jones11 mean %f std %f" % (src.name, abs(a), abs(jones[0]).std()))
                                print("%s: jones22 mean %f std %f" % (src.name, abs(d), abs(jones[3]).std()))
                            if options.enable_plots:
                                pylab.plot(abs(jones[0]), label="|J11| " + src.name)
                        # new-style averaging of Mueller matrix
                        else:
                            muellers = [Jones2Mueller(np.matrix([[a, b], [c, d]])) for a, b, c, d in
                                        np.broadcast(*jones)]
                            mueller = sum(muellers) / len(muellers)
                            if options.enable_plots:
                                pylab.plot([m[0, 0] for m in muellers], label='M11 ' + src.name)
                            if options.verbose > 1:
                                print("%s: jones11 mean %f std %f" % (
                                    src.name, abs(jones[0].mean()), abs(jones[0]).std()))
                                print("%s: jones22 mean %f std %f" % (
                                    src.name, abs(jones[3].mean()), abs(jones[3]).std()))
                                print("%s: mueller11 mean %f std %f" % (
                                    src.name, mueller[0, 0], np.std([m[0, 0] for m in muellers])))
                        bg = mueller[0, 0]
                        ##          OMS 6/7/2015: let's do full inversion now to correct all four polarizations
                        if options.app_to_int:
                            if options.beam_nopol:
                                mueller = 1 / bg
                            else:
                                mueller = np.linalg.inv(mueller)
                        else:
                            if options.beam_nopol:
                                mueller = bg
                        ##            #for now, ignore full Stokes and just use Stokes' I
                        #            src.setAttribute('beamgain',bg)
                        nobeam = (bg < options.beam_clip)
                        spi = freqgrid = spiBg = None
                        # if no beam gain at this position, set appropriate tag
                        if nobeam:
                            src.setAttribute('nobeam', True)
                            src.setAttribute('Iapp', src.flux.I)
                        else:
                            src.removeAttribute('nobeam')
                            src.setAttribute('beamgain', bg)
                            iquv0 = np.matrix([[getattr(src.flux, stokes, 0.)] for stokes in "IQUV"])
                            iquv = mueller * iquv0
                            if options.verbose > 1:
                                print("%s: from %s to %s" % (src.name, iquv0.T, iquv.T))
                            if options.app_to_int and hasattr(src.flux, "I"):
                                src.setAttribute("Iapp", src.flux.I)
                            for i, stokes in enumerate("IQUV"):
                                if hasattr(src.flux, stokes):
                                    setattr(src.flux, stokes, iquv[i, 0])
                                # add spectral index of position in the beam
                            src_spectrum = getattr(src, 'spectrum', None)
                            if options.beam_spi and (src_spectrum or options.force_beam_spi_wo_spectrum):
                                # determine spectral index by determining bg across the freqs (using only Stokes' I)
                                import scipy.optimize

                                bw = options.beam_spi * 1e+6 / 2
                                # make a frequency grid of 10 points across the band
                                # freqgrid = np.arange(beamRefFreq-bw,beamRefFreq+bw,bw/5)
                                freqgrid = np.arange(beamRefFreq - bw, beamRefFreq + bw * 1.01, bw / 5)
                                gxx = vbs[0].interpolate(l, m, freq=freqgrid, freqaxis=2)
                                gyy = vbs[3].interpolate(l, m, freq=freqgrid, freqaxis=2)
                                spiBg = (gxx * gxx.conj() + gyy * gyy.conj()).real
                                spiBg = spiBg[:, 0, :]
                                # power law fit
                                logbg1 = np.log10(spiBg)
                                logbg = np.log10(spiBg.mean(axis=0))
                                logfreq = np.log10(freqgrid)
                                fitfunc = lambda p, x: p[0] + p[1] * x
                                errfunc = lambda p, x, y: (y - fitfunc(p, x))
                                pinit = [10 ** logbg[0], 0.]
                                if np.isinf(logbg).sum() > 0:
                                    spi = 0.
                                    amp0 = spiBg[0, 0]
                                else:
                                    out = scipy.optimize.leastsq(errfunc, pinit, args=(logfreq, logbg))
                                    spi = out[0][1]
                                    amp0 = 10. ** out[0][0]

                                # look for Spectral Index in spi attribute
                                # if no spectrum: add a SpectralIndex class to the source
                                # else: add spectral index from PB to SI (int-to-app), subtract (app-to-int)
                                if src_spectrum is None:
                                    setattr(src, 'spectrum', Tigger.Models.ModelClasses.SpectralIndex(spi, beamRefFreq))
                                else:
                                    ispiVal = getattr(src_spectrum, 'spi', None)
                                    setattr(src, 'spectrum', Tigger.Models.ModelClasses.SpectralIndex(ispiVal - spi,
                                                                                                      beamRefFreq) if options.app_to_int else Tigger.Models.ModelClasses.SpectralIndex(
                                        ispiVal + spi, beamRefFreq))

                        if options.verbose:
                            print(("%s: beamgain" % src.name), bg, "spi", spi, "clipped" if nobeam else "")
                        #            if spiBg is not None:
                        #              print src.name,repr(freqgrid),repr(spiBg.mean(0))

                    else:
                        bg = pbexp(r, fq)
                        src.setAttribute('beamgain', bg)
                        if hasattr(src.flux, 'I'):
                            src.setAttribute('Iapp', src.flux.I if options.app_to_int else src.flux.I * bg)
                        for stokes in "IQUV":
                            x = getattr(src.flux, stokes, None)
                            if x is not None:
                                setattr(src.flux, stokes, x / bg if options.app_to_int else x * bg)
                    nsrc += 1
            if options.enable_plots:
                pylab.legend()
                pylab.savefig("beamgains.png")
                print("Saved plot beamgains.png")
            print("Converted between apparent/intrinsic flux for %d model sources" % nsrc)
        if len(model.sources) != nsrc:
            print("  (%d sources were skipped for whatever reason, probably they didn't have an 'r' attribute)" % (
                    len(model.sources) - nsrc))

    # rename using COPART
    if options.rename:
        print("Renaming sources using the COPART convention")
        typecodes = dict(Gau="G", FITS="F")
        # sort sources by decreasing flux
        from past.builtins import cmp
        from functools import cmp_to_key
        sources = sorted(sources, key=cmp_to_key(lambda a, b: cmp(b.brightness(), a.brightness())))
        projection = Coordinates.Projection.SinWCS(*model.fieldCenter())
        # work out source clusters
        l = np.zeros(len(sources), float)
        m = np.zeros(len(sources), float)
        for i, src in enumerate(sources):
            l[i], m[i] = projection.lm(src.pos.ra, src.pos.dec)
        if options.cluster_dist:
            # now, convert to dist[i,j]: distance between sources i and j
            dist = np.sqrt(
                (l[:, np.newaxis] - l[np.newaxis, :]) ** 2 + (m[:, np.newaxis] - m[np.newaxis, :]) ** 2)
            # cluster[i] is (N,R), where N is cluster number for source #i, and R is rank of that source in the cluster
            # place source 0 into cluster 0,#0
            cluster = [(0, 0)]
            clustersize = [1]
            clusterflux = [sources[0].brightness()]
            dist0 = options.cluster_dist * DEG / 3600
            for i in range(1, len(sources)):
                src = sources[i]
                # find closest brighter source, and assign to its cluster if close enough
                imin = dist[i, :i].argmin()
                if dist[i, imin] <= dist0:
                    iclust, rank = cluster[imin]
                    cluster.append((iclust, clustersize[iclust]))
                    clustersize[iclust] += 1
                    clusterflux[iclust] += src.brightness()
                # else start new cluster from source
                else:
                    cluster.append((len(clustersize), 0))
                    clustersize.append(1)
                    clusterflux.append(src.brightness())
        else:
            cluster = [(i, 0) for i, src in enumerate(sources)]
        # now go over and rename the sources
        # make array of source names
        chars = [chr(x) for x in range(ord('a'), ord('z') + 1)]
        names = morenames = list(chars)
        while len(names) < len(sources):
            morenames = [ch + name for ch in chars for name in morenames]
            names += morenames
        # make a second version where the single-char names are capitalized
        Names = list(names)
        Names[:26] = [n.upper() for n in chars]
        # now go over and rename the sources
        clustername = {}
        for i, src in enumerate(sources):
            iclust, rank = cluster[i]
            # for up name of cluster based on rank-0 source
            if not rank:
                # lookup radius, in units of arcmin
                rad_min = math.sqrt(l[i] ** 2 + m[i] ** 2) * (60 / DEG)
                # divide by radial step
                rad = min(int(rad_min / options.radial_step), 10)
                radchr = '0123456789x'[rad]
                if rad_min > options.radial_step * 0.01:
                    # convert p.a. to tens of degrees
                    pa = math.atan2(l[i], m[i])
                    if pa < 0:
                        pa += math.pi * 2
                    pa = round(pa / (DEG * 10)) % 36
                    # make clustername
                    clusname = clustername[iclust] = "%s%02d%s" % (Names[iclust], pa, radchr)
                else:
                    clusname = clustername[iclust] = "%s0" % (Names[iclust])
                src.name = "%s%s" % (clusname, typecodes.get(src.typecode, ''))
                if options.cluster_dist:
                    src.setAttribute('cluster_lead', True)
            else:
                clusname = clustername[iclust]
                src.name = "%s%s%s" % (clusname, names[rank - 1], typecodes.get(src.typecode, ''))
            if options.cluster_dist:
                src.setAttribute('cluster', clusname)
                src.setAttribute('cluster_size', clustersize[iclust])
                src.setAttribute('cluster_flux', clusterflux[iclust])
    # check for duplicate names (if renaming, duplicate names cannot happen anyway, unless the naming algorithm above is broken)
    else:
        names = dict()
        sources0 = sources
        sources = []
        for i, src in enumerate(sources0):
            if src.name in names:
                print("Duplicate source '%s' at #%d (first found at #%d), removing" % (src.name, i, names[src.name]))
            else:
                names[src.name] = i
                sources.append(src)
    # assign prefix to source names
    if options.prefix:
        print("Prefixing source names with '%s'" % options.prefix)
        for src in sources:
            src.name = options.prefix + src.name
    # merge clusters
    if options.merge_clusters:
        tags = set(options.merge_clusters.split(',')) if options.merge_clusters != "ALL" else None
        # build up dict of clusters
        clusters = dict()
        for src in sources:
            clusname = getattr(src, 'cluster', '')
            clusters.setdefault(clusname, {})[src.name] = src
        # unclustered sources copied over as-is
        new_sources = list(clusters.pop('', {}).values())
        # next, deal with each cluster
        for clusname, srcdict in clusters.items():
            # leading source has the same name as the cluster
            src0 = srcdict.get(clusname)
            # if no leading source, or leading source not tagged, or length 1, then copy cluster as-is
            if not src0 or len(srcdict) < 2 or (tags is not None and
                                                not any([getattr(src0, tag, None) for tag in tags])):
                new_sources += list(srcdict.values())
            else:
                # sum fluxes
                for x in 'IQUV':
                    if hasattr(src0.flux, x):
                        setattr(src0.flux, x, sum([getattr(s.flux, x, 0) for s in srcdict.values()]))
                if hasattr(src0, 'Iapp'):
                    src0.Iapp = sum([getattr(s, 'Iapp', 0) for s in srcdict.values()])
                new_sources.append(src0)
                print("Merged cluster %s (%d sources)" % (src0.name, len(srcdict)))
        sources = new_sources
        model.setSources(sources)
    # save output
    print("Saving model containing %d sources to %s (%s)" % (len(sources), output, output_doc))
    export_func(model, output, sources=sources, format=options.output_format or None)
