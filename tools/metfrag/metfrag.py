from __future__ import absolute_import, print_function

import ConfigParser
import argparse
import csv
import glob
import multiprocessing
import os
import re
import shutil
import sys
import tempfile
from collections import defaultdict

import six

print(sys.version)

parser = argparse.ArgumentParser()
parser.add_argument('--input_pth')
parser.add_argument('--result_pth', default='metfrag_result.csv')

parser.add_argument('--temp_dir')
parser.add_argument('--polarity', default='pos')
parser.add_argument('--minMSMSpeaks', default=1)

parser.add_argument('--MetFragDatabaseType', default='PubChem')
parser.add_argument('--LocalDatabasePath', default='')
parser.add_argument('--LocalMetChemDatabaseServerIp', default='')

parser.add_argument('--DatabaseSearchRelativeMassDeviation', default=5)
parser.add_argument('--FragmentPeakMatchRelativeMassDeviation', default=10)
parser.add_argument('--FragmentPeakMatchAbsoluteMassDeviation', default=0.001)
parser.add_argument('--NumberThreads', default=1)
parser.add_argument('--UnconnectedCompoundFilter', action='store_true')
parser.add_argument('--IsotopeFilter', action='store_true')

parser.add_argument('--FilterMinimumElements', default='')
parser.add_argument('--FilterMaximumElements', default='')
parser.add_argument('--FilterSmartsInclusionList', default='')
parser.add_argument('--FilterSmartsExclusionList', default='')
parser.add_argument('--FilterIncludedElements', default='')
parser.add_argument('--FilterExcludedElements', default='')
parser.add_argument('--FilterIncludedExclusiveElements', default='')

parser.add_argument('--score_thrshld', default=0)
parser.add_argument('--pctexplpeak_thrshld', default=0)
parser.add_argument('--schema')
parser.add_argument('--cores_top_level', default=1)
parser.add_argument('--chunks', default=1)
parser.add_argument('--meta_select_col', default='name')
parser.add_argument('--skip_invalid_adducts', action='store_true')

parser.add_argument('--ScoreSuspectLists', default='')
parser.add_argument('--MetFragScoreTypes',
                    default="FragmenterScore,OfflineMetFusionScore")
parser.add_argument('--MetFragScoreWeights', default="1.0,1.0")

args = parser.parse_args()
print(args)

config = ConfigParser.ConfigParser()
config.read(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.ini'))

if os.stat(args.input_pth).st_size == 0:
    print('Input file empty')
    exit()

# Create temporary working directory
if args.temp_dir:
    wd = args.temp_dir
else:
    wd = tempfile.mkdtemp()

if os.path.exists(wd):
    shutil.rmtree(wd)
    os.makedirs(wd)
else:
    os.makedirs(wd)

######################################################################
# Setup regular expressions for MSP parsing dictionary
######################################################################
regex_msp = {}
regex_msp['name'] = [r'^Name(?:=|:)(.*)$']
regex_msp['polarity'] = [r'^ion.*mode(?:=|:)(.*)$',
                         r'^ionization.*mode(?:=|:)(.*)$',
                         r'^polarity(?:=|:)(.*)$']
regex_msp['precursor_mz'] = [r'^precursor.*m/z(?:=|:)\s*(\d*[.,]?\d*)$',
                             r'^precursor.*mz(?:=|:)\s*(\d*[.,]?\d*)$']
regex_msp['precursor_type'] = [r'^precursor.*type(?:=|:)(.*)$',
                               r'^adduct(?:=|:)(.*)$',
                               r'^ADDUCTIONNAME(?:=|:)(.*)$']
regex_msp['num_peaks'] = [r'^Num.*Peaks(?:=|:)\s*(\d*)$']
regex_msp['msp'] = [r'^Name(?:=|:)(.*)$']  # Flag for standard MSP format

regex_massbank = {}
regex_massbank['name'] = [r'^RECORD_TITLE:(.*)$']
regex_massbank['polarity'] = [r'^AC\$MASS_SPECTROMETRY:\s+ION_MODE\s+(.*)$']
regex_massbank['precursor_mz'] = [
    r'^MS\$FOCUSED_ION:\s+PRECURSOR_M/Z\s+(\d*[.,]?\d*)$']
regex_massbank['precursor_type'] = [
    r'^MS\$FOCUSED_ION:\s+PRECURSOR_TYPE\s+(.*)$']
regex_massbank['num_peaks'] = [r'^PK\$NUM_PEAK:\s+(\d*)']
regex_massbank['cols'] = [r'^PK\$PEAK:\s+(.*)']
regex_massbank['massbank'] = [
    r'^RECORD_TITLE:(.*)$']  # Flag for massbank format

if args.schema == 'msp':
    meta_regex = regex_msp
elif args.schema == 'massbank':
    meta_regex = regex_massbank
elif args.schema == 'auto':
    # If auto we just check for all the available paramter names and then
    # determine if Massbank or MSP based on the name parameter
    meta_regex = {}
    meta_regex.update(regex_massbank)
    meta_regex['name'].extend(regex_msp['name'])
    meta_regex['polarity'].extend(regex_msp['polarity'])
    meta_regex['precursor_mz'].extend(regex_msp['precursor_mz'])
    meta_regex['precursor_type'].extend(regex_msp['precursor_type'])
    meta_regex['num_peaks'].extend(regex_msp['num_peaks'])
    meta_regex['msp'] = regex_msp['msp']
else:
    sys.exit("No schema selected")

adduct_types = {
    '[M+H]+': 1.007276,
    '[M+NH4]+': 18.034374,
    '[M+Na]+': 22.989218,
    '[M+K]+': 38.963158,
    '[M+CH3OH+H]+': 33.033489,
    '[M+ACN+H]+': 42.033823,
    '[M+ACN+Na]+': 64.015765,
    '[M+2ACN+H]+': 83.06037,
    '[M-H]-': -1.007276,
    '[M+Cl]-': 34.969402,
    '[M+HCOO]-': 44.99819,
    '[M-H+HCOOH]-': 44.99819,
    # same as above but different style of writing adduct
    '[M+CH3COO]-': 59.01385,
    '[M-H+CH3COOH]-': 59.01385
    # same as above but different style of writing adduct
}
inv_adduct_types = {int(round(v, 0)): k for k, v in adduct_types.iteritems()}


# function to extract the meta data using the regular expressions
def parse_meta(meta_regex, meta_info=None):
    if meta_info is None:
        meta_info = {}
    for k, regexes in six.iteritems(meta_regex):
        for reg in regexes:
            m = re.search(reg, line, re.IGNORECASE)
            if m:
                meta_info[k] = '-'.join(m.groups()).strip()
    return meta_info


######################################################################
# Setup parameter dictionary
######################################################################
def init_paramd(args):
    paramd = defaultdict()

    paramd["MetFragDatabaseType"] = args.MetFragDatabaseType

    if args.MetFragDatabaseType == "LocalCSV":
        paramd["LocalDatabasePath"] = args.LocalDatabasePath
    elif args.MetFragDatabaseType == "MetChem":
        paramd["LocalMetChemDatabase"] = \
            config.get('MetChem', 'LocalMetChemDatabase')
        paramd["LocalMetChemDatabasePortNumber"] = \
            config.get('MetChem', 'LocalMetChemDatabasePortNumber')
        paramd["LocalMetChemDatabaseServerIp"] =  \
            args.LocalMetChemDatabaseServerIp
        paramd["LocalMetChemDatabaseUser"] = \
            config.get('MetChem', 'LocalMetChemDatabaseUser')
        paramd["LocalMetChemDatabasePassword"] = \
            config.get('MetChem', 'LocalMetChemDatabasePassword')

    paramd["FragmentPeakMatchAbsoluteMassDeviation"] =  \
        args.FragmentPeakMatchAbsoluteMassDeviation
    paramd["FragmentPeakMatchRelativeMassDeviation"] =  \
        args.FragmentPeakMatchRelativeMassDeviation
    paramd["DatabaseSearchRelativeMassDeviation"] = \
        args.DatabaseSearchRelativeMassDeviation
    paramd["SampleName"] = ''
    paramd["ResultsPath"] = os.path.join(wd)

    if args.polarity == "pos":
        paramd["IsPositiveIonMode"] = True
        paramd["PrecursorIonModeDefault"] = "1"
        paramd["PrecursorIonMode"] = "1"
        paramd["nm_mass_diff_default"] = 1.007276
    else:
        paramd["IsPositiveIonMode"] = False
        paramd["PrecursorIonModeDefault"] = "-1"
        paramd["PrecursorIonMode"] = "-1"
        paramd["nm_mass_diff_default"] = -1.007276

    paramd["MetFragCandidateWriter"] = "CSV"
    paramd["NumberThreads"] = args.NumberThreads

    if args.ScoreSuspectLists:
        paramd["ScoreSuspectLists"] = args.ScoreSuspectLists

    paramd["MetFragScoreTypes"] = args.MetFragScoreTypes
    paramd["MetFragScoreWeights"] = args.MetFragScoreWeights

    dct_filter = defaultdict()
    filterh = []

    if args.UnconnectedCompoundFilter:
        filterh.append('UnconnectedCompoundFilter')

    if args.IsotopeFilter:
        filterh.append('IsotopeFilter')

    if args.FilterMinimumElements:
        filterh.append('MinimumElementsFilter')
        dct_filter['FilterMinimumElements'] = args.FilterMinimumElements

    if args.FilterMaximumElements:
        filterh.append('MaximumElementsFilter')
        dct_filter['FilterMaximumElements'] = args.FilterMaximumElements

    if args.FilterSmartsInclusionList:
        filterh.append('SmartsSubstructureInclusionFilter')
        dct_filter[
            'FilterSmartsInclusionList'] = args.FilterSmartsInclusionList

    if args.FilterSmartsExclusionList:
        filterh.append('SmartsSubstructureExclusionFilter')
        dct_filter[
            'FilterSmartsExclusionList'] = args.FilterSmartsExclusionList

    # My understanding is that both 'ElementInclusionExclusiveFilter'
    # and 'ElementExclusionFilter' use 'FilterIncludedElements'
    if args.FilterIncludedExclusiveElements:
        filterh.append('ElementInclusionExclusiveFilter')
        dct_filter[
            'FilterIncludedElements'] = args.FilterIncludedExclusiveElements

    if args.FilterIncludedElements:
        filterh.append('ElementInclusionFilter')
        dct_filter['FilterIncludedElements'] = args.FilterIncludedElements

    if args.FilterExcludedElements:
        filterh.append('ElementExclusionFilter')
        dct_filter['FilterExcludedElements'] = args.FilterExcludedElements

    if filterh:
        fcmds = ','.join(filterh) + ' '
        for k, v in six.iteritems(dct_filter):
            fcmds += "{0}={1} ".format(str(k), str(v))

        paramd["MetFragPreProcessingCandidateFilter"] = fcmds

    return paramd


######################################################################
# Function to run metfrag when all metainfo and peaks have been parsed
######################################################################
def run_metfrag(meta_info, peaklist, args, wd, spectrac, adduct_types):
    # Get sample details (if possible to extract) e.g. if created as part of
    # the msPurity pipeline) choose between getting additional details to add
    # as columns as either all meta data from msp, just details from the
    # record name (i.e. when using msPurity and we have the columns coded into
    # the name) or just the spectra index (spectrac)].
    # Returns the parameters used and the command line call

    paramd = init_paramd(args)
    if args.meta_select_col == 'name':
        # have additional column of just the name
        paramd['additional_details'] = {'name': meta_info['name']}
    elif args.meta_select_col == 'name_split':
        # have additional columns split by "|" and
        # then on ":" e.g. MZ:100.2 | RT:20 | xcms_grp_id:1
        paramd['additional_details'] = {
            sm.split(":")[0].strip(): sm.split(":")[1].strip() for sm in
            meta_info['name'].split("|")}
    elif args.meta_select_col == 'all':
        # have additional columns based on all the meta information
        # extracted from the MSP
        paramd['additional_details'] = meta_info
    else:
        # Just have and index of the spectra in the MSP file
        paramd['additional_details'] = {'spectra_idx': spectrac}

    paramd["SampleName"] = "{}_metfrag_result".format(spectrac)

    # =============== Output peaks to txt file  ==============================
    paramd["PeakListPath"] = os.path.join(wd,
                                          "{}_tmpspec.txt".format(spectrac))

    # write spec file
    with open(paramd["PeakListPath"], 'w') as outfile:
        for p in peaklist:
            outfile.write(p[0] + "\t" + p[1] + "\n")

    # =============== Update param based on MSP metadata ======================
    # Replace param details with details from MSP if required
    if 'precursor_type' in meta_info and \
            meta_info['precursor_type'] in adduct_types:
        adduct = meta_info['precursor_type']
        nm = float(meta_info['precursor_mz']) - adduct_types[
            meta_info['precursor_type']]
        paramd["PrecursorIonMode"] = \
            int(round(adduct_types[meta_info['precursor_type']], 0))
    elif not args.skip_invalid_adducts:
        adduct = inv_adduct_types[int(paramd['PrecursorIonModeDefault'])]
        paramd["PrecursorIonMode"] = paramd['PrecursorIonModeDefault']
        nm = float(meta_info['precursor_mz']) - paramd['nm_mass_diff_default']
    else:
        print('Skipping {}'.format(paramd["SampleName"]))
        return '', ''

    paramd['additional_details']['adduct'] = adduct
    paramd["NeutralPrecursorMass"] = nm

    # ============== Create CLI cmd for metfrag ===============================
    cmd = "metfrag"
    for k, v in six.iteritems(paramd):
        if k not in ['PrecursorIonModeDefault', 'nm_mass_diff_default',
                     'additional_details']:
            cmd += " {}={}".format(str(k), str(v))

    # ============== Run metfrag ==============================================
    # print(cmd)
    # Filter before process with a minimum number of MS/MS peaks
    if plinesread >= float(args.minMSMSpeaks):

        if int(args.cores_top_level) == 1:
            os.system(cmd)

    return paramd, cmd


def work(cmds):
    return [os.system(cmd) for cmd in cmds]


######################################################################
# Parse MSP file and run metfrag CLI
######################################################################
# keep list of commands if performing in CLI in parallel
cmds = []
# keep a dictionary of all params
paramds = {}
# keep count of spectra (for uid)
spectrac = 0
# this dictionary will store the meta data results form the MSp file
meta_info = {}

with open(args.input_pth, "r") as infile:
    # number of lines for the peaks
    pnumlines = 0
    # number of lines read for the peaks
    plinesread = 0
    for line in infile:
        line = line.strip()

        if pnumlines == 0:
            # ============== Extract metadata from MSP ========================
            meta_info = parse_meta(meta_regex, meta_info)

            if ('massbank' in meta_info and 'cols' in meta_info) or (
                    'msp' in meta_info and 'num_peaks' in meta_info):
                pnumlines = int(meta_info['num_peaks'])
                plinesread = 0
                peaklist = []

        elif plinesread < pnumlines:
            # ============== Extract peaks from MSP ==========================
            # .split() will split on any empty space (i.e. tab and space)
            line = tuple(line.split())
            # Keep only m/z and intensity, not relative intensity
            save_line = tuple(line[0].split() + line[1].split())
            plinesread += 1
            peaklist.append(save_line)

        elif plinesread and plinesread == pnumlines:
            # ======= Get sample name and additional details for output =======
            spectrac += 1
            paramd, cmd = run_metfrag(meta_info, peaklist, args, wd, spectrac,
                                      adduct_types)

            if paramd:
                paramds[paramd["SampleName"]] = paramd
                cmds.append(cmd)

            meta_info = {}
            pnumlines = 0
            plinesread = 0

    # end of file. Check if there is a MSP spectra to run metfrag on still
    if plinesread and plinesread == pnumlines:

        paramd, cmd = run_metfrag(meta_info, peaklist, args, wd, spectrac + 1,
                                  adduct_types)

        if paramd:
            paramds[paramd["SampleName"]] = paramd
            cmds.append(cmd)

# Perform multiprocessing on command line call level
if int(args.cores_top_level) > 1:
    cmds_chunks = [cmds[x:x + int(args.chunks)] for x in
                   list(range(0, len(cmds), int(args.chunks)))]
    pool = multiprocessing.Pool(processes=int(args.cores_top_level))
    pool.map(work, cmds_chunks)
    pool.close()
    pool.join()

######################################################################
# Concatenate and filter the output
######################################################################
# outputs might have different headers. Need to get a list of all the
# headers before we start merging the files
# outfiles = [os.path.join(wd, f) for f in glob.glob(os.path.join(wd,
# "*_metfrag_result.csv"))]
outfiles = glob.glob(os.path.join(wd, "*_metfrag_result.csv"))

if len(outfiles) == 0:
    print('No results')
    sys.exit()

headers = []
c = 0
for fn in outfiles:
    with open(fn, 'r') as infile:
        reader = csv.reader(infile)
        if sys.version_info >= (3, 0):
            headers.extend(next(reader))
        else:
            headers.extend(reader.next())
        # check if file has any data rows
        for i, row in enumerate(reader):
            c += 1
            if i == 1:
                break

# if no data rows (e.g. matches) then do not save an
# output and leave the program
if c == 0:
    print('No results')
    sys.exit()

additional_detail_headers = ['sample_name']
for k, paramd in six.iteritems(paramds):
    additional_detail_headers = list(set(
        additional_detail_headers + list(paramd['additional_details'].keys())))

# add inchikey if not already present (missing in metchem output)
if 'InChIKey' not in headers:
    headers.append('InChIKey')

headers = additional_detail_headers + sorted(list(set(headers)))

# Sort files nicely
outfiles.sort(
    key=lambda s: int(re.match(r'^.*/(\d+)_metfrag_result.csv', s).group(1)))

print(outfiles)

# merge outputs
with open(args.result_pth, 'a') as merged_outfile:
    dwriter = csv.DictWriter(merged_outfile, fieldnames=headers,
                             delimiter='\t', quotechar='"')
    dwriter.writeheader()

    for fn in outfiles:

        with open(fn) as infile:
            reader = csv.DictReader(infile, delimiter=',', quotechar='"')
            for line in reader:
                bewrite = True
                for key, value in line.items():
                    # Filter when no MS/MS peak matched
                    
                    if key == "ExplPeaks":
                        if float(args.pctexplpeak_thrshld) > 0 \
                                and value and "NA" in value:
                            bewrite = False
                    # Filter with a score threshold
                    elif key == "Score":
                        if value and float(value) <= float(args.score_thrshld):
                            bewrite = False
                    elif key == "NoExplPeaks":
                        nbfindpeak = float(value)
                    elif key == "NumberPeaksUsed":
                        totpeaks = float(value)
                # Filter with a relative number of peak matched
                try:
                    pctexplpeak = nbfindpeak / totpeaks * 100
                except ZeroDivisionError:
                    bewrite = False
                else:
                    if pctexplpeak < float(args.pctexplpeak_thrshld):
                        bewrite = False

                # Write the line if it pass all filters
                if bewrite:
                    bfn = os.path.basename(fn)
                    bfn = bfn.replace(".csv", "")
                    line['sample_name'] = paramds[bfn]['SampleName']
                    ad = paramds[bfn]['additional_details']

                    if args.MetFragDatabaseType == "MetChem":
                        # for some reason the metchem database option does
                        # not report the full inchikey (at least in the Bham
                        # setup. This ensures we always get the fully inchikey
                        line['InChIKey'] = '{}-{}-{}'.format(line['InChIKey1'],
                                                             line['InChIKey2'],
                                                             line['InChIKey3'])

                    line.update(ad)
                    dwriter.writerow(line)
