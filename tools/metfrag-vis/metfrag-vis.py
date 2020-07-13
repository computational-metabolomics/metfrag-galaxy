#!/usr/bin/env python

# Load modules
import argparse
import base64
import csv
import os
import re
import time
import urllib.parse

import matplotlib.pyplot as plt

import pubchempy

import requests

# Parse arguments
parser = argparse.ArgumentParser(
    description='Visualise MetFrag results in html.')
parser.add_argument('-v', '--version', action='version',
                    version='MetFrag-vis Version 0.9',
                    help='show version')
parser.add_argument('-i', '--input', metavar='metfrag_results.tsv',
                    dest="input_tsv", required=True,
                    help='MetFrag results as input')
parser.add_argument('-o', '--output', metavar='metfrag_results.html',
                    dest="output_html", required=True,
                    help='Write MetFrag results into this output file')
parser.add_argument('-m', '--max-candidates', metavar='10',
                    dest="max_candidates", default=10, type=int,
                    required=False,
                    help='Maximum number of candidates per compound [1-1000]')
parser.add_argument('-s', '--synonyms', dest='synonyms', action='store_true',
                    required=False,
                    help='Fetch synonyms from PubChem [disabled by default]')
parser.add_argument('-c', '--classyfire', dest='classyfire',
                    action='store_true', required=False,
                    help='Fetch compound classes from ClassyFire'
                         ' [disabled by default]')

args = parser.parse_args()

# Input CSV with MetFrag results
input_tsv = args.input_tsv

# Output html of MetFrag results
output_html = args.output_html

# Max number of candidates per compound
max_candidates = args.max_candidates

# PubChem synonyms
pubchem_synonyms_enabled = args.synonyms

# ClassyFire classes
classyfire_classes_enabled = args.classyfire


# ---------- cdk_inchi_to_svg ----------
def cdk_inchi_to_svg(inchi):
    if "cdk-inchi-to-svg" in os.environ:
        JAVA_CMD = 'cdk-inchi-to-svg' + ' ' + str(
            '\'' + inchi + '\'') + ' ' + 'cdk-inchi-to-svg-output.svg'
    else:
        JAVA_BINARY = '/usr/local/bin/java'
        CDK_INCHI_TO_SVG_JAR = '/usr/local/bin/' \
                               'cdk-inchi-to-svg-0.0.1-' \
                               'SNAPSHOT-jar-with-dependencies.jar'
        JAVA_CMD = str(
            JAVA_BINARY + ' ' + '-jar' + ' ' + CDK_INCHI_TO_SVG_JAR + ' '
            + str('\'' + inchi + '\'') + ' ' + 'cdk-inchi-to-svg-output.svg')

    # Exec cdk-inchi-to-svg JAVA binary
    exitcode = os.system(JAVA_CMD)

    # Check whether binary has successfully been run
    if (exitcode == 0):
        with open("cdk-inchi-to-svg-output.svg", "r") as svg_file:
            svg_string = []
            for line in svg_file:
                if not ('<?xml' in line) and not ('<!DOCTYPE' in line):
                    if (' fill=\'#FFFFFF\'' in line):
                        line = re.sub(' fill=\'#FFFFFF\'',
                                      ' fill=\'#FFFFFF\' fill-opacity=\'0.0\'',
                                      line)
                    svg_string.append(line)
        svg_file.close()
        os.remove("cdk-inchi-to-svg-output.svg")
        return (str(''.join(svg_string)))
    else:
        return ('&nbsp;')


# ---------- pubchem_link ----------
def pubchem_link(compound_name):
    return (str('https://pubchem.ncbi.nlm.nih.gov/#query=' + compound_name))


# ---------- kegg_link ----------
def kegg_link(compound_name):
    return (str(
        'https://www.genome.jp/dbget-bin/'
        'www_bfind_sub?mode=bfind&max_hit=1000&dbkey=kegg&keywords=' +
        compound_name))


# ---------- biocyc_link ----------
def biocyc_link(compound_name):
    biocyc_url = urllib.parse.urlparse(
        str(
            'https://www.biocyc.org/'
            'substring-search?type=NIL&object=' +
            compound_name + '&quickSearch=Quick+Search'))
    return (biocyc_url.geturl())


# ---------- hmdb_link ----------
def hmdb_link(compound_name):
    hmdb_url = urllib.parse.urlparse(
        str(
            'https://hmdb.ca/unearth/q?utf8=\xe2&query=' +
            compound_name + '&searcher=metabolites&button='))
    return (hmdb_url.geturl())


# ---------- hmdb_link ----------
def chebi_link(inchi):
    return (str(
        'https://www.ebi.ac.uk/chebi/advancedSearchFT.do?searchString='
        + inchi))


# ---------- PubChem Synonyms ----------
def fetch_pubchem_synonyms(inchi):
    if not ('InChI=' in inchi):
        return ('&nbsp;')

    # Fetch CID from InChI
    print('Retrieving PubChem CID from InChI...')
    compound = pubchempy.get_compounds(identifier=inchi, namespace='inchi')
    compound_cid = re.sub(r'\).*', '', re.sub(r'.*\(', '', str(compound)))
    if len(compound_cid) <= 1:
        print(str('Warning. No match for InChI \"' + str(inchi) + '\".'))
        return ('&nbsp;')

    # Retrieve compound
    print('Retrieving PubChem compound information...')
    compound = pubchempy.Compound.from_cid(compound_cid)
    if ('synonyms' in dir(compound)):
        return ('; '.join(compound.synonyms))
    else:
        print(str('Warning. No synonyms found for CID \"' + str(
            compound_cid) + '\".'))
        return ('&nbsp;')


# ---------- ClassyFire ----------
def fetch_classyfire_classes(inchi):
    if not ('InChI=' in inchi):
        return ('&nbsp;')

    # Send POST request to ClassyFire
    print('Sending request to ClassyFire...')
    classyfire_url = 'http://classyfire.wishartlab.com/queries.json'
    classyfire_post = str(
        '{\"label\":\"metfrag\",\"query_input\":\"' + inchi +
        '\",\"query_type\":\"STRUCTURE\"}')
    classyfire_headers = {'Content-Type': 'application/json'}
    classyfire_request = requests.post(classyfire_url, data=classyfire_post,
                                       headers=classyfire_headers)

    # Only continue when request has been successfully sent
    if (classyfire_request.status_code != 201):
        print('Error! Could not send request to ClassyFire. \"',
              str(classyfire_request.status_code) + ': ' + str(
                  classyfire_request.reason), '\". Skipping entry.')
        return ('&nbsp;')

    # Get ClassyFire Query ID
    classyfire_request.json()
    classyfire_query_id = classyfire_request.json()['id']

    # Query ClassyFire in max. 20 attempts
    classyfire_request_loop = 0
    while (classyfire_request_loop < 20):
        print(str(
            'Sending query ' + str(classyfire_query_id) + ' to ClassyFire...'))
        time.sleep(10)
        classyfire_query = requests.get(
            str('http://classyfire.wishartlab.com/queries/' + str(
                classyfire_query_id) + '.json'))

        if (classyfire_query.status_code == 200) and (
                classyfire_query.json()['classification_status'] == 'Done'):
            classyfire_request_loop = 999
            break
        else:
            classyfire_request_loop += 1

    if classyfire_request_loop == 999:
        # Direct parent
        direct_parent_name = classyfire_query.json()[
            'entities'][0]['direct_parent']['name']
        direct_parent_url = classyfire_query.json()[
            'entities'][0]['direct_parent']['url']
        direct_parent = str(
            '<a target="_blank" href="' + direct_parent_url + '">' +
            direct_parent_name + '</a>')

        # Alternative parents
        alt_parents = []
        for i in range(0, len(classyfire_query.json()['entities'][0][
                                  'alternative_parents'])):
            alt_parent_name = classyfire_query.json()[
                'entities'][0]['alternative_parents'][i]['name']
            alt_parent_url = classyfire_query.json()[
                'entities'][0]['alternative_parents'][i]['url']
            alt_parent = str(
                '<a target="_blank" href="' + alt_parent_url + '">' +
                alt_parent_name + '</a>')
            alt_parents.append(alt_parent)

        # Concat classes
        classes = str('<b>' + direct_parent + '</b>, <br>' + str(
            ', <br>'.join(alt_parents)))
    else:
        print('Warning. Timout sending query to ClassyFire. Skipping entry.')
        classes = '&nbsp;'

    return (classes)


# ---------- Plot Spectrum ----------
def plot_spectrum(spectrum, spectrum_explained, spectrum_explained_formulas):
    # Plot
    plt.figure(figsize=[5.5, 4.4])
    plt.xlabel('m/z')
    plt.ylabel('intensity')

    # Plot spectrum
    x = []
    y = []
    for i in spectrum.split(';'):
        t = i.split('_')
        x.append(t[0])
        y.append(t[1])

    for i in range(0, len(x)):
        plt.plot([float(x[i]), float(x[i])], [0, float(y[i])], linewidth=1,
                 color='black')
        plt.plot(float(x[i]), float(y[i]), 'o', color='black', markersize=4)

    if not (spectrum_explained == 'NA') and not (
            spectrum_explained_formulas == 'NA'):
        # Plot explained peaks
        ex = []
        ey = []
        for i in spectrum_explained.split(';'):
            t = i.split('_')
            ex.append(t[0])
            ey.append(y[x.index(t[0])])

        for i in range(0, len(ex)):
            plt.plot([float(ex[i]), float(ex[i])], [0, float(ey[i])],
                     linewidth=3, color='#2b8126')
            plt.plot(float(ex[i]), float(ey[i]), 'o', color='#2b8126',
                     markersize=8)

        # Plot formulas on explained peaks
        ex = []
        ey = []
        ez = []
        for i in spectrum_explained_formulas.split(';'):
            t = i.split(':')
            ex.append(t[0])
            ey.append(y[x.index(t[0])])
            ez.append(t[1])

        for i in range(0, len(ex)):
            plt.text(float(ex[i]), float(ey[i]) + 1000, ez[i], color='#2b8126',
                     fontsize=8,
                     horizontalalignment='center', verticalalignment='bottom')

    # Save SVG
    plt.savefig("metfrag-vis-spectrum.svg", format="svg", transparent=True)
    plt.close()

    # Import SVG
    with open("metfrag-vis-spectrum.svg", "r") as svg_file:
        svg_string = []
        for line in svg_file:
            if not ('<?xml' in line) and not ('<!DOCTYPE' in line) and not (
                    '  "http://www.w3.org/Graphics' in line):
                svg_string.append(line)
    svg_file.close()
    os.remove("metfrag-vis-spectrum.svg")
    return (str(''.join(svg_string)))


# #################### MAIN ####################
if pubchem_synonyms_enabled:
    print('Fetching of PubChem Synonyms enabled.')
if classyfire_classes_enabled:
    print('Fetching of ClassyFire Classes enabled.')

# Open output html file
try:
    metfrag_html = open(output_html, "w")
except Exception as e:
    print("Error writing output file. {}".format(e))
    exit(1)

# Write html header
metfrag_html.write('<!DOCTYPE html>\n')
metfrag_html.write('<html>\n')
metfrag_html.write('<head>\n')
metfrag_html.write('<title>' + 'msPurity MetFrag results' + '</title>\n')
metfrag_html.write('<style type="text/css">\n')
metfrag_html.write('svg { width: 200px; height: 100%; }\n')
metfrag_html.write(
    'body { font-family: Lucida, Verdana, Arial, Helvetica, sans-serif; '
    'font-size: 13px; text-align: left; '
    'color: #000000; margin: 8px 8px 8px 8px; }\n')
metfrag_html.write(
    'A { color: #2b8126; text-decoration: none; background: transparent; }\n')
metfrag_html.write(
    'A:visited { '
    'color: #19681a; text-decoration: none; background: transparent; '
    '}\n')
metfrag_html.write(
    'A:hover { '
    'color: #8fc180; text-decoration: underline; background: transparent; '
    '}\n')
metfrag_html.write(
    'h1 { font-size: 32px; font-weight: bold; text-align: center; '
    'padding: 0px 0px 4px 0px; margin: 26px 0px 0px 0px; }\n')
metfrag_html.write(
    'h2 { font-size: 24px; font-weight: bold; text-align: left; '
    'padding: 0px 0px 4px 0px; margin: 26px 0px 0px 0px; }\n')
metfrag_html.write(
    'table { font-family: Lucida, Verdana, Arial, Helvetica, sans-serif; '
    'font-size: 10px; text-align: left; '
    'line-height: 10px; border: 1px solid #e3efdf; '
    'background-color: #ecf5ea; margin-bottom: 8px; '
    'min-width: 1600px; max-width: 2400px; }\n')
metfrag_html.write(
    '#tablediv { width: 100%; min-width: 20px; max-width: 200px; }\n')
metfrag_html.write('.tdmax { min-width: 200px; max-width: 200px; }\n')
metfrag_html.write('.tdvar { min-width: 200px; max-width: 600px; }\n')
metfrag_html.write('tr:nth-child(even) { background-color: #f6faf5; }\n')
metfrag_html.write('</style>\n')
metfrag_html.write('</head>\n')
metfrag_html.write('<body>\n')

# Read input csv file
with open(input_tsv, "r") as metfrag_file:
    metfrag_results = csv.DictReader(metfrag_file, delimiter='\t')
    # Parse each line
    line_count = 0
    compound = ""
    candidates = 0
    for row in metfrag_results:

        # Start new document
        if (line_count == 0):
            if os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'metfrag_logo.png'):
                logo_pth = os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    'metfrag_logo.png')
            else:
                logo_pth = '/usr/local/share/metfrag/metfrag_logo.png'

            with open(logo_pth, "rb") as png_file:
                png_encoded = base64.b64encode(png_file.read())
            metfrag_html.write(str(
                '\n<h1><img style="vertical-align:bottom" '
                'src="data:image/png;base64,' +
                png_encoded.decode('utf-8') +
                '" alt="metfrag-logo" width="150"></img><text '
                'style="line-height:2.0">&nbsp;&nbsp;results</text></h1>\n'
            ))
        else:
            # Parameter list at beginning of document
            if (line_count == 1):
                metfrag_html.write('\n<h2>Parameter list</h2>\n')
                metfrag_html.write(str('MetFragDatabaseType=' +
                                       re.sub(' .*', '',
                                              re.sub('.*MetFragDatabaseType=',
                                                     '',
                                                     row[
                                                         "MetFragCLIString"]))
                                       + '<br>\n')
                                   )
                metfrag_html.write(str('PrecursorIonMode=' +
                                       re.sub(' .*', '',
                                              re.sub('.*PrecursorIonMode=', '',
                                                     row[
                                                         "MetFragCLIString"]))
                                       + '<br>\n')
                                   )
                metfrag_html.write(str('DatabaseSearchRelativeMassDeviation=' +
                                       re.sub(' .*', '',
                                              re.sub(
                                                  '.*DatabaseSearchRelative'
                                                  'MassDeviation=',
                                                  '',
                                                  row[
                                                      "MetFragCLIString"])) +
                                       '<br>\n')
                                   )
                metfrag_html.write(
                    str('FragmentPeakMatchAbsoluteMassDeviation=' +
                        re.sub(' .*', '',
                               re.sub(
                                   '.*FragmentPeakMatchAbsoluteMassDeviation=',
                                   '',
                                   row["MetFragCLIString"])) + '<br>\n')
                )
                metfrag_html.write(
                    str('FragmentPeakMatchRelativeMassDeviation=' +
                        re.sub(' .*', '',
                               re.sub(
                                   '.*FragmentPeakMatchRelativeMassDeviation=',
                                   '',
                                   row["MetFragCLIString"])) + '<br>\n')
                )
                metfrag_html.write(str('FilterExcludedElements=' +
                                       re.sub(' .*', '',
                                              re.sub(
                                                  '.*FilterExcludedElements=',
                                                  '',
                                                  row[
                                                      "MetFragCLIString"])) +
                                       '<br>\n')
                                   )
                metfrag_html.write(str('FilterIncludedElements=' +
                                       re.sub(' .*', '',
                                              re.sub(
                                                  '.*FilterIncludedElements=',
                                                  '',
                                                  row[
                                                      "MetFragCLIString"])) +
                                       '<br>\n')
                                   )
                metfrag_html.write(str('MetFragScoreTypes=' +
                                       re.sub(' .*', '',
                                              re.sub('.*MetFragScoreTypes=',
                                                     '',
                                                     row[
                                                         "MetFragCLIString"]))
                                       + '<br>\n')
                                   )
            # New compound in list
            if (row["name"] != compound):
                compound = row["name"]
                candidates = 0
                identifier = row["name"]
                monoisotopic_mass = row["MonoisotopicMass"]
                precursor_mz = row["precursor_mz"]

                if "retention_time" in row:
                    precursor_rt = row["retention_time"]
                    try:
                        precursor_rt = round(float(precursor_rt), 4)
                    except ValueError:
                        continue
                else:
                    precursor_rt = ''

                if "precursor_type" in row:
                    precursor_type = row["precursor_type"]
                elif "adduct" in row:
                    precursor_type = row["adduct"]
                else:
                    precursor_type = ''

                if line_count > 1:
                    metfrag_html.write(str('</table>\n'))

                metfrag_html.write(str('\n' + '<h2>' + identifier + '</h2>\n'))
                metfrag_html.write(str('<p><b>Precursor Type:</b> ' + str(
                    precursor_type) + '<br>'))
                metfrag_html.write(str('<b>Precursor Mass:</b> ' + str(
                    round(float(precursor_mz), 4)) + '<br>'))
                metfrag_html.write(
                    str('<b>Precursor Retention Time:</b> ' + str(
                        precursor_mz) + '<br></p>'))
                metfrag_html.write(str('\n' + '<table>\n'))
                metfrag_html.write(str(
                    '<tr style="vertical-align:bottom; '
                    'background-color:#e3efdf;">'
                    + '<td class="tdmax">' + '<b>Spectrum</b>' + '</td>'
                    + '<td class="tdmax">' + '<b>Structure</b>' + '</td>'
                    + '<td>' + '<b>Monoisotopic Mass</b>' + '</td>'
                    + '<td>' + '<b>Molecular Formula</b>' + '</td>'
                    + '<td>' + '<b>Compound Name</b>' + '</td>'
                    + '<td class="tdvar">' + '<b>PubChem Synonyms</b>'+'</td>'
                    + '<td>' + '<b>Compound Classes</b>' + '</td>'
                    + '<td>' + '<b>MetFrag Score</b>' + '</td>'
                    + '<td>' + '<b>MetFusion Score</b>' + '</td>'
                    + '<td>' + '<b>Fragmenter Score</b>' + '</td>'
                    + '<td>' + '<b>Suspectlist Score</b>' + '</td>'
                    + '<td>' + '<b>Explained Peaks</b>' + '</td>'
                    + '<td>' + '<b>MetFrag Web</b>' + '</td>'
                    + '<td>' + '<b>External Links</b>' + '</td>'
                    + '<td class="tdmax">' + '<b>InChI</b>' + '</td>'
                    + '</tr>\n'))

            # Compound candidate
            if (candidates < max_candidates):
                # Column variables
                inchi = row["InChI"]
                smiles = row["SMILES"]
                mol_formula = row["MolecularFormula"]
                compound_name = row["IUPACName"]
                frag_score = row["FragmenterScore"]
                metfusion_score = row["OfflineMetFusionScore"]
                score = row["Score"]
                if "SuspectListScore" in row:
                    suspectlist_score = row["SuspectListScore"]
                else:
                    suspectlist_score = 0
                peaks_explained = row["NoExplPeaks"]
                peaks_used = row["NumberPeaksUsed"]
                spectrum_explained = row["ExplPeaks"]
                spectrum_explained_formulas = row["FormulasOfExplPeaks"]
                identifier = row["Identifier"]

                # PubChem Synonyms
                if pubchem_synonyms_enabled:
                    pubchem_synonyms = fetch_pubchem_synonyms(inchi)
                else:
                    pubchem_synonyms = '&nbsp;'

                # Compound Classes
                if classyfire_classes_enabled:
                    compound_classes = fetch_classyfire_classes(inchi)
                else:
                    compound_classes = '&nbsp;'

                # Draw Spectrum
                spectrum = re.sub(' .*', '', re.sub('.*PeakListString=', '',
                                                    row["MetFragCLIString"]))
                spectrum_string = plot_spectrum(spectrum, spectrum_explained,
                                                spectrum_explained_formulas)

                # Draw SVG
                svg_string = cdk_inchi_to_svg(str(inchi))

                # External links
                external_links = str(
                    '<a target="_blank" href="' + pubchem_link(
                        compound_name) + '">PubChem</a>' + ', ' +
                    '<a target="_blank" href="' + kegg_link(
                        compound_name) + '">KEGG</a>' + ', ' +
                    '<a target="_blank" href="' + hmdb_link(
                        compound_name) + '">HMDB</a>' + ', ' +
                    '<a target="_blank" href="' + biocyc_link(
                        compound_name) + '">BioCyc</a>' + ', ' +
                    '<a target="_blank" href="' + chebi_link(
                        inchi) + '">ChEBI</a>')

                # MetFragWeb
                FragmentPeakMatchAbsoluteMassDeviation = str(
                    '' +
                    re.sub(' .*', '',
                           re.sub(
                               '.*FragmentPeakMatchAbsoluteMassDeviation=',
                               'FragmentPeakMatchAbsoluteMassDeviation=',
                               row["MetFragCLIString"]))
                )
                FragmentPeakMatchRelativeMassDeviation = str(
                    '' +
                    re.sub(' .*', '',
                           re.sub(
                               '.*FragmentPeakMatchRelativeMassDeviation=',
                               'FragmentPeakMatchRelativeMassDeviation=',
                               row["MetFragCLIString"]))
                )
                DatabaseSearchRelativeMassDeviation = str(
                    '' +
                    re.sub(' .*', '',
                           re.sub(
                               '.*DatabaseSearchRelativeMassDeviation=',
                               'DatabaseSearchRelativeMassDeviation=',
                               row["MetFragCLIString"]))
                )
                IonizedPrecursorMass = str(
                    'IonizedPrecursorMass=' + str(row["precursor_mz"]))
                NeutralPrecursorMass = str(
                    '' + re.sub(' .*', '',
                                re.sub(
                                    '.*NeutralPrecursorMass=',
                                    'NeutralPrecursorMass=',
                                    row[
                                        "MetFragCLIString"]))
                )
                NeutralPrecursorMolecularFormula = str(
                    'NeutralPrecursorMolecularFormula=' + str(
                        row["MolecularFormula"]))
                PrecursorIonMode = str(
                    '' + re.sub(' .*', '', re.sub('.*PrecursorIonMode=',
                                                  'PrecursorIonMode=',
                                                  row["MetFragCLIString"])))
                PeakList = str(
                    '' + re.sub(' .*', '',
                                re.sub('.*PeakListString=', 'PeakList=',
                                       row["MetFragCLIString"])))
                MetFragDatabaseType = str(
                    '' + re.sub(' .*', '',
                                re.sub(
                                    '.*MetFragDatabaseType=',
                                    'MetFragDatabaseType=',
                                    row["MetFragCLIString"])))

                metfrag_web = str(
                    'https://msbi.ipb-halle.de/MetFrag/landing.xhtml?' +
                    FragmentPeakMatchAbsoluteMassDeviation + '&' +
                    FragmentPeakMatchRelativeMassDeviation + '&' +
                    DatabaseSearchRelativeMassDeviation + '&' +
                    IonizedPrecursorMass + '&' +
                    NeutralPrecursorMass + '&' +
                    # NeutralPrecursorMolecularFormula + '&' +
                    PrecursorIonMode + '&' +
                    PeakList + '&' +
                    MetFragDatabaseType)

                # Write html code
                metfrag_html.write(str('<tr style="vertical-align:center">' +
                                       '<td class="tdmax">' +
                                       spectrum_string +
                                       '</td>' +
                                       '<td class="tdmax">' +
                                       svg_string +
                                       '</td>' +
                                       '<td>' +
                                       monoisotopic_mass +
                                       '</td>' +
                                       '<td>' +
                                       mol_formula +
                                       '</td>' +
                                       '<td>' +
                                       compound_name +
                                       '</td>' +
                                       '<td class="tdvar">' +
                                       pubchem_synonyms +
                                       '</td>' +
                                       '<td>' +
                                       compound_classes + '</td>' +
                                       '<td>' +
                                       str(round(float(score), 3)) +
                                       '</td>' +
                                       '<td>' +
                                       str(round(float(metfusion_score), 3)) +
                                       '</td>' +
                                       '<td>' +
                                       str(round(float(frag_score), 3)) +
                                       '</td>' +
                                       '<td>' +
                                       str(
                                           round(float(suspectlist_score), 3)
                                       ) +
                                       '</td>' +
                                       '<td>' +
                                       peaks_explained +
                                       ' / ' +
                                       peaks_used +
                                       '</td>' +
                                       '<td>' +
                                       '<a  target="_blank" href="' +
                                       metfrag_web +
                                       '">MetFragWeb</a>' +
                                       '</td>' +
                                       '<td>' +
                                       external_links +
                                       '</td>' +
                                       '<td class="tdmax">' +
                                       inchi +
                                       '</td>' +
                                       '</tr>\n'))

        line_count += 1
        candidates += 1

    # Finish candidate list
    metfrag_html.write(str('</table>\n'))

# Write html footer
metfrag_html.write('\n</body>\n')
metfrag_html.write('</html>\n')

# Close output html file
metfrag_html.close()
