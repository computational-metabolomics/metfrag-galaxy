import argparse
import csv
import os
import sys
from collections import defaultdict

print(sys.version)

parser = argparse.ArgumentParser()
parser.add_argument('--input')
parser.add_argument('--db_local')
parser.add_argument('--db_online')
parser.add_argument('--ppm')
parser.add_argument('--ppm_frag')
parser.add_argument('--fragmasstol')
parser.add_argument('--polarity')
parser.add_argument('--minMSMSpeaks')
parser.add_argument('--results')
parser.add_argument('--threads')
parser.add_argument('--unconnectcompnd', action='store_true')
parser.add_argument('--isotopefilter', action='store_true')
parser.add_argument('--minelem')
parser.add_argument('--maxelem')
parser.add_argument('--subinclusion')
parser.add_argument('--subexclusion')
parser.add_argument('--eleminclusion')
parser.add_argument('--elemexcluinclusion')
parser.add_argument('--elemexclusion')
parser.add_argument('--score_thrshld')
parser.add_argument('--pctexplpeak_thrshld')

args = parser.parse_args()
print args



#Stock filters in dict (names have been choose to add them directly in MetFragPreProcessingCandidateFilter parameter)
dct_filter = defaultdict(list)
dct_filter['UnconnectedCompoundFilter'].append(args.unconnectcompnd)
dct_filter['IsotopeFilter'].append(args.isotopefilter)
dct_filter['MinimumElementsFilter'].append(args.minelem)
dct_filter['MinimumElementsFilter'].append("FilterMinimumElements")
dct_filter['MaximumElementsFilter'].append(args.maxelem)
dct_filter['MaximumElementsFilter'].append("FilterMaximumElements")
dct_filter['SmartsSubstructureInclusionFilter'].append(args.subinclusion)
dct_filter['SmartsSubstructureInclusionFilter'].append("FilterSmartsInclusionList")
dct_filter['SmartsSubstructureExclusionFilter'].append(args.subexclusion)
dct_filter['SmartsSubstructureExclusionFilter'].append("FilterSmartsExclusionList")
dct_filter['ElementInclusionFilter'].append(args.eleminclusion)
dct_filter['ElementInclusionFilter'].append("FilterIncludedElements")
dct_filter['ElementInclusionExclusiveFilter'].append(args.elemexcluinclusion)
dct_filter['ElementInclusionExclusiveFilter'].append("FilterIncludedElements")
dct_filter['ElementExclusionFilter'].append(args.elemexclusion)
dct_filter['ElementExclusionFilter'].append("FilterExcludedElements")

print "\nAll filters :\n"
for key,val in dct_filter.items():
    print key, "=>", val
#Keep only used filters
dct_keepf = defaultdict(list)
for n in dct_filter:
    if dct_filter[n][0] != '':
        dct_keepf[n] = dct_filter[n]
print "\nFilter(s) to use :\n"
for key,val in dct_keepf.items():
    print key, "=>", val
print "\n"

#Create repo to work
if not os.path.exists("./tmet"):
    os.makedirs("./tmet")
else:
    #Clean the repo if already exist
    for filename in os.listdir("./tmet"):
        os.remove("./tmet" + "/" + filename)

with open(args.input,"r") as infile:
    numlines = 0
    for line in infile:
        line = line.strip()
        if numlines == 0:
            if "CH$NAME:" in line:
                featid = line.split("CH$NAME: ")[1]
                featid = featid.split(" ")
                featid = featid[5].split(":")[1] + "-" + featid[7].split(":")[1] + "-" + featid[3].split(":")[1]
            if "MS$FOCUSED_ION:" in line:
                mz = float(line.split("MS$FOCUSED_ION: PRECURSOR_M/Z ")[1])
                if args.polarity=="pos":
                    mz2 = mz-1.007276
                else:
                    mz2 = mz+1.007276
            if "PK$NUM_PEAK:" in line:
                numlines = int(line.split("PK$NUM_PEAK: ")[1])
                linesread = 0
                peaklist = []
        else:
            if linesread == numlines:
                numlines = 0
                #write spec file
                with open('./tmpspec.txt', 'w') as outfile:
                    for p in peaklist:
                        outfile.write(p[0]+"\t"+p[1]+"\n")
                #create commandline input
                cmd_command = "PeakListPath=./tmpspec.txt "
                if args.db_local != "None":
                    cmd_command += "MetFragDatabaseType=LocalCSV "
                    cmd_command += "LocalDatabasePath={0} ".format(args.db_local)
                else:
                    cmd_command += "MetFragDatabaseType={0} ".format(args.db_online)
                cmd_command += "FragmentPeakMatchAbsoluteMassDeviation={0} ".format(args.fragmasstol)
                cmd_command += "FragmentPeakMatchRelativeMassDeviation={0} ".format(args.ppm_frag)
                cmd_command += "DatabaseSearchRelativeMassDeviation={0} ".format(args.ppm)
                cmd_command += "NeutralPrecursorMass={0} ".format(mz2)
                cmd_command += "SampleName={0}_metfrag ".format(featid)
                cmd_command += "ResultsPath=./tmet/ "
                if args.polarity == "pos":
                    cmd_command += "IsPositiveIonMode=True "
                else:
                    cmd_command += "IsPositiveIonMode=False "
                if args.polarity == "pos": ### Annotation information. Create a dict for the PrecurorIonModes??
                    cmd_command += "PrecursorIonMode=1 "
                else:
                    cmd_command += "PrecursorIonMode=-1 "
                cmd_command += "MetFragCandidateWriter=CSV " ## TSV not available
                cmd_command += "NumberThreads={} ".format(args.threads)

                #All pre processing filters
                param = []
                defparam = defaultdict(list)
                for f in dct_keepf:
                    if f == "UnconnectedCompoundFilter":
                        if dct_keepf[f][0] == True:
                            param .append("UnconnectedCompoundFilter")
                    elif f == "IsotopeFilter":
                        if dct_keepf[f][0] == True:
                            param.append("IsotopeFilter")
                    else:
                        param.append(f)
                        defparam[dct_keepf[f][1]].append(dct_keepf[f][0])

                if param:
                    cmd_command += "MetFragPreProcessingCandidateFilter={} ".format(','.join(param))
                    for key in defparam:
                        cmd_command += "{0}={1} ".format(str(key),str(defparam[key][0]))
                

                #Filter before process with a minimum number of MS/MS peaks
                if linesread >= float(args.minMSMSpeaks):
                    # run Metfrag
                    print "metfrag {0}".format(cmd_command)
                    os.system("metfrag {0}".format(cmd_command))
            else:
                #One line not need between numpeak and peaklist
                if not "PK$PEAK:" in line:
                    line = tuple(line.split("\t"))
                    #Keep only m/z and intensity, not relative intensity
                    save_line = tuple(line[0].split("\t") + line[1].split("\t"))
                    linesread += 1
                    peaklist.append(save_line)


#outputs might have different headers. Need to get a list of all the headers before we start merging the files
outfiles = sorted(os.listdir("./tmet"))

headers = []
c = 0
for fname in outfiles:
    with open("./tmet/"+fname) as infile:
        reader = csv.reader(infile)
        headers.extend(reader.next())
        # check if file has any data rows 
        for i, row in enumerate(reader):
            c+=1
            if i==1:
                break

# if no data rows (e.g. matches) then do not save an output and leave the program        
if c==0:
    sys.exit()

headers = ['UID'] + sorted(list(set(headers)))

#merge outputs
with open(args.results, 'a') as merged_outfile:

    dwriter = csv.DictWriter(merged_outfile, fieldnames=headers, delimiter='\t')
    dwriter.writeheader()

    for fname in outfiles:
        fileid = os.path.basename(fname).split("_")[0]
        with open("./tmet/"+fname) as infile:
            reader = csv.DictReader(infile, delimiter=',', quotechar='"')
            for line in reader:
                bewrite = True
                for key, value in line.items():
                    #Filter when no MS/MS peak matched
                    if key == "ExplPeaks":
                        if "NA" in value:
                            bewrite = False
                    #Filter with a score threshold
                    elif key == "Score":
                        if value <= args.score_thrshld:
                            bewrite = False
                    elif key == "NoExplPeaks":
                        nbfindpeak = float(value)
                    elif key == "NumberPeaksUsed":
                        totpeaks = float(value)
                #Filter with a relative number of peak matched
                pctexplpeak = nbfindpeak / totpeaks * 100
                if pctexplpeak < float(args.pctexplpeak_thrshld):
                    bewrite = False
                #Write the line if it pass all filters
                if bewrite:
                    line['UID'] = fileid
                    dwriter.writerow(line)


