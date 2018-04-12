import argparse
import csv
import os

parser = argparse.ArgumentParser()
parser.add_argument('--input')
parser.add_argument('--db_local')
parser.add_argument('--db_online')
parser.add_argument('--ppm')
parser.add_argument('--ppm_frag')
parser.add_argument('--fragmasstol')
parser.add_argument('--polarity')
parser.add_argument('--results')

args = parser.parse_args()
print args

os.makedirs("tmet")

with open(args.input,"r") as infile:
    numlines = 0
    for line in infile:
        line = line.strip()
        if numlines == 0:
            if "NAME" in line:
                featid = line.split("NAME: ")[1]
            if "PRECURSORMZ" in line:
                mz = float(line.split("PRECURSORMZ: ")[1])
                if args.polarity=="pos":
                    mz2 = mz-1.007276
                else:
                    mz2 = mz+1.007276
            if "Num Peaks" in line:
                numlines = int(line.split("Num Peaks: ")[1])
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
                cmd_command = "PeakListPath=tmpspec.txt "
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
                cmd_command += "NumberThreads=1 "
                # run Metfrag
                print "metfrag {0}".format(cmd_command)
                os.system("metfrag {0}".format(cmd_command))
            else:
                line = tuple(line.split("\t"))
                linesread += 1
                peaklist.append(line)


#merge outputs
outfiles = os.listdir("./tmet")
with open(args.results, 'a') as outfile:
    first_write = True
    for fname in outfiles:
        fileid = os.path.basename(fname)
        fileid = fileid.split("_")[0]
        with open("./tmet/"+fname) as infile:
            reader = csv.DictReader(infile, delimiter=',', quotechar='"')
            for line in reader:
                if first_write:
                    headers = ['UID']
                    headers.extend(line.keys())
                    dwriter = csv.DictWriter(outfile, fieldnames=headers, delimiter='\t')
                    dwriter.writerow(dict((fn,fn) for fn in headers))
                    first_write = False 
                line['UID'] = fileid
                dwriter.writerow(line)


