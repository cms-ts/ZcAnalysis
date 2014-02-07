#!/usr/bin/python

import sys
import subprocess
import math


#channel = "DYJetsToLL"
#channel = "TTbar"
#channel = 'WW'

path = "/gpfs/cms/users/casarsa/analysis/Zc/work/output/v01/"

w = {'Wj'        : 31200./57709905.,
     'WW'        : 54.838/10000431.,
     'WZ'        : 33.21/10000283. ,
     'ZZ'        : 8.059/9799908.  ,
     'TTbar'     : 225.197/6923750.,
     'DYJetsToLL': 3503.71/30459503.}


lumi_ee = 19789.0
lumi_mm = 19751.0
lumi_em = 19780.0

def main():

    if len(sys.argv) < 2:
        print "usage: ./parse_logs.py <channel>"
        print "   channel = DYJetsToLL, Ztautau, TTbar, WW, WZ, ZZ, Wj"  
        return

    channel = sys.argv[1]

    #pattern = "Weighted yield (inclusive) ="
    pattern = "Weighted yield (tagged)    ="

    is_tautau = 0
    if ( channel=="Ztautau"):
        is_tautau = 1
        #pattern = "Z-->tautau weighted yield (inclusive) ="
        pattern = "Z-->tautau weighted yield (tagged)    ="
        channel = "DYJetsToLL"
        
    scale = w[channel]

    p = subprocess.Popen(['ls', path+channel],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = p.communicate()

    n_ee, err_ee, n_mm, err_mm, n_em, err_em = 0., 0., 0., 0., 0., 0.

    for idir in output.split():
        print idir
        f = open(path+channel+'/'+idir+'/job.log', 'r')

        counter = 0

        for line in f.readlines():

            if line.find(pattern) > -1:
                counter += 1
                if counter==1:
                    print counter, line.split()
                    if not is_tautau:
                        n_ee   += float(line.split()[4])
                        err_ee += float(line.split()[5])
                    else:
                        n_ee   += float(line.split()[4])
                        err_ee += float(line.split()[5])

                if counter==2:
                    print counter, line.split()
                    if not is_tautau:
                        n_mm   += float(line.split()[4])
                        err_mm += float(line.split()[5])
                    else:
                        n_mm   += float(line.split()[4])
                        err_mm += float(line.split()[5])

                if counter==3:
                    print counter, line.split()
                    if not is_tautau:
                        n_em   += float(line.split()[4])
                        err_em += float(line.split()[5])
                    else:
                        n_em   += float(line.split()[4])
                        err_em += float(line.split()[5])
                        

        f.close()

    #print  n_ee,err_ee,n_mm,err_mm
    print "\n"
    print "dielectron channel:    %f +- %f" % (n_ee*scale*lumi_ee,math.sqrt(err_ee)*scale*lumi_ee)
    print "dimuon channel:        %f +- %f" % (n_mm*scale*lumi_mm,math.sqrt(err_mm)*scale*lumi_mm)
    print "electron-muon channel: %f +- %f" % (n_em*scale*lumi_em,math.sqrt(err_em)*scale*lumi_em)
        

if __name__ == "__main__":
    main()
