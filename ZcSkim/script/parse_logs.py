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


lumi_e  = 19789.0
lumi_mu = 19751.0


def main():

    if len(sys.argv) < 2:
        print "usage: ./parse_logs.py <channel>"
        print "   channel = DYJetsToLL, Ztautau, TTbar, WW, WZ, ZZ, Wj"  
        return

    channel = sys.argv[1]

    pattern = "Weighted yield (tagged)    ="
    #pattern = "Weighted yield (inclusive) ="

    is_tautau = 0
    if ( channel=="Ztautau"):
        is_tautau = 1
        #pattern = "Z-->tautau weighted yield (inclusive) ="
        pattern = "Z-->tautau weighted yield (tagged)    ="
        channel = "DYJetsToLL"
        
    scale = w[channel]

    p = subprocess.Popen(['ls', path+channel],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = p.communicate()
    #print output

    n_ele,err_ele,n_muo,err_muo = 0.,0.,0.,0.

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
                        n_ele   += float(line.split()[4])
                        err_ele += float(line.split()[5])
                    else:
                        n_ele   += float(line.split()[5])
                        err_ele += float(line.split()[6])

                if counter==2:
                    print counter, line.split()
                    if not is_tautau:
                        n_muo   += float(line.split()[6])
                        err_muo += float(line.split()[7])
                    else:
                        n_muo   += float(line.split()[7])
                        err_muo += float(line.split()[8])
                        

        f.close()

    print "\n"
    print "electron channel: %f +- %f" % (n_ele*scale*lumi_e,math.sqrt(err_ele)*scale*lumi_e)
    print "muon channel:     %f +- %f" % (n_muo*scale*lumi_mu,math.sqrt(err_muo)*scale*lumi_mu)
        

if __name__ == "__main__":
    main()
