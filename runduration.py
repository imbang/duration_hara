import os

# Measure The Duration of High Frequency Energy Radiation
# =======================================================
# based on Measurement of the duration of high-frequency energy radiation
#          and its application to determination of the magnitudes of large shallow earthquake, Tatsuko Hara,2007
#
# slightly modified on :
#       - type window of smoothing process
# ================================================================================================================


__author__ = 'bayu imbang laksono'
from obspy.core import read
from glob import glob
import numpy as np

from common import smooth


class durStation:
    def __init__(self,nmfile):
        self.namafile = nmfile
        st = read(nmfile)
        self.trace = st[0]
        self.t10 = None
        self.t25 = None
        self.t33 = None
        self.t50 = None
        self.t75 = None
        self.t80 = None
        self.t100 = None
        self.idxmax = -9999

    def namafile(self):
        return self.namafile
    
    def trace(self):
        return self.trace

    def getDuration(self):
        return "%5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f" % (self.t10, self.t25, self.t33, self.t50,\
        self.t75, self.t80, self.t100)
    
    def setDuration(self):
        pass

    def getIdxMax(self,value):
        return self.idxmax

    def setIdxMax(self,value):
        self.idxmax = value

class durProcessing:
    def __init__(self,directSAC='/home/sysop/TITIPANSUGENG/1_bku2000/data',dirOut='/home/sysop/TITIPANSUGENG'):
        self.dirSAC = directSAC
        self.dirout = dirOut
        self.sac = []
        self.readDirSAC()

    def readDirSAC(self):
        listfiles = glob(self.dirSAC+'/*.SAC')
        for fl in listfiles:
            #print "reading %s" % fl
            print ".",
            tmp = durStation(fl)
            self.sac.append(tmp)
        print "OK"

    def infoSAC(self,value):
        for i in range(0,len(self.sac)):
            if i==value:
                return self.sac[i].trace

    def plot(self,value):
        for i in range(0,len(self.sac)):
            if i==value:
                return self.sac[i].trace

    def stage1(self):
        pass

    def bandpass(self,fn=2,fx=4):
        # filter bandpass
        for i in range(0,len(self.sac)):
            self.sac[i].trace.filter('bandpass',freqmin=fn,freqmax=fx,corners=2,zerophase=True)

    def power2(self):
        # square trace
        for i in range(0,len(self.sac)):
            tr = self.sac[i].trace
            tmp = np.power(tr.data,2)
            tr.data = tmp.copy()

    def putP_PP(self):
        # detect P wave
        self.getTravelTimes()

    def peaks(self):
        # find peak
        for i in range(0,len(self.sac)):
            tr = self.sac[i].trace
            idx = np.argmax(tr.data)
            self.sac[i].setIdxMax(idx)
            tr.stats.sac.t8 = float(tr.stats.delta * idx)

    def smooth(self,value):
        # smoothing
        for i in range(0,len(self.sac)):
            tr = self.sac[i].trace
            data = tr.data.copy()
            datasmooth = smooth(data,value,'hanning')
            tr.data = datasmooth.copy()

    def duration(self):
        for i in range(0,len(self.sac)):
            idx = self.sac[i].getIdxMax()
            valmax = self.sac[i].data[idx]
            val25 = 75.0 / 100.0 * valmax

    def stage8(self):
        pass

    def save(self,newdir):
        newdir = os.path.join(self.dirout,newdir)
        if not os.path.isdir(newdir):
            os.makedirs(newdir)

        for i in range(0,len(self.sac)):
            nmfile = os.path.basename(self.sac[i].namafile)
            #nmfile = nmfile +'.1'
            nmfile = os.path.join(newdir,nmfile)
            tr = self.sac[i].trace
            tr.write(nmfile,format='SAC')

    def getTravelTimes(self):
        for i in range(0,len(self.sac)):
            hipo = float(self.sac[i].trace.stats.sac.evdp) / 1000
            delta = float(self.sac[i].trace.stats.sac.gcarc)
            fid = open('data.inp','w')
            fid.write('n\n')
            fid.write('/home/sysop/iaspei-tau/tables/iasp91\n')
            fid.write('P\n')
            fid.write('PP\n')
            fid.write('\n')
            fid.write('n\n')
            fid.write('%f\n' % hipo)
            fid.write('%f\n' % delta)
            fid.write('-5\n')
            fid.write('-5\n')
            fid.close()
            os.system('ttimes < data.inp > /dev/null 2>&1')

            fid = open('ttimes.lst','r')
            fid.readline()
            fid.readline()
            fid.readline()
            fid.readline()
            fid.readline()
            while 1:
                tmp = fid.readline()
                tmp = tmp.split()
                if len(tmp)==0:
                    break
                if tmp[1]=="P":
                    self.sac[i].trace.stats.sac.a = float(self.sac[i].trace.stats.sac.o) + float(tmp[2])
                elif tmp[1]=="PP":
                    self.sac[i].trace.stats.sac.t0 = float(self.sac[i].trace.stats.sac.o) + float(tmp[2])
            


if __name__ == "__main__":

    hd = durProcessing()
    hd.bandpass()
    hd.putP_PP()
    hd.save('coba')
    hd.power2()
    hd.save('coba1')
    hd.stage6(400)
    hd.peaks()
    hd.save('coba2')

