__author__ = 'bayu imbang laksono'
# Measure The Duration of High Frequency Energy Radiation
# =======================================================
# based on Measurement of the duration of high-frequency energy radiation
#          and its application to determination of the magnitudes of large shallow earthquake, Tatsuko Hara,2007
#
# slightly modified on :
#       - type window of smoothing process
# ================================================================================================================
import os
from obspy.core import read
from glob import glob
import numpy as np

from common import smooth


class durStation:
    def __init__(self,nmfile):
        self.namafile = nmfile
        st = read(nmfile)
        self.trace = st[0]
        self.channel = self.trace.stats.channel
        self.network = self.trace.stats.network
        self.station = self.trace.stats.station
        self.idxmax = -9999
        self.idxPP = -9999
        self.duration = -9999
        self.p_pp = -9999

    def __str__(self):
        return "%4s %2s %3s %5.1f" % (self.station,self.network,self.channel,self.duration)

    def namafile(self):
        return self.namafile
    
    def trace(self):
        return self.trace

    def setTrace(self,trace):
        self.trace = trace

    def getDuration(self):
        return self.duration
    
    def setDuration(self,value):
        self.duration = value

    def getIdxMax(self):
        return self.idxmax

    def setIdxMax(self,value):
        self.idxmax = value

    def getIdxPP(self):
        return self.idxPP

    def setIdxPP(self,value):
        self.idxPP = value

    def getP_PP(self):
        return self.p_pp

    def setP_PP(self,dt):
        self.p_pp = dt


class Duration:
    def __init__(self,directSAC='/home/sysop/TITIPANSUGENG/1_bku2000/data',dirOut='/home/sysop/TITIPANSUGENG'):
        self.dirSAC = directSAC
        self.dirout = dirOut
        self.sac = []
        self._readDirSAC()
        self.stacnt = 0

    def _readDirSAC(self):
        listfiles = glob(self.dirSAC+'/*.SAC')
        for fl in listfiles:
            #print "reading %s" % fl
            print ".",
            tmp = durStation(fl)
            self.sac.append(tmp)
        self.stacnt = len(self.sac)
        print "OK"

    def getStation(self,value):
        return self.sac[value]
    
    def plot(self,idx):
        try:
            id = int(idx)
            self._plot(id)
        except:
            for i in range(0,len(self.sac)):
                self._plot(i)

    def _plot(self,i):
        self.saveToDir('tmp')
        newdir = os.path.join(self.dirout,'tmp')
        nmfile = os.path.basename(self.sac[i].namafile)
        nmfile = os.path.join(newdir,nmfile)
        fid = open('plot.m','w')
        fid.write('r %s\n' % nmfile)
        fid.write('qdp off\n')
        fid.write('ppk\n')
        fid.write('wh\n')
        fid.write('q\n')
        fid.close()
        os.system('sac < plot.m')
        os.unlink('plot.m')

    def delete(self,idx):
        pass

    def list(self):
        pass
    
    def rmean(self):
        pass

    def adjustP(self,idx):
        try:
            id = int(idx)
            self._adjustP(id)
        except:
            for i in range(0,len(self.sac)):
                self._adjustP(i)


    def _adjustP(self,i):
        self.sac[i].trace.stats.sac.t0 = float(self.sac[i].trace.stats.sac.a)  + self.sac[i].getP_PP()

    def reload(self):
        self._readDirSAC()

    def clear(self):
        pass

    def stage1(self):
        self.bandpass()
        self.putP_PP()
        self.saveToDir('stage1')

    def stage2(self):
        self.power2()
        self.normalize('all')
        self.peaks()
        self.saveToDir('stage2')

    def stage3(self,smt):
        self.smooth(smt)
        self.duration()
        self.saveToDir('stage3')

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

    def putP_PP(self,idx):
        try:
            id = int(idx)
            self._getTravelTimes(id)
        except:
            for i in range(0,len(self.sac)):
                self._getTravelTimes(i)

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
            valmax = self.sac[i].trace.data[idx]
            val25 = 0.75 * valmax
            idx25 = (np.abs(self.sac[i].trace.data[idx:] - val25)).argmin()
            idx25 = idx + idx25
            tmp = float(self.sac[i].trace.stats.delta * idx25)
            self.sac[i].trace.stats.sac.t7 = tmp
            dur = float(self.sac[i].trace.stats.sac.t0) - tmp
            self.sac[i].setDuration(dur)

    def normalize(self,idx):
        try:
            id = int(idx)
            self._normalize(id)
        except:
            for i in range(0,len(self.sac)):
                self._normalize(i)

    def _normalize(self,i):
        tr = self.sac[i].trace
        tr.normalize()
        self.sac[i].setTrace(tr)


    def saveToDir(self,newdir):
        newdir = os.path.join(self.dirout,newdir)
        if not os.path.isdir(newdir):
            os.makedirs(newdir)

        for i in range(0,len(self.sac)):
            nmfile = os.path.basename(self.sac[i].namafile)
            nmfile = os.path.join(newdir,nmfile)
            tr = self.sac[i].trace
            tr.write(nmfile,format='SAC')


    def _getTravelTimes(self,i):
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
                tmp1 = float(tmp[2])
                self.sac[i].trace.stats.sac.a = float(self.sac[i].trace.stats.sac.o) + tmp1
            elif tmp[1]=="PP":
                tmp2 = float(tmp[2])
                self.sac[i].trace.stats.sac.t0 = float(self.sac[i].trace.stats.sac.o) + tmp2
                self.sac[i].setP_PP(tmp2 - tmp1)

    def report(self):
        print "============================================================="
        for i in range(0,len(self.sac)):
            print i,str(self.sac[i])
        print "============================================================="    


if __name__ == "__main__":

    hd = durProcessing()
    hd.bandpass()
    hd.putP_PP('all')
    hd.saveToDir('coba')
    hd.power2()
    hd.normalize('all')
    hd.peaks()
    hd.saveToDir('coba1')
    hd.smooth(200)
    hd.duration()
    hd.saveToDir('coba2')
    hd.report()

