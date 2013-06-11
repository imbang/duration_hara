__author__ = 'bayu imbang laksono'
# Measure The Duration of High Frequency Energy Radiation
# =======================================================
# based on Measurement of the duration of high-frequency energy radiation
#          and its application to determination of the magnitudes of large shallow earthquake, Tatsuko Hara,2007
#
# slightly modified on :
#       - type window of smoothing process

# next tasks :
# 1. RMEAN
# 2. Error when DELETE
# 3.
# ================================================================================================================
import os
from obspy.core import read,UTCDateTime
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
        self.idxP = -9999
        if self.trace.stats.sac.t7<>self.trace.stats.sac.t0:
            self.duration = self.trace.stats.sac.t0 - self.trace.stats.sac.t7
        else:
            self.duration = -9999
        self.p_pp = -9999

    def __str__(self):
        return "%4s %2s %3s %.3f %5.1f" % (self.station,self.network,self.channel,self.trace.stats.sac.gcarc,self.duration)

    def getDelta(self):
        return self.trace.stats.sac.gcarc

    def basename(self):
        return os.path.basename(self.namafile)
    
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

    def getIdxP(self):
        return self.idxP

    def setIdxP(self,value):
        self.idxP = value

    def getP_PP(self):
        return self.p_pp

    def setP_PP(self,dt):
        self.p_pp = dt


class Duration:
    def __init__(self,directSAC=None,dirOut=None):
        self.dirSAC = directSAC
        self.dirout = dirOut
        self.sac = []
        self.stack = None
        self.stacnt = 0
        self.isLoad = False
        self.isSorted = False
        self.isStacked = False


    def setDirSAC(self,dir):
        self.dirSAC = dir

    def setDirOut(self,dir):
        self.dirout = dir
        
    def loadSAC(self):
        self._readDirSAC(self.dirSAC)
        self.isLoad = True
        self.isSorted = False
        self.isStacked = False
        self.stack = None

    def reloadStage1(self):
        dir = os.path.join(self.dirout,'stage1')
        del self.sac
        self.sac = []
        self._readDirSAC(dir)

    def reloadStage2(self):
        dir = os.path.join(self.dirout,'stage2')
        del self.sac
        self.sac = []
        self._readDirSAC(dir)

    def reloadStage3(self):
        dir = os.path.join(self.dirout,'stage3')
        del self.sac
        self._readDirSAC(dir)

    def _readDirSAC(self,dir):
        listfiles = glob(dir+'/*.SAC')
        for fl in listfiles:
            #print "reading %s" % fl
            print ".",
            tmp = durStation(fl)
            self.sac.append(tmp)
        self.stacnt = len(self.sac)
        print "OK"

    def _readFileSAC(self,dir,idx):
        fl = os.path.join(dir,self.sac[idx].basename())
        tmp = durStation(fl)
        self.sac[idx] = tmp


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
        print i,str(self.sac[i])
        os.system('sac < plot.m')
        os.unlink('plot.m')
        self._readFileSAC(newdir,i)

    def delete(self,idx):
        del self.sac[idx]

    def list(self):
        pass
    
    def rmean(self):
        pass

    def sort(self):
        ok = False
        self.isSorted = True
        while not ok:
            ok = True
            for j in range(len(self.sac)-1):
                tr1 = self.sac[j].getDelta()
                tr2 = self.sac[j+1].getDelta()
                if tr1>tr2:
                    tmp = self.sac[j]
                    self.sac[j] = self.sac[j+1]
                    self.sac[j+1] = tmp
                    ok = False


    def stacking(self):
        if not self.isSorted:
            print "sorting traces..."
            self.sort()

        dtPPP = self.sac[len(self.sac)-1].getP_PP() + 10.0
        delta = self.sac[len(self.sac)-1].trace.stats.delta
        nptsPPP = round(dtPPP / delta)
        t = np.linspace(0,nptsPPP * delta, nptsPPP)
        yy1 = np.zeros(nptsPPP)

        if not self.isStacked:
            self.isStacked = True
            self.stack = self.sac[0].trace.copy()
            self.stack.stats.starttime = ('2000-01-01T0:0:0')

        for i in range(0,len(self.sac)):
            a5 = self.sac[i].trace.stats.sac.a - 5.0
            tdur = a5 + t
            deltax = self.sac[i].trace.stats.delta
            npts = self.sac[i].trace.stats.npts
            x = np.linspace(0,npts * deltax, npts)
            y = self.sac[i].trace.data.copy()
            yy = np.interp(tdur,x,y)
            yy1 = np.add(yy1,yy)
            self.sac[i].trace.data = yy.copy()
            self.sac[i].trace.stats.npts = len(yy)
            self.sac[i].trace.stats.sac.a = 5.0
            self.sac[i].trace.stats.sac.t0 = self.sac[i].trace.stats.sac.t0 - a5
            self.sac[i].trace.stats.sac.t8 = self.sac[i].trace.stats.sac.t8 - a5
            self.sac[i].trace.stats.sac.t7 = self.sac[i].trace.stats.sac.t7 - a5

        yy1 = yy1 / len(self.sac)
        self.stack.data = yy1.copy()
        self.stack.stats.sac.a = 5.0
        self.stack.stats.sac.t0 = -12345
        self.stack.stats.sac.t8 = -12345
        self.stack.stats.sac.t7 = -12345
        self.stack.stats.network = "XX"
        self.stack.stats.station = "XXXX"

        idx = np.argmax(self.stack.data)
        self.stack.stats.sac.t8 = float(self.stack.stats.delta * idx)

        valmax = self.stack.data[idx]
        val25 = 0.75 * valmax
        idx25 = (np.abs(self.stack.data[idx:] - val25)).argmin()
        idx25 = idx + idx25
        tmp = float(self.stack.stats.delta * idx25)
        self.stack.stats.sac.t7 = tmp

        self.stack.stats.update()
        newdir = os.path.join(self.dirout,'stack')
        if not os.path.isdir(newdir):
            os.makedirs(newdir)
        nmfile = os.path.join(newdir,'stack.SAC')
        self.stack.write(nmfile,format='SAC')
        print "save to dir stack --------"

    def adjustP(self,idx):
        try:
            id = int(idx)
            self._adjustP(id)
        except:
            for i in range(0,len(self.sac)):
                self._adjustP(i)


    def _adjustP(self,i):
        self.sac[i].trace.stats.sac.t0 = float(self.sac[i].trace.stats.sac.a)  + self.sac[i].getP_PP()

    def clear(self):
        del self.sac
        self.sac=[]

    def stage1(self):
        self.bandpass()
        self.putP_PP('all')
        self.saveToDir('stage1')

    def stage2(self):
        self.power2()
        self.normalize('all')
        self.saveToDir('stage2')

    def stage3(self,smt):
        self.smooth(smt)
        self.peaks()
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
            ta = self.sac[i].getIdxP()
            #if i==0: print ta
            tpp = self.sac[i].getIdxPP()
            #if i==0: print tpp

            t = np.linspace(0, tr.stats.delta * tr.stats.npts, tr.stats.npts)
            tx = t[(t>=ta)&(t<=tpp)]
            #if i==0: print tr.stats.sac.a,tr.stats.sac.t0,tx[0],tx[len(tx)-1]
            idxa = np.where(t==tx[0])[0][0]
            idxpp = np.where(t==tx[len(tx)-1])[0][0]
            #if i==0: print "idxa,idxpp,",idxa,idxpp
            #idx = np.argmax(tr.data[idxa:idxpp])
            idx = np.argmax(tr.data[idxa:idxpp])
            #if i==0: print idx
            self.sac[i].setIdxMax(idxa+idx)
            tr.stats.sac.t8 = float(tr.stats.delta * (idxa+idx))

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
#            dur = float(self.sac[i].trace.stats.sac.t0) - tmp
            dur = tmp - float(self.sac[i].trace.stats.sac.a)
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
            nmfile = self.sac[i].basename()
            nmfile = os.path.join(newdir,nmfile)
            tr = self.sac[i].trace
            tr.write(nmfile,format='SAC')


    def _getTravelTimes(self,i):
        hipo = float(self.sac[i].trace.stats.sac.evdp) / 1000
        delta = float(self.sac[i].trace.stats.sac.gcarc)
        fid = open('data.inp','w')
        fid.write('n\n')
        nmtable='/home/sysop/iaspei-tau/tables/iasp91'
        if not os.path.exists(nmtable):
            print "nama file model kecepatan tidak ditemukan..."
            return
        fid.write(nmtable+'\n')
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
                idx1 = int (tmp1 / self.sac[i].trace.stats.delta)
                #if i==0: print "P",idx1
                #self.sac[i].setIdxP(idx1)
                self.sac[i].setIdxP(float(self.sac[i].trace.stats.sac.o) + tmp1)
                self.sac[i].trace.stats.sac.a = float(self.sac[i].trace.stats.sac.o) + tmp1
            elif tmp[1]=="PP":
                tmp2 = float(tmp[2])
                idx2 = int (tmp2 / self.sac[i].trace.stats.delta)
                #if i==0: print "PP",idx2
                #self.sac[i].setIdxPP(idx2)
                self.sac[i].setIdxPP(float(self.sac[i].trace.stats.sac.o) + tmp2)
                self.sac[i].trace.stats.sac.t0 = float(self.sac[i].trace.stats.sac.o) + tmp2
                self.sac[i].setP_PP(tmp2 - tmp1)
                break

    def report(self):
        print "============================================================="
        for i in range(0,len(self.sac)):
            print i,str(self.sac[i])
        print "============================================================="    


if __name__ == "__main__":

    hd = durProcessing()

    hd.report()

