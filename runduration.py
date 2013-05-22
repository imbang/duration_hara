import os

__author__ = 'bayu imbang laksono'
from obspy.core import read
from glob import glob
import numpy as np

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    return y


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

    def namafile(self):
        return self.namafile
    
    def trace(self):
        return self.trace

    def getDuration(self):
        return "%5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f" % (self.t10, self.t25, self.t33, self.t50,\
        self.t75, self.t80, self.t100)

class durProcessing:
    def __init__(self,directSAC='/home/sysop/TITIPANSUGENG/1_bku2000/data'):
        self.dirSAC = directSAC
        self.sac = []
        self.readDirSAC()

    def readDirSAC(self):
        listfiles = glob(self.dirSAC+'/*.SAC')
        for fl in listfiles:
            print "reading %s" % fl
            tmp = durStation(fl)
            self.sac.append(tmp)

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

    def stage5(self):
        # find peak
        pass

    def stage6(self):
        pass

    def stage7(self):
        pass

    def stage8(self):
        pass

    def saveto(self,newdir):
        newdir = os.path.join(self.dirSAC,newdir)
        if not os.path.isdir(newdir):
            os.makedirs(newdir)

        for i in range(0,len(self.sac)):
            nmfile = os.path.basename(self.sac[i].namafile)
            nmfile = nmfile +'.1'
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
            os.system('ttimes < data.inp')

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
                    #print tmp[2]
                    self.sac[i].trace.stats.sac.a = float(self.sac[i].trace.stats.sac.o) + float(tmp[2])
                elif tmp[1]=="PP":
                    #print tmp[2]
                    self.sac[i].trace.stats.sac.t0 = float(self.sac[i].trace.stats.sac.o) + float(tmp[2])
            


if __name__ == "__main__":
    #hd = durProcessing('')
    hd = durProcessing('/home/sysop/TITIPANSUGENG/1_bku2000/data')
    hd.stage2()

    tmp = hd.infoSAC(0)
    tmp.plot()
