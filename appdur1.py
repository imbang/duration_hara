import numpy as np
import duration as dr
import matplotlib as mpl
mpl.use('WXAgg')
from matplotlib.figure import Figure
import wx
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx

class CanvasPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()

    def draw(self,tr):
        t = np.linspace(0, tr.stats.delta * tr.stats.npts, tr.stats.npts)
        s = tr.data
        self.axes.plot(t, s)

if __name__ == "__main__":
    hd = dr.Duration()
    hd.setDirSAC('/home/sysop/TITIPANSUGENG/1_bku2000/data')
    hd.setDirOut('/home/sysop/dur')
    hd.loadSAC()
    hd.bandpass()
    hd.sort()
    hd.report()
    tmp = hd.getStation(10)

    app = wx.PySimpleApp()
    fr = wx.Frame(None, title='test')
    panel = CanvasPanel(fr)
    panel.draw(tmp.trace)
    fr.Show()
    app.MainLoop()

#hd.bandpass()
#hd.putP_PP('all')
#hd.saveToDir('coba')
#hd.power2()
#hd.normalize('all')
#hd.peaks()
#hd.saveToDir('coba1')
#hd.smooth(200)
#hd.duration()
#hd.saveToDir('coba2')
