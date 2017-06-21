from pylab import *
from shapely.geometry import *
from read_mist_models import *

class HessDiagram:
    def __init__(self,xmin,ymin,xmax,ymax,nx,ny):
        self.XMIN=xmin; self.XMAX=xmax; self.YMIN=ymin; self.YMAX=ymax
        self.NX=nx; self.NY=ny
        self.Xgrid=linspace(self.XMIN,self.XMAX,self.NX)
        self.Ygrid=linspace(self.YMIN,self.YMAX,self.NY)
        self.boxes=[]
        for i in range(self.NX-1):
            x0,x1=self.Xgrid[i],self.Xgrid[i+1]
            for j in range(self.NY-1):
                y0,y1=self.Ygrid[j],self.Ygrid[j+1]
                self.boxes.append(HessBox(x0,y0,x1,y1))

        self.iso_string=None; self.isoX=None; self.isoY=None
        
    def add_iso(self,iso,DM=0):
        self.iso = iso
        self.isoX = iso['Bessell_V'] - iso['Bessell_I']
        self.isoY = iso['Bessell_V'] + DM
        self.iso_string=LineString(zip(self.isoX,self.isoY))

    def do_mass_intervals(self):
        for box in self.boxes:
            box.iso_intersection(self.iso_string)
            box.set_mass_intervals(self.iso,self.isoX,self.isoY)

    def do_IMF(self,alpha=2.35):
        for box in self.boxes:
            box.IMF_calc(alpha)
            
class HessBox:
    def __init__(self,xmin,ymin,xmax,ymax):
        self.xmin=xmin
        self.ymin=ymin
        self.xmax=xmax
        self.ymax=ymax
        self.center = 0.5*array([self.xmin+self.xmax, self.ymin+self.ymax])
        self.coords=array([xmin,ymin,xmax,ymax])
        self.mass_intervals=array([])
        self.IMF_weight=0.0
        self.intersects=False
        self.intersection=empty((1,2))
        self.lines = self._LS()

    def __call__(self,isoX,isoY):
        return self.get_mass_intervals(isoX,isoY)
        
    def line_plot(self):
        plot([self.xmin,self.xmax],[self.ymin,self.ymin],color='Red')
        plot([self.xmin,self.xmax],[self.ymax,self.ymax],color='Red')
        plot([self.xmin,self.xmin],[self.ymin,self.ymax],color='Red')
        plot([self.xmax,self.xmax],[self.ymin,self.ymax],color='Red')
        
    def _LS(self):
        p0=(self.xmin,self.ymin)
        p1=(self.xmax,self.ymin)
        p2=(self.xmax,self.ymax)
        p3=(self.xmin,self.ymax)
        p4=(self.xmin,self.ymin)
        return LineString((p0,p1,p2,p3,p4))

    def iso_intersection(self,iso_string):
        if iso_string.intersects(self.lines):
            self.intersects = True
            self.intersection = array(iso_string.intersection(self.lines))
        else:
            self.intersects = False
    
    def set_mass_intervals(self,iso,isoX,isoY):
        def between(x,x1,x2):
            test1 = x1 >= x > x2
            test2 = x1 <= x < x2
            return test1 or test2            
        
        def find_point_in_iso(point,isoX,isoY):
            if isscalar(point): return -1, 0.0
            x=point[0]
            y=point[1]
            for i in range(len(isoX)-1):
                if between(x,isoX[i], isoX[i+1]) and between(y,isoY[i],isoY[i+1]):
                    alfa=(isoY[i+1]-y)/(isoY[i+1]-isoY[i])
                    return i, alfa
            return -1, 0.0

        def check(m): return len(masses)>0 and len(masses)%2!=0            
        
        guess=0
        masses=[]
        for point in self.intersection:
            i,alfa=find_point_in_iso(point,isoX,isoY)
            if i >= 0: 
                mass=alfa*iso['initial_mass'][i] + (1-alfa)*iso['initial_mass'][i+1]
                masses.append(mass)
                guess=i

        if check(masses):
            self.mass_intervals=empty((1,2))
        else:
            self.mass_intervals=reshape(masses,(-1,2))

    def IMF_calc(self,alpha=2.35):
        def N(M,alpha): return pow(M,1-alpha)/(1-alpha)
        IMF_weight = 0.0
        for point in self.mass_intervals:
            Mlo=min(point)
            Mhi=max(point)
            IMF_weight += N(Mhi,alpha) - N(Mlo,alpha)
        self.IMF_weight = IMF_weight







#example of how to run the code:
        
if __name__ == '__main__':
    p000=ISOCMD('/data/MIST/MIST_v1.0_UBVRIplus/MIST_v1.0_feh_p0.00_afe_p0.0_vvcrit0.4_UBVRIplus.iso.cmd')

    #isochrone logAge=9.6, excluding TP- and post-AGB
    iso=p000.isocmds[92]
    iso=iso[iso['EEP']<=808]

    H=HessDiagram(0.5,8,2.5,23,100,100)
    H.add_iso(iso,DM=10)
    H.do_mass_intervals()
    H.do_IMF(alpha=2.35)

    #now make a plot of the Hess diagram
    figure(1,figsize=(8,8))
    for box in H.boxes:
        if box.intersects:
            scatter(box.center[0],box.center[1], c=box.IMF_weight, marker='s', s=20, vmin=0, vmax=0.27, edgecolors='none', cmap='viridis')
    ylim(22,8)
    xlabel('V-I')
    ylabel('V')
    colorbar()
    show()
