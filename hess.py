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
        

p000=ISOCMD('/data/MIST/MIST_v1.0_UBVRIplus/MIST_v1.0_feh_p0.00_afe_p0.0_vvcrit0.4_UBVRIplus.iso.cmd')

#isochrone example 9.6 
tmp_iso=p000.isocmds[92]
iso=tmp_iso[tmp_iso['EEP']<=808]


if False:
    color=iso['Bessell_B']-iso['Bessell_V']
    mag=iso['Bessell_V']
    iso_string=LineString(zip(color,mag))

    #define CMD grid
    NX=100
    NY=100
    XMIN,XMAX = -0.5,2.5
    YMIN,YMAX = -6,20
    color_grid = linspace(XMIN,XMAX,NX)
    mag_grid = linspace(YMIN,YMAX,NY)

    #list of hess boxes
    hessboxes=[]
    for x in range(NX-1):
        xmin, xmax = color_grid[x], color_grid[x+1]
        for y in range(NY-1):
            ymin, ymax =mag_grid[y], mag_grid[y+1]
            p0=(xmin,ymin)
            p1=(xmax,ymin)
            p2=(xmax,ymax)
            p3=(xmin,ymax)
            p4=(xmin,ymin)
            l=LineString((p0,p1,p2,p3,p4))
            if iso_string.intersects(l):
                c=array(iso_string.intersection(l))
                h=HessBox(xmin,ymin,xmax,ymax,c)
                h.mass_intervals=get_mass_intervals(h)
                hessboxes.append(h)

    print len(hessboxes)

    #for i,hb in enumerate(hessboxes):
    #    masses=get_mass_intervals(hb)
    #    print i, masses
