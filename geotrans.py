###################################################
#MODULES REQUIRED
###################################################
from matplotlib import pyplot as plt
from matplotlib import patches as pat
from cmath import sqrt as csqrt,phase
from scipy import constants as const
from scipy.optimize import newton,brentq,fsolve
from scipy.integrate import quad as integrate
from sys import argv,exit
from os import system
from numpy import *
from lmfit import minimize,Parameters,Parameter,report_fit

###################################################
#MACROS
###################################################

#//////////////////////////////
#MACROS
#//////////////////////////////
RAD=180/pi
DEG=pi/180
RAND=random.rand
TAB="\t"
BARL="*"*50+"\n"
RBAR="\n"+"*"*50

MAG=lambda P:sqrt(sum(P**2))
AR=lambda x,y:array([x,y])
AR3=lambda x,y,z:array([x,y,z])
ARCTAN=lambda num,den:mod(arctan2(num,den),2*pi)
EQUAL=lambda x,y:abs(x-y)<=ZERO
DENSITY=lambda M,R:M/(4*pi/3*R**3)
def VERB(routine):print BARL,routine,RBAR

#//////////////////////////////
#NUMERICAL CONSTANTS
#//////////////////////////////
#MACHINE PRECISION
ZERO=finfo(float).eps

#TOLERANCE FOR INTERSECTION POINTS
INTERTOL=1E-10

#TOLERANCE FOR SELECTING INTERSECTION POINTS
INTERFUNTOL=1E-3

#TOLERANCE FOR ANGLE DETECTION IN INTERSECTION
ANGLETOL=3*DEG

#TOLERANCE FOR CALCULATING DISTANCE
DISTANCETOL=1E-15

#TOLERANCE FOR REJECTING RINGS FROM AREA CALCULATION
NORINGTOL=1E-5

#OTHER
IMAGTOL=1E-5
FIGTOL=8E-3

#//////////////////////////////
#ALIASES
#//////////////////////////////
INSIDE=+1
OUTSIDE=-1

INGRESS=-1
EGRESS=+1

CLOSEST=-1
FARTHEST=+1

#//////////////////////////////
#DEBUGGING
#//////////////////////////////
VERBOSE=[0]*10
VERBOSE[0]=0
VERBOSE[1]=1
VERBOSE[2]=1
VERBOSE[3]=1
for i in xrange(10):
    if VERBOSE[i]==0:
        VERBOSE[i:]=[0]*(10-i)
        break

#########################################
#CONSTANTS
#########################################
#////////////////////
#UNITS AND MULTIPLES
#////////////////////
KM=1E3
DAY=const.day
YEAR=const.year
HOUR=const.hour
MINUTE=const.minute
LIGHTY=const.light_year
PARSEC=3.2616*LIGHTY

NANO=1E-9
KILO=1E3
MEGA=1E6
GIGA=1E9

#////////////////////
#PLANETARY
#////////////////////
RJUP=69911.0*KM
MJUP=1.898E27 #kg
RSAT=58232.0*KM
MSAT=5.6846E26 #kg
RSAT_BRING=92000.0*KM
RSAT_ARING=136775.0*KM
RSUN=696342.0*KM
MSUN=1.98855E30 #kg
LSUN=3.846E26 #W
TSUN=5778.0 #K
RHOSUN=MSUN/(4*pi/3*RSUN**3)

#////////////////////
#PHYSICAL CONSTANTS
#////////////////////
AU=149597871.0*KM
GCONST=const.G
HP=const.h
KB=const.k
CSPEED=const.c

###################################################
#DATA TYPES
###################################################
class Point(object):
    def __init__(self,pos,fig1,fig2):
        self.pos=pos
        self.fig1=fig1
        self.fig2=fig2
        self.name='P'
    def __str__(self):
        s="("
        s+="Pos = "+str(self.pos)+","
        s+="Fig1 = "+str(self.fig1)+","
        s+="Fig2 = "+str(self.fig2)
        s+=")"
        return s

class Figure(object):
    def __init__(self,C,a,b,ct,st,name):
        self.C=C
        self.a=a
        self.b=b
        self.cost=ct
        self.sint=st
        self.name=name
    def __str__(self):
        s="("
        s+="C="+str(self.C)+","
        s+="a="+str(self.a)+","
        s+="b="+str(self.b)+","
        s+="cos(t)="+str(self.cost)+","
        s+="sin(t)="+str(self.sint)+","
        s+="name="+str(self.name)
        s+=")"
        return s
    def area(self):
        return pi*self.a*self.b

FNULL=Figure(AR(0,0),0,0,1,0,'')
FONES=Figure(AR(0,0),1.0,1.0,1.0,0.0,'')

class Orbit(object):
    def __init__(self,a,e,P,Mos):
        self.a=a
        self.e=e
        self.P=P
        self.n=2*pi/P
        self.b=a*sqrt(1-e*e)
        self.Mos=Mos
    def __str__(self):
        s="("
        s+="a = "+str(self.a)+","
        s+="e = "+str(self.e)+","
        s+="P = "+str(self.P)
        s+=")"
        return s

class dict2obj(object):
    """
    Class that allows the conversion from dictionary to class-like
    objects
    """
    def __init__(self,dic={}):self.__dict__.update(dic)
    def __add__(self,other):
        self.__dict__.update(other.__dict__)
        return self

def copyObject(obj):
    """
    Create a copy of an object
    """
    new=dict()
    dic=obj.__dict__
    for key in dic.keys():
        var=dic[key]
        new[key]=var
        try:
            var.__dict__
            inew=copyObject(var)
            new[key]=inew
        except:
            pass
    nobj=dict2obj(new)
    return nobj

###################################################
#CONFIGURATION
###################################################
#//////////////////////////////
#DATA RELATED
#//////////////////////////////
def commonFigs(P1,P2):
    """
    Get the figures common to points P1 and P2
    """
    Fs=list(set([P1.fig1,P1.fig2])&set([P2.fig1,P2.fig2]))
    return Fs

def figNames(Fs):
    """
    Get a list of the names of figures Fs
    """
    s=""
    for F in Fs:s+=F.name+","
    s=s.strip(",")
    return s

def pointNames(Ps):
    """
    Get a list of the names of points Ps
    """
    s=""
    for P in Ps:s+=P.name+","
    s=s.strip(",")
    return s

def toPoint(v): 
    """
    Convert array v to Point P
    """
    P=Point(v,FNULL,FNULL)
    return P

#//////////////////////////////
#BASIC
#//////////////////////////////
def randomVal(xmin,xmax):
    x=xmin+(xmax-xmin)*RAND()
    return x

def rotMat(axis,theta):
    """
    Rotation in 3d using Euler-Rodrigues formula
    """
    axis=asarray(axis)
    axis=axis/sqrt(dot(axis,axis))
    a=cos(theta/2)
    b,c,d=-axis*sin(theta/2)
    aa,bb,cc,dd=a*a,b*b,c*c,d*d
    bc,ad,ac,ab,bd,cd=b*c,a*d,a*c,a*b,b*d,c*d
    M=array([[aa+bb-cc-dd,2*(bc+ad),2*(bd-ac)],
             [2*(bc-ad),aa+cc-bb-dd,2*(cd+ab)],
             [2*(bd+ac),2*(cd-ab),aa+dd-bb-cc]])
    return M

def rotTrans(r,cost,sint,b):
    """
    Rotate and translate vector r by angle t and displ. b
    """
    M=array([[cost,-sint],[sint,cost]])
    rp=dot(M,r)+b
    return rp

def realArray(zs):
    """
    Get real parts of zs if they are mostly (IMAGTOL) real
    """
    ims=imag(zs)
    return real(zs[abs(ims)<=IMAGTOL])

def roundComplex(zs):
    """
    Round complex numbers zs to IMAGTOL tolerance
    **OPTIMIZE
    """
    for i in xrange(len(zs)):
        if abs(zs[i].imag)<=IMAGTOL:zs[i]=real(zs[i])
    return zs

def sortPolygonVertices(Ps):
    """
    Sort the vertices of a polygon counter clockwise.
    **OPTIMIZE
    """
    Np=len(Ps)
    if Np==0:return Ps
    C=AR(0,0)
    for P in Ps:C=C+P.pos
    C=C/Np
    i=0
    qs=[]
    for P in Ps:
        dr=P.pos-C
        q=ARCTAN(dr[1],dr[0])
        qs+=[q]
        i+=1
    ies=argsort(qs)
    return Ps[ies]

def FIGUREAREA(F):
    """
    Compute the area of figure F
    """
    if F.b/F.a<NORINGTOL:return 0
    else:return pi*F.a*F.b

def ellipseCoefficients(F):
    """
    Coefficients of the ellipse polynomia
    See Zuluaga et al. (2014)
    """
    if VERBOSE[3]:VERB("ellipseCoefficients")
    C=F.C
    x=C[0];x2=x**2
    y=C[1];y2=y**2
    a=F.a;a2=a**2
    b=F.b;b2=b**2

    v1=1-b2/a2
    v2=-2*x*(b2/a2)
    v3=-2*y
    v4=x2*(b2/a2)+y2-b2+(b2/a2)
    v8=-2*x*(b2/a2)
    v10=-2*x*(b2/a2)
    if VERBOSE[3]:vs=array([v1,v2,v3,v4,v8,v10])
    if VERBOSE[3]:print "Auxiliar variables: ",vs
    
    e=v2*v10-v4**2
    d=-2*v3*v4
    c=-v2*v8-v3**2-2*v1*v4
    b=-2*v1*v3
    a=-v1**2
    if VERBOSE[3]:print "Coefficients: ",a,b,c,d,e

    return a,b,c,d,e

def thirdPolyRoot(a,b,c,d):
    """
    Roots of: a x^3 + b x^2 + c x + d = 0
    See: Wikipedia/Cubic Equation
    """
    if VERBOSE[3]:VERB("thirdPolyRoot")

    #AUXILIAR VARIABLES
    D=18*a*b*c*d-4*b**3*d+b**2*c**2-4*a*c**3-27*a**2*d**2
    D0=b**2-3*a*c
    D1=2*b**3-9*a*b*c+27*a**2*d
    D2=D1**2-4*D0**3
    C0=((D1+csqrt(D2))/2)**(1./3)
    u1=1;u2=(-1+3**0.5*1j)/2;u3=(-1-3**0.5*1j)/2
    
    #ROOTS
    xs=[]
    xs+=[-1/(3*a)*(b+u1*C0+D0/(u1*C0))]
    xs+=[-1/(3*a)*(b+u2*C0+D0/(u2*C0))]
    xs+=[-1/(3*a)*(b+u3*C0+D0/(u3*C0))]

    if VERBOSE[3]:print "Roots poly-3:",xs
    return array(xs)

def fourthPolyRoots(a,b,c,d,e):
    """
    Roots of: a x^4 + b x^3 + c x^2 + d x + e = 0
    See Wolfram MathWorld/Quartic Equation
    **OPTIMIZE
    """
    if VERBOSE[3]:VERB("fourthPolyRoots")

    #NORMALIZE
    a0=float(e)/a
    a1=float(d)/a
    a2=float(c)/a
    a3=float(b)/a
    a4=1.0
    if VERBOSE[3]:
        print "Polynomial (%f) x^4 + (%f) x^3 + (%f) x^2 + (%f) x + (%f) = 0"%(a4,a3,a2,a1,a0)

    #RELATED THIRD POLYNOMIAL
    a=1.0
    b=-a2
    c=a1*a3-4*a0
    d=4*a2*a0-a1**2-a3**2*a0
    xs=thirdPolyRoot(a,b,c,d)
    xs=realArray(xs)
    if VERBOSE[3]:print "Real roots poly-3= ",xs

    #CHOOSE REAL ROOTS
    Rs=array([csqrt(0.25*a3**2-a2+x) for x in xs])

    if VERBOSE[3]:print "Auxiliar Rs = ",Rs
    Rs=roundComplex(Rs)
    if VERBOSE[3]:print "Rounded Auxiliar Rs = ",Rs
    
    #SOLVE POLYNOMIAL
    Ds=[];Es=[]

    #**REMOVE THIS FOR LOOP
    for x,R in zip(xs,Rs):
        if abs(R)<IMAGTOL:
            if VERBOSE[3]:print "Auxiliar R is effectively 0."
            Ds+=[csqrt(0.75*a3**2-2*a2+2*csqrt(x**2-4*a0))]
            Es+=[csqrt(0.75*a3**2-2*a2-2*csqrt(x**2-4*a0))]
        else:
            if VERBOSE[3]:print "Auxiliar R is effectively different from 0."
            Ds+=[csqrt(0.75*a3**2-R**2-2*a2+0.25*(4*a3*a2-8*a1-a3**3)/R)]
            Es+=[csqrt(0.75*a3**2-R**2-2*a2-0.25*(4*a3*a2-8*a1-a3**3)/R)]
    Ds=array(Ds);Es=array(Es)

    z1s=-0.25*a3+0.5*Rs+0.5*Ds
    z2s=-0.25*a3+0.5*Rs-0.5*Ds
    z3s=-0.25*a3-0.5*Rs+0.5*Es
    z4s=-0.25*a3-0.5*Rs-0.5*Es

    zs=array(zip(z1s,z2s,z3s,z4s))
    if VERBOSE[3]:print "Root sets poly-4: ",zs

    return zs

def ellipseRadius(F,cost,sint):
    """
    Radius of the ellipse at angle t from center
    See Wikipedia/Ellipse
    """
    return F.a*F.b/((F.b*cost)**2+(F.a*sint)**2)**0.5

def ellipseRadiusE(a,e,f):
    return a*(1-e**2)/(1+e*cos(f))

def ellipsePoint(F,cost,sint):
    """
    Coordinates of an ellipse point at angle t from center
    """
    r=ellipseRadius(F,cost,sint)
    return AR(r*cost,r*sint)+F.C

def ellipsePointEcc(F,cosE,sinE):
    """
    Coordinates of an ellipse point at angle t from center
    """
    return AR(F.a*cosE,F.b*sinE)+F.C

def pointInFigure(F,P):
    """
    Evaluate if point P is in figure
    Positive: Inside
    Negative: Outside
    ** OPTIMIZE
    """
    if EQUAL(F.a,F.b):
        d=P.pos-F.C
        r=MAG(d)
        rf=F.a
        return (rf-r)/r
    ct=F.cost
    st=F.sint
    d=P.pos-F.C
    d=rotTrans(d,ct,-st,AR(0,0))
    r=MAG(d)
    if EQUAL(r,ZERO):return +F.a
    cost=d[0]/r
    sint=d[1]/r
    rf=ellipseRadius(F,cost,sint)
    return (rf-r)/r

#//////////////////////////////
#INTERSECTION
#//////////////////////////////
def qCircleCircle(F1,F2):
    """
    Determine if two circles intersect
    Returns: 
      +1: Disjoint and contained
       0: Intersection
      -1: Disjoint and not contained
    """
    R1=F1.a;R2=F2.a
    D=MAG(F2.C-F1.C)
    if D>=(R1+R2):return -1
    if D<=abs(R1-R2):return 1
    return 0

def cIc(F1,F2):
    """
    Computes the 2 points of intersection of 2 (nonconcentric) circles
    It is assumed that F1<F2
    """
    q=qCircleCircle(F1,F2) or 0
    if q:return \
            Point(AR(123*q,123*q),F1,F2),\
            Point(AR(123*q,123*q),F1,F2)

    r=F1.a
    C=F1.C-F2.C
    R=F2.a
    x=C[0];y=C[1]
    D=MAG(C)
    r2=r**2;R2=R**2;D2=D**2
    dis=sqrt(\
        -D2**2-2*D2*r2+2*D2*R2-\
             r2**2+2*r2*R2+4*r2*x**2+\
             4*r2*y**2-R2**2)
    det=-2*r*y
    den=D2+r2-2*r*C[0]-R2
    t1=mod(2*arctan2(+dis+det,den),2*pi)
    t2=mod(2*arctan2(-dis+det,den),2*pi)
    return \
        Point(C+AR(r*cos(t1),r*sin(t1)),F1,F2),\
        Point(C+AR(r*cos(t2),r*sin(t2)),F1,F2)

def eIcAnalytical(F1,F2):
    """
    Computes the 4 intersection points between an ellipse
    (F1) and an unitary circle (F2).
    """

    #NON-CONCENTRIC ELLIPSE AND UNITARY CIRCLE (F1<F2)
    if VERBOSE[3]:VERB("eIc")
    if VERBOSE[3]:print "Figures: ",F1,F2
    
    C=F1.C
    x=C[0]
    y=C[1]
    a=F1.a
    b=F1.b
    qin=sign(pointInFigure(F2,toPoint(C)))
    if VERBOSE[3]:print "qin = ",qin

    #CHECK IF ELLIPSE F1 IS A CIRCLE 
    if EQUAL(a,b):
        if VERBOSE[3]:print "Ellipse is a Circle."
        return cIc(F1,F2)

    #GET THE ELLIPSE COEFFICIENTS
    ac,bc,cc,dc,ec=ellipseCoefficients(F1)

    #SOLVE THE CORRESPONDING FOURTH ORDER EQUATION
    ys=fourthPolyRoots(ac,bc,cc,dc,ec)
    ys=realArray(ys)
    if VERBOSE[3]:print "Intersetion in y (Real) = ",ys
    ys=ys[abs(ys)<1]
    if VERBOSE[3]:print "Intersection in y (In range) = ",ys
    #IF NO SOLUTION, NO INTERSECTION
    if len(ys)==0:
        return \
            Point(AR(123*qin,123*qin),F1,F2),\
            Point(AR(123*qin,123*qin),F1,F2),\
            Point(AR(123*qin,123*qin),F1,F2),\
            Point(AR(123*qin,123*qin),F1,F2)
    
    #SOLUTION
    qtrad=0
    try:
        xs=[]
        alpha1=0
        alpha2=1
        beta1=-2*x/a**2
        beta2=1/a**2
        det=(alpha1*beta2-alpha2*beta1)
        if VERBOSE[3]:print "Determinant: ",det
        if abs(det)>=1E-13:
            for yp in ys:
                alpha0=yp**2-1
                beta0=yp**2/b**2-2*y*yp/b**2+x**2/a**2+y**2/b**2-1
                xs+=[(alpha2*beta0-alpha0*beta2)/det]
            xs=array(xs)
        else:qtrad=1
    except ValueError as error:
        if VERBOSE[3]:print "Error:",error
        qtrad=2

    if qtrad:
        if VERBOSE[3]:print "Using traditional formula for x (reason = %d)"%qtrad
        xs=[]
        #**OPTIMIZE
        for y in ys:
            x=sqrt(1-y**2)
            xs+=[x,-x]
        xs=array(xs)
        
    if VERBOSE[3]:print "Intersection in x (qtrad = %d) = "%qtrad,xs

    #CHOOSE THE ACTUAL SOLUTION
    lys=len(ys)
    if lys>2:
        #IF MORE THAN 2 POINTS
        qps=[]
        i=0
        for x,y in zip(xs,ys):
            if VERBOSE[3]:print "Testing couple: ",x,y
            qps+=[pointInFigure(F1,toPoint(AR(x,y)))]
            if VERBOSE[3]:print "Value: ",qps[i]
            i+=1

        #GET THE POINTS CLOSER TO THE FIGURE
        ysave=ys
        qps=abs(array(qps))
        cond=(qps<=FIGTOL)
        qps=qps[cond];xs=xs[cond];ys=ys[cond]
        iargs=qps.argsort()

        #MORE THAN TWO POINTS ARE REALLY CLOSE TO THE FIGURE
        if len(qps)>2:
            Ps=\
                Point(AR(xs[iargs[0]],ys[iargs[0]]),F1,F2),\
                Point(AR(xs[iargs[1]],ys[iargs[1]]),F1,F2),\
                Point(AR(xs[iargs[2]],ys[iargs[2]]),F1,F2),\
                Point(AR(xs[iargs[3]],ys[iargs[3]]),F1,F2)

        #ONLY 2 POINTS ARE REALLY CLOSE TO THE FIGURE
        else:
            Ps=\
                Point(AR(xs[iargs[0]],ys[iargs[0]]),F1,F2),\
                Point(AR(xs[iargs[1]],ys[iargs[1]]),F1,F2),\
                Point(AR(123,123),F1,F2),\
                Point(AR(123,123),F1,F2)
    else:
        #ONLY 2 POINTS ARE REALLY CLOSE TO THE FIGURE
        Ps=\
            Point(AR(xs[0],ys[0]),F1,F2),\
            Point(AR(xs[1],ys[1]),F1,F2),\
            Point(AR(-123,-123),F1,F2),\
            Point(AR(-123,-123),F1,F2)
    
    if VERBOSE[3]:print "Points (%d): "%lys,Ps
    return Ps

def cIe(F1,F2):
    """
    Intersection of an ellipse (F1) with a circle (F2)
    See: Wolfram MathWorld/Circle Ellipse Intersection
    """
    C=F1.C;a=F1.a;b=F1.b;R=F2.a

    if R<b:
        #CIRCLE IS COMPLETELY CONTAINED WITHIN ELLIPSE
        return \
            Point(AR(-123,-123),F1,F2),\
            Point(AR(-123,-123),F1,F2),\
            Point(AR(-123,-123),F1,F2),\
            Point(AR(-123,-123),F1,F2)

    #COMPUTE SOLUTION
    a2=a**2;b2=b**2;R2=R**2
    x=a*sqrt((R2-b2)/(a2-b2))
    y=b*sqrt((a2-R2)/(a2-b2))
    
    return\
        Point(AR(+x,+y)+C,F1,F2),\
        Point(AR(-x,+y)+C,F1,F2),\
        Point(AR(-x,-y)+C,F1,F2),\
        Point(AR(+x,-y)+C,F1,F2)
    
#//////////////////////////////
#AREAS
#//////////////////////////////
def ellipseSegment(F,P1,P2):
    """
    Compute the area of an ellipse segment between points P1 and P2
    """
    C=F.C
    a=F.a
    b=F.b
    dr1=P1.pos-C
    dr2=P2.pos-C
    t1=ARCTAN(dr1[1]/b,dr1[0]/a)
    t2=ARCTAN(dr2[1]/b,dr2[0]/a)
    dt=abs(t2-t1)
    if dt>pi:dt=2*pi-dt
    A=a*b/2*(dt-sin(dt))
    return A

def ellipseSegmentOriented(F,P1,P2,sgn=+1):
    """
    Compute the area of an ellipse segment between points P1 and P2

    sgn:

      +1: This is a minor area figure.  The largest area contained by
          the two complimentary segments will be returned.

      -1: This is a major area figure.  The largest area contained by
          the two complimentary segments will be returned.

    """
    
    if VERBOSE[2]:VERB("ellipseSegmentOriented")
    if VERBOSE[2]:print "Ellipse Segment between points (sgn = %d): "%sgn,P1,P2
    q1=pointInFigure(F,P1)
    q2=pointInFigure(F,P2)
    if VERBOSE[2]:print "Condition point 1:",q1
    if VERBOSE[2]:print "Condition point 2:",q2
    if q1<-FIGTOL or q2<-FIGTOL:
        return 0

    C=F.C
    a=F.a
    b=F.b
    
    #ANGLES
    dr1=P1.pos-C
    dr2=P2.pos-C
    t1=ARCTAN(dr1[1]/b,dr1[0]/a)
    t2=ARCTAN(dr2[1]/b,dr2[0]/a)
    dt=abs(t2-t1)
    if VERBOSE[3]:print "t1,t2,dt = ",t1*RAD,t2*RAD,dt*RAD

    #FAKE SEGMENT
    if t1==t2:
        if pointInFigure(F,P1)<ZERO:
            if VERBOSE[3]:print "Segment closed."
            if sgn<0:
                if VERBOSE[3]:print "Minor curve: All area."
                return pi*a*b
            if sgn>0:
                if VERBOSE[3]:print "Major curve: No area."
                return 0.0
        if t1>pi:return -123
        if t1<pi:return +123

    if dt>pi:dt=2*pi-dt
    A=a*b/2*(dt-sin(dt))

    if sgn<0:
        #SMALL FIGURE CONDITION
        R=P2.pos-P1.pos
        det=(C[0]*R[1]-C[1]*R[0])
        if EQUAL(det,ZERO):det=1
        lam=(P1.pos[0]*R[1]-P1.pos[1]*R[0])/det
        if lam>1:A=pi*a*b-A

    return A

def planeTriangle(Ps):
    """
    Area of a triangle via Heron's Formula
    """
    P1,P2,P3=Ps
    S=[MAG(P1.pos-P2.pos),MAG(P2.pos-P3.pos),MAG(P1.pos-P3.pos)]
    s=sum(S)/2
    A=sqrt(s*prod(s-S))
    return A

def planeQuad(Ps):
    """
    Area of a quadrangle via BRETSCHNEIDER'S FORMULA
    """
    P1,P2,P3,P4=Ps
    a=MAG(P1.pos-P2.pos)
    b=MAG(P2.pos-P3.pos)
    c=MAG(P3.pos-P4.pos)
    d=MAG(P4.pos-P1.pos)
    S=array([a,b,c,d])
    p=MAG(P2.pos-P4.pos)
    q=MAG(P1.pos-P3.pos)
    s=sum(S)/2
    Aq=sqrt(prod(s-S)-0.25*(a*c+b*d+p*q)*(a*c+b*d-p*q))
    return Aq

def montecarloArea(Fs,oper,Npoints=1E3):
    """
    Calculate an area using montecarlo integration.

    Area calculated is that contained in figures Fs according to
    inlcusion/exclusion operations in oper.
    """
    #GET MINIMUM AND MAXIMUM VALUES FOR RANDOM POINTS
    Es=linspace(0,2*pi,100)
    cosEs=cos(Es)
    sinEs=sin(Es)
    xs=array([])
    ys=array([])
    for F in Fs[1:]:
        C=F.C
        cost=F.cost;sint=F.sint
        x=C[0];y=C[1]
        a=F.a
        b=F.b
        xs=concatenate((xs,(a*cosEs*cost-b*sinEs*sint)+x))
        ys=concatenate((ys,(a*cosEs*sint+b*sinEs*cost)+y))
    xmin=xs.min();xmax=xs.max()
    ymin=ys.min();ymax=ys.max()
    if VERBOSE[0]:print "Montecarlo xmin,xmax = ",xmin,xmax
    if VERBOSE[0]:print "Montecarlo ymin,ymax = ",ymin,ymax

    #COMPUTE AREA
    xs=[];ys=[]
    i=0
    c=0
    while i<Npoints:
        #GENERATE POINT
        P=Point(AR(xmin+(xmax-xmin)*RAND(),ymin+(ymax-ymin)*RAND()),FNULL,FNULL)
        cond=True
        o=0
        for F in Fs:
            cond=(cond and ((oper[o]*pointInFigure(F,P))>=0))
            o+=1
        if cond:
            c+=1
            xs+=[P.pos[0]]
            ys+=[P.pos[1]]
        i+=1
        
    A=(xmax-xmin)*(ymax-ymin)*(float(c)/i)
    if c>0:dA=(1/sqrt(c))*A
    else:dA=0
    return A,dA,xs,ys

def convexTriangle(Ps,shapes=[+1,+1,+1]):
    """
    Area of a convex triangle limited by vertices Ps.  
    
    Shapes Determine which sides should be plane.  If shapes[i]=0 side
    between vertex Pi and Pi+1 should be treated as plane.  If
    shapes[i]=-1 this side wil be concave and substracted from the
    total area.

    **OPTIMIZE
    """
    P1,P2,P3=Ps

    #Side 12
    if shapes[0]==0:
        if VERBOSE[3]:print "P1,P2 is plane."
        A1=0
    else:
        Fs=commonFigs(P1,P2)
        if VERBOSE[3]:print "P1,P2 common figures: ",figNames(Fs)
        try:A1=min([ellipseSegment(F,P1,P2) for F in Fs])
        except IndexError:A1=0
        if VERBOSE[3]:print "A1 = ",A1

    #Side 23
    if shapes[1]==0:
        if VERBOSE[3]:print "P2,P3 is plane."
        A2=0
    else:
        Fs=commonFigs(P2,P3)
        if VERBOSE[3]:print "P2,P3 common figures: ",figNames(Fs)
        try:A2=min([ellipseSegment(F,P2,P3) for F in Fs])
        except IndexError:A2=0
        if VERBOSE[3]:print "A2 = ",A2

    #Side 31
    if shapes[2]==0:
        if VERBOSE[3]:print "P3,P1 is plane."
        A3=0
    else:
        Fs=commonFigs(P3,P1)
        if VERBOSE[3]:print "P3,P1 common figures: ",figNames(Fs)
        try:A3=min([ellipseSegment(F,P3,P1) for F in Fs])
        except IndexError:A3=0
        if VERBOSE[3]:print "A2 = ",A3
        
    Ac=A1+A2+A3
    if VERBOSE[3]:print "Curved Area = ",Ac

    At=planeTriangle(Ps)
    if VERBOSE[3]:print "Plane Area = ",At

    return At+Ac
    
def convexQuad(Ps,shapes=[+1,+1,+1,+1]):
    """
    Area of a convex quadrangle limited by vertices Ps.  
    
    Shapes Determine which sides should be plane.  If shapes[i]=0 side
    between vertex Pi and Pi+1 should be treated as plane.  If
    shapes[i]=-1 this side wil be concave and substracted from the
    total area.

    **OPTIMIZE
    """

    P1,P2,P3,P4=Ps

    if VERBOSE[3]:print "Convex Quadrangle:"

    #Side 12
    if shapes[0]==0:
        if VERBOSE[3]:print "P1,P2 is plane."
        A12=0
    else:
        Fs=commonFigs(P1,P2)
        if VERBOSE[3]:print "P1,P2 common figures: ",figNames(Fs)
        try:A12=min([ellipseSegment(F,P1,P2) for F in Fs])
        except IndexError:A12=0
        if VERBOSE[3]:print "A12 = ",A12

    #Side 23
    if shapes[1]==0:
        if VERBOSE[3]:print "P2,P3 is plane."
        A23=0
    else:
        Fs=commonFigs(P2,P3)
        if VERBOSE[3]:print "P2,P3 common figures: ",figNames(Fs)
        try:A23=min([ellipseSegment(F,P2,P3) for F in Fs])
        except IndexError:A23=0
        if VERBOSE[3]:print "A23 = ",A23

    #Side 34
    if shapes[2]==0:
        if VERBOSE[3]:print "P3,P4 is plane."
        A34=0
    else:
        Fs=commonFigs(P3,P4)
        if VERBOSE[3]:print "P3,P1 common figures: ",figNames(Fs)
        try:A34=min([ellipseSegment(F,P3,P4) for F in Fs])
        except IndexError:A34=0
        if VERBOSE[3]:print "A34 = ",A34

    #Side 41
    if shapes[3]==0:
        if VERBOSE[3]:print "P4,P1 is plane."
        A41=0
    else:
        Fs=commonFigs(P4,P1)
        if VERBOSE[3]:print "P4,P1 common figures: ",figNames(Fs)
        try:A41=min([ellipseSegment(F,P4,P1) for F in Fs])
        except IndexError:A41=0
        if VERBOSE[3]:print "A41 = ",A41
        
    Ac=A12+A23+A34+A41
    if VERBOSE[3]:print "Curved Area = ",Ac
    Aq=planeQuad(Ps)

    if VERBOSE[3]:print "Plane Area = ",Aq
    return Aq+Ac

def convexPolygon(Ps):
    """
    Area of a convex polygon limited by vertices Ps.  Up to an
    hexagon.
    
    **OPTIMIZE
    """
    if VERBOSE[3]:VERB("convexPolygon")
    
    nP=len(Ps)
    if VERBOSE[3]:print "Sides of the polygon: ",nP
    if nP<=1:
        if VERBOSE[3]:print "No points."
        A=0.0
    elif nP==2:
        if VERBOSE[3]:print "2 points: a leaf."
        A=leafArea(Ps)
    elif nP==3:
        if VERBOSE[3]:print "3 points: a triangle."
        A=convexTriangle(Ps)
    elif nP==4:
        if VERBOSE[3]:print "4 points: a quadrangle."
        A=convexQuad(Ps)
    elif nP==5:
        if VERBOSE[3]:print "5 points: a pentagon."
        A1=convexQuad(Ps[:4],shapes=[+1,+1,+1,0])
        A2=convexTriangle((Ps[0],Ps[3],Ps[4]),shapes=[0,+1,+1])
        A=A1+A2
    elif nP==6:
        if VERBOSE[3]:print "6 points: a hexagon."
        A1=convexQuad((Ps[:4]),shapes=[+1,+1,+1,0])
        A2=convexQuad((Ps[0],Ps[3],Ps[4],Ps[5]),shapes=[0,+1,+1,+1])
        A=A1+A2
    else:
        print "An excessive number of polygon sides."
    return A
        
def leafArea(Ps):
    """
    Area of a leaf between points Ps and limited by figures F1 and F2
    common to both boints.
    """
    if VERBOSE[1]:VERB("leafArea")

    #BY DEFAULT FIGURE 1 > FIGURE2
    F1=Ps[0].fig1;AF1=FIGUREAREA(F1)
    F2=Ps[1].fig2;AF2=FIGUREAREA(F2)
    sgn=+1

    if AF1<AF2:
        sgn=-1
        AF=AF1
    else:AF=AF2

    if VERBOSE[1]:print "Figures: ",F1,F2
    if VERBOSE[1]:print "Larger figure (+1 if F1): ",sgn

    #CHECK IF POINTS ARE SPECIAL
    if Ps[0].pos[0]==123:
        #SMALLER FIGURE INSIDE
        if VERBOSE[1]:print "Small figure completely inside."
        return AF
    elif Ps[0].pos[0]==-123:
        #SMALLER FIGURE OUTSIDE
        if VERBOSE[1]:print "Small figure completely outside."
        return 0

    Pspos=[P.pos for P in Ps]
    if VERBOSE[1]:print "Points: ",pointNames(Ps)
    if VERBOSE[1]:print "Point positions: ",Pspos

    A1=ellipseSegmentOriented(F1,Ps[0],Ps[1],sgn=sgn)
    A2=ellipseSegmentOriented(F2,Ps[0],Ps[1],sgn=-sgn)
    if A1==0 and A2==0:return 0
    if VERBOSE[1]:print "Segment 1 Area: ",A1
    if VERBOSE[1]:print "Segment 2 Area: ",A2

    if sgn>0:
        if A2==-123:return FIGUREAREA(F2)
        if A2==123:return 0
    else:
        if A1==-123:return FIGUREAREA(F1)
        if A1==123:return 0

    Al=A1+A2
    return Al

#//////////////////////////////
#RINGED PLANET AREA
#//////////////////////////////
def ringedPlanetArea(S):
    """
    Area of a planet with its rings
    """
    Planet=S.Planet
    Ringe=S.Ringext
    Ringi=S.Ringint

    C=Planet.C
    ct=Ringe.cost
    st=Ringe.sint
    Rp=Planet.a
    Rea=Ringe.a
    Reb=Ringe.b
    Ria=Ringi.a
    Rib=Ringi.b

    #ONLY PLANET
    if Rea==0 or Reb/Rea<NORINGTOL:return pi*Rp**2
    
    #PUTTING RINGS HORIZONTAL
    Ca=rotTrans(C,ct,-st,AR(0,0))

    #NORMALIZING POSITION
    Cg=abs(Ca)

    #FIGURES
    Star=Figure(AR(0.0,0.0),1.0,1.0,1.0,0.0,'Star')
    Planet=Figure(Cg,Rp,Rp,1.0,0.0,'Planet')
    Ringe=Figure(Cg,Rea,Reb,1.0,0.0,'Ringe')
    Ringi=Figure(Cg,Ria,Rib,1.0,0.0,'Ringi')

    #AREAS
    Asp=pi*Rp**2
    Asre=pi*Rea*Reb
    Asri=pi*Ria*Rib

    #INTERSECTION POINTS
    Ppre1,Ppre2,Ppre3,Ppre4=cIe(Ringe,Planet)
    Ppri1,Ppri2,Ppri3,Ppri4=cIe(Ringi,Planet)

    #EXTERNAL RING
    Ps=array([Ppre1,Ppre2,Ppre3,Ppre4])
    Fs=[Planet,Planet,Planet,Planet]
    Qs=array([pointInFigure(F,P) for F,P in zip(Fs,Ps)])
    Pine=sortPolygonVertices(Ps[Qs>-FIGTOL])

    #INTERNAL RING
    Ps=array([Ppri1,Ppri2,Ppri3,Ppri4]) 
    Fs=[Planet,Planet,Planet,Planet]
    Qs=array([pointInFigure(F,P) for F,P in zip(Fs,Ps)])
    Pini=sortPolygonVertices(Ps[Qs>-FIGTOL])

    #AREAS
    Pcusp=toPoint(AR(Planet.C[0],Planet.C[1]+Planet.b))
    if len(Pine)==0:
        qcusp=pointInFigure(Ringe,Pcusp)
        if qcusp<-FIGTOL:Asrec=0.0
        else:Asrec=Asp
    else:Asrec=convexPolygon(Pine)
    if len(Pini)==0:
        qcusp=pointInFigure(Ringi,Pcusp)
        if qcusp<-FIGTOL:Asric=0.0
        else:Asric=Asp
    else:Asric=convexPolygon(Pini)

    #RINGED PLANET AREA
    try:
        #WITH FINITE OPACITY
        Aringed=Asp+S.block*((Asre-Asrec)-(Asri-Asric))
    except:
        #ASSUMING INFINITE OPACITY
        Aringed=Asp+Asre-Asri-Asrec+Asric

    return Aringed

#//////////////////////////////
#TRANSIT AREA
#//////////////////////////////
def transitArea(S):
    """
    Compute transit area of planet with his rings (Ringe, Ringi) over
    a star with unitary radius.
    ** OPTIMIZE A LOT
    """
    #////////////////////////////////////////
    #INTERSECTION METHOD
    #////////////////////////////////////////
    #eIc=eIcAnalytical
    Planet=S.Planet
    Ringe=S.Ringext
    Ringi=S.Ringint

    #////////////////////////////////////////
    #BASIC PROPERTIES
    #////////////////////////////////////////
    C=Planet.C
    ct=Ringe.cost
    st=Ringe.sint
    Rp=Planet.a
    Rea=Ringe.a
    Reb=Ringe.b
    Ria=Ringi.a
    Rib=Ringi.b

    #////////////////////////////////////////
    #FIGURES ORIENTATION FIX
    #////////////////////////////////////////
    #RINGS HORIZONTAL
    Ca=rotTrans(C,ct,-st,AR(0,0))
    #NORMALIZING POSITION
    Cg=abs(Ca)
    #FIGURES IN NORMAL POSITION
    Star=Figure(AR(0.0,0.0),1.0,1.0,1.0,0.0,'Star')
    Planet=Figure(Cg,Rp,Rp,1.0,0.0,'Planet')
    Ringe=Figure(Cg,Rea,Reb,1.0,0.0,'Ringe')
    Ringi=Figure(Cg,Ria,Rib,1.0,0.0,'Ringi')
    Feqs=[Planet,Ringe,Ringi]

    #INTERSECTION POINTS
    Psa=[]

    #////////////////////////////////////////
    #STAR AND PLANET
    #////////////////////////////////////////
    Psp1,Psp2=cIc(Planet,Star)
    Psp1.name='Psp1';Psp2.name='Psp2';
    Psa+=[Psp1,Psp2]
    Asp=leafArea((Psp1,Psp2))

    #////////////////////////////////////////
    #IF NO RINGS (i=90) USE ONLY PLANET
    #////////////////////////////////////////
    if Ringe.a==0 or Reb/Rea<NORINGTOL:
        return Asp,Asp,0,0,0,0,Psa,Feqs

    #////////////////////////////////////////
    #PLANET AND RINGS
    #////////////////////////////////////////
    Ppre1,Ppre2,Ppre3,Ppre4=cIe(Ringe,Planet)
    Ppre1.name='Ppre1';Ppre2.name='Ppre2';
    Ppre3.name='Ppre3';Ppre4.name='Ppre4';
    Psa+=[Ppre1,Ppre2,Ppre3,Ppre4]
    Ppri1,Ppri2,Ppri3,Ppri4=cIe(Ringi,Planet)
    Ppri1.name='Ppri1';Ppri2.name='Ppri2';
    Ppri3.name='Ppri3';Ppri4.name='Ppri4';
    Psa+=[Ppri1,Ppri2,Ppri3,Ppri4]

    #////////////////////////////////////////
    #STAR AND EXTERNAL RINGS
    #////////////////////////////////////////
    qine=0;qoute=0
    Psre1,Psre2,Psre3,Psre4=eIc(Ringe,Star)
    Psre1.name='Psre1';Psre2.name='Psre2';
    Psre3.name='Psre3';Psre4.name='Psre4';
    Psa+=[Psre1,Psre2,Psre3,Psre4]
    Psn=[str(P) for P in Psa]
    Psre=array([Psre1,Psre2,Psre3,Psre4])
    qsre=array([P.pos[0] for P in Psre])
    if len(qsre[qsre==123])==4:qine=1
    if len(qsre[qsre==-123])==4:qoute=1
    if (qine+qoute)==0:
        Fsre=[Ringe,Ringe,Ringe,Ringe]
        Qsre=array(\
            [pointInFigure(F,P)\
                 for F,P in zip(Fsre,Psre)])
        Psre=sortPolygonVertices(Psre[Qsre>-FIGTOL])
        Asre=convexPolygon(Psre)
    else:
        Asre=qine*FIGUREAREA(Ringe)
    
    #////////////////////////////////////////
    #STAR AND INTERNAL RINGS
    #////////////////////////////////////////
    qini=0;qouti=0
    Psri1,Psri2,Psri3,Psri4=eIc(Ringi,Star)
    Psri1.name='Psri1';Psri2.name='Psri2';
    Psri3.name='Psri3';Psri4.name='Psri4';
    Psa+=[Psri1,Psri2,Psri3,Psri4]
    Psn=[str(P) for P in Psa]
    Psri=array([Psri1,Psri2,Psri3,Psri4])
    qsri=array([P.pos[0] for P in Psri])
    if len(qsri[qsri==123])==4:qine=1
    if len(qsri[qsri==-123])==4:qoute=1
    if (qini+qouti)==0:
        Fsri=[Ringe,Ringe,Ringe,Ringe]
        Qsri=array(\
            [pointInFigure(F,P)\
                 for F,P in zip(Fsri,Psri)])
        Psri=sortPolygonVertices(Psri[Qsri>-FIGTOL])
        Asri=convexPolygon(Psri)
    else:
        Asri=qini*FIGUREAREA(Ringi)

    #////////////////////////////////////////
    #COMMON POINTS
    #////////////////////////////////////////
    #EXTERNAL RING
    Ps=array([Ppre1,Psre1,Ppre2,Psp2,Ppre3,
              Psre2,Ppre4,Psp1,Psre3,Psre4])
    Fs=[Star,Planet,Star,Ringe,Star,
        Planet,Star,Ringe,Planet,Planet]
    Qs=array([pointInFigure(F,P) for F,P in zip(Fs,Ps)])
    Pine=sortPolygonVertices(Ps[Qs>0])
    line=len(Pine)

    #INTERNAL RING
    Ps=array([Ppri1,Psri1,Ppri2,Psp2,Ppri3,
              Psri2,Ppri4,Psp1,Psri3,Psri4])
    Fs=[Star,Planet,Star,Ringi,Star,
        Planet,Star,Ringi,Planet,Planet]
    Qs=array([pointInFigure(F,P) for F,P in zip(Fs,Ps)])
    Pini=sortPolygonVertices(Ps[Qs>0])

    #////////////////////////////////////////
    #COMMON AREAS
    #////////////////////////////////////////
    if len(Pine)==0 or len(Pini)==0:
        Pcusp=toPoint(AR(Planet.C[0],Planet.C[1]+Planet.b))
        qcusp=pointInFigure(Star,Pcusp)
    if len(Pine)==0:
        if qcusp<0:Asrec=0.0
        else:Asrec=Asp
    else:Asrec=convexPolygon(Pine)
    if len(Pini)==0:
        if qcusp<0:Asric=0.0
        else:Asric=Asp
    else:Asric=convexPolygon(Pini)

    #////////////////////////////////////////
    #TRANSIT AREA
    #////////////////////////////////////////
    try:
        #WITH FINITE OPACITY
        Atrans=Asp+S.block*((Asre-Asrec)-(Asri-Asric))
    except:
        #ASSUMING INFINITE OPACITY
        Atrans=Asp+Asre-Asri-Asrec+Asric

    return Atrans,Asp,Asre,Asri,Asrec,Asric,Psa,Feqs

def transitAreaTime(t,S):
    #UPDATE POSITION
    updatePosition(S,t)
    #AREAS
    Es=transitArea(S)
    At=Es[0]
    return At

def transitAreaTimeFast(t,tcs,Ar,S):
    #ONLY COMPUTE AREA IF AT INGRESS OR EGRESS PHASE
    if t<=tcs[1]:At=0
    elif t<=tcs[2]:At=transitAreaTime(t,S)
    elif t<=tcs[3]:At=Ar
    elif t<=tcs[4]:At=transitAreaTime(t,S)
    else:At=0
    return At

#//////////////////////////////
#TRANSIT OBLATE AREA
#//////////////////////////////
def transitAreaOblate(S):
    """
    Compute transit area of an oblate planet
    """
    #////////////////////////////////////////
    #INTERSECTION METHOD
    #////////////////////////////////////////
    Planet=S.Planet

    #////////////////////////////////////////
    #BASIC PROPERTIES
    #////////////////////////////////////////
    C=Planet.C
    ct=Planet.cost
    st=Planet.sint
    Rpa=Planet.a
    Rpb=Planet.b

    #////////////////////////////////////////
    #FIGURES ORIENTATION FIX
    #////////////////////////////////////////
    #RINGS HORIZONTAL
    Ca=rotTrans(C,ct,-st,AR(0,0))
    #NORMALIZING POSITION
    Cg=abs(Ca)
    #FIGURES IN NORMAL POSITION
    Star=Figure(AR(0.0,0.0),1.0,1.0,1.0,0.0,'Star')
    Planet=Figure(Cg,Rpa,Rpb,1.0,0.0,'Planet')
    Feqs=[Planet]

    #INTERSECTION POINTS
    Psa=[]

    #////////////////////////////////////////
    #STAR AND PLANET
    #////////////////////////////////////////
    qine=0;qoute=0
    Psp1,Psp2,Psp3,Psp4=eIc(Planet,Star)
    Psp1.name='Psp1';Psp2.name='Psp2';
    Psp3.name='Psp3';Psp4.name='Psp4';
    Psa+=[Psp1,Psp2,Psp3,Psp4]
    Psn=[str(P) for P in Psa]
    Psp=array([Psp1,Psp2,Psp3,Psp4])
    qsp=array([P.pos[0] for P in Psp])
    if len(qsp[qsp==123])==4:qine=1
    if len(qsp[qsp==-123])==4:qoute=1
    if (qine+qoute)==0:
        Fsp=[Planet,Planet,Planet,Planet]
        Qsp=array(\
            [pointInFigure(F,P)\
                 for F,P in zip(Fsp,Psp)])
        Psp=sortPolygonVertices(Psp[Qsp>-FIGTOL])
        Asp=convexPolygon(Psp)
    else:
        Asp=qine*Planet.area()

    #////////////////////////////////////////
    #TRANSIT AREA
    #////////////////////////////////////////
    return Asp,Psa,Feqs

def transitOblateAreaTime(t,S):
    #UPDATE POSITION
    updatePosition(S,t)
    #AREAS
    Es=transitOblateArea(S)
    At=Es[0]
    return At

def transitOblateAreaTimeFast(t,tcs,Ar,S):
    #ONLY COMPUTE AREA IF AT INGRESS OR EGRESS PHASE
    if t<=tcs[1]:At=0
    elif t<=tcs[2]:At=transitAreaOblateTime(t,S)
    elif t<=tcs[3]:At=Ar
    elif t<=tcs[4]:At=transitAreaOblateTime(t,S)
    else:At=0
    return At

def contactFunction(t,S,sgn):
    updatePosition(S,t)
    d=extremePointMultiple((S.Planet,
                            S.Ringext,
                            S.Ringint),
                           sgn=sgn)
    return d-1

def contactTime(S,Phase=EGRESS,Side=OUTSIDE,
                tola=1E-3,tolx=0.0,maxfun=10):
    """
    Calculate contact time for system S
    """
    tmin=S.tcen+(Phase-1)*S.dtstar
    tmax=S.tcen+(Phase+1)*S.dtstar
    t=brentq(contactFunction,tmin,tmax,args=(S,Side))
    return t

def contactTimes(S):
    """
    Calculate contact times for system S
    """
    #CALCULATE EXTREME TIMES
    t1=contactTime(S,Phase=INGRESS,Side=OUTSIDE)
    t4=contactTime(S,Phase=EGRESS,Side=OUTSIDE)

    #CHECK GRAZING TRANSIT
    if S.grazing:
        t2=t3=S.tcen
    else:
        t2=contactTime(S,Phase=INGRESS,Side=INSIDE)
        t3=contactTime(S,Phase=EGRESS,Side=INSIDE)
        
    return S.tcen,t1,t2,t3,t4
    
def transitAreaMontecarlo(Planet,Ringe,Ringi,NP=1E3):
    """
    Compute transit area via Montecarlo integration
    """
    Star=Figure(AR(0.0,0.0),1.0,1.0,1.0,0.0,'Star')

    mA1,dA1,xs1,ys1=montecarloArea([Star,Ringe,Planet,Ringi],
                                   [+1,+1,-1,-1],Npoints=NP)
    mA2,dA2,xs2,ys2=montecarloArea([Star,Planet],
                                   [+1,+1],Npoints=NP)
    mA=mA1+mA2
    dA=sqrt(dA1**2+dA2**2)
    if dA>0:
        nd=abs(log10(dA))+1
        fmA="%%.%df"%nd
        fmE="%%.%df"%(nd+1)
        mA=float(fmA%mA)
        dA=float(fmE%dA)
    return mA,dA,concatenate((xs1,xs2)),concatenate((ys1,ys2))

def transitAreaOblateMontecarlo(Planet,NP=1E3):
    """
    Compute transit area via Montecarlo integration
    """
    Star=Figure(AR(0.0,0.0),1.0,1.0,1.0,0.0,'Star')
    mA,dA,xs,ys=montecarloArea([Star,Planet],
                               [+1,+1],Npoints=NP)
    if dA>0:
        nd=abs(log10(dA))+1
        fmA="%%.%df"%nd
        fmE="%%.%df"%(nd+1)
        mA=float(fmA%mA)
        dA=float(fmE%dA)
    return mA,dA,xs,ys

#//////////////////////////////
#PLOT
#//////////////////////////////
def plotEllipse(ax,F,patch=False,**args):
    """
    Plot ellipse.  Optional arguments are for "plot" command.
    """
    C=F.C
    a=F.a
    b=F.b
    if patch:
        cir=pat.Circle((F.C[0],F.C[1]),F.a,**args)
        ax.add_patch(cir)
    else:
        Es=linspace(0,2*pi,1000)
        xs=a*cos(Es)
        ys=b*sin(Es)
        rs=array([rotTrans(AR(x,y),F.cost,F.sint,C) for x,y in zip(xs,ys)])
        ax.plot(rs[:,0],rs[:,1],'-',**args)
        #ax.plot([C[0]],[C[1]],'ko')
        
def plotPoint(ax,P,label=False,**args):
    """
    Plot a point.  Optional arguments are for "plot" command.
    """
    ax.plot([P.pos[0]],[P.pos[1]],'o',
            markeredgecolor='none',**args)
    if label:
        x=P.pos[0]+0.005
        y=P.pos[1]+0.005
        ax.text(x,y,P.name,fontsize=10)

###################################################
#NOTES FOR THIS RELEASE
###################################################
"""
Convention for interior or exterior:

- If intersection points between two curves does not exist we will
  assume that intersection points are 123 when one curve contains the
  other and -123 otherwise.
"""

###################################################
#ADDITIONAL
###################################################
###################################################
#NUMERICAL INTERSECTION
###################################################
def parsE(F):
    a=F.a;a2=a**2
    b=F.b;b2=b**2
    x=F.C[0];x2=x**2
    y=F.C[1];y2=y**2
    return x,x2,y,y2,a,a2,b,b2

def funsE(E):
    sE=sin(E)
    cE=cos(E)
    s2E=2*sE*cE
    c2E=cE**2-sE**2
    return sE,cE,s2E,c2E
    
def trigFunc(E,F):
    x,x2,y,y2,a,a2,b,b2=parsE(F)
    sE,cE,s2E,c2E=funsE(E)
    f=a2*cE**2+b2*sE**2+x2+y2+2*x*a*cE+2*y*b*sE-1
    return f

def trigFuncD(E,F):
    x,x2,y,y2,a,a2,b,b2=parsE(F)
    sE,cE,s2E,c2E=funsE(E)
    f=(b2-a2)*s2E-2*x*a*sE+2*y*b*cE
    return f

def trigFuncD2(E,F):
    x,x2,y,y2,a,a2,b,b2=parsE(F)
    sE,cE,s2E,c2E=funsE(E)
    f=2*(b2-a2)*c2E-2*x*a*cE-2*y*b*sE
    return f

def uniqueReals(xs,xtol=1E-3):
    es=array(xs)
    ilen=len(es)
    ns=[]
    while ilen>0:
        ies=array(range(ilen))
        ns+=[es[0]]
        dx=abs(es[0]-es)
        cond=dx<=xtol
        ires=ies[cond]
        es=delete(es,ires)
        ilen=len(es)
    return ns

def uniqueRoots(Es,F,xtol=1E-3):
    Er=sort(Es)
    fEr=trigFunc(Er,F)
    dEr=trigFuncD(Er,F)
    zEs=array(zip(Er,fEr,dEr))
    i=0
    clusters=[]

    while len(zEs)>0:
        ilen=len(zEs)
        ies=array(range(ilen))
        #if VERB:print "ies = ",ies
        zE=zEs[0]
        #if VERB:print "zE = ",zE*RAD
        dE=abs(zE[0]-zEs[:,0])
        #if VERB:print "zEs = ",zEs[:,0]*RAD
        #if VERB:print "dE = ",dE*RAD
        #if VERB:print "xtol = ",xtol*RAD
        cond=dE<=xtol
        ires=ies[cond]
        #if VERB:print "ires = ",ires
        cluster=zEs[ires]
        #if VERB:print "Cluster = ",cluster
        clusters+=[cluster]
        zEs=delete(zEs,ires,axis=0)
        #if VERB:print "zEs = ",zEs
        #if VERB:raw_input()

    roots=[]
    for cluster in clusters:
        #if VERB:print "*"*80
        isort=argsort(abs(cluster[:,2]))
        #if VERB:print "Sorted indexes: ",isort
        cluster=cluster[isort]
        #cluster[-1,2]*=-1
        #if VERB:print "Sorted cluster: ",cluster

        while len(cluster)>0:
            ilen=len(cluster)
            root=cluster[0,0]
            droot=cluster[0,2]
            roots+=[root]
            #if VERB:print 
            #if VERB:print "Size of cluster = ",ilen
            #if VERB:print "Main root = ",root
            #if VERB:print "Derivative = ",droot
            #if VERB:print "Cummulative roots = ",roots
            if ilen>1:
                ies=array(range(ilen))
                cond=droot*cluster[:,2]>0
                ires=ies[cond]
            else:
                #if VERB:print "Only one element."
                ires=ilen*[0]

            #if VERB:print "Common roots = ",ires
            if len(ires)==ilen:
                #if VERB:print "All values in cluster are the same."
                #if VERB:raw_input()
                break
            else:
                cluster=delete(cluster,ires,axis=0)
                #if VERB:print "Remaining cluster = ",cluster
            #if VERB:raw_input()

    #if VERB:print "="*80
    #if VERB:print "Final roots = ",array(roots)*RAD
    return roots

def ellipseDistanceDerivative(s,F,sgn):
    """
    Calculate the derivative of an ellipse point to origin
    """
    c=sgn*(1-s*s)**0.5
    f=(F.b*F.b-F.a*F.a)*s*c-F.a*F.C[0]*s+F.b*F.C[1]*c    
    return f

def extremePoint(Ellipse,sgn=-1):
    C=array(Ellipse.C)

    #IF ELLIPSE IS A CIRCLE
    if EQUAL(Ellipse.a,Ellipse.b):
        return MAG(C)+sgn*Ellipse.a

    #IF AN ELLIPSE
    Ellipse.C=abs(rotTrans(Ellipse.C,
                           +Ellipse.cost,
                           -Ellipse.sint,
                           AR(0,0)))
    sE=brentq(ellipseDistanceDerivative,0,sgn,
              args=(Ellipse,sgn),xtol=DISTANCETOL)
    cE=sgn*(1-sE**2)**0.5
    d=MAG(ellipsePointEcc(Ellipse,cE,sE))
    Ellipse.C=array(C)
    return d

def extremePoints(Ellipse):
    C=array(Ellipse.C)

    #IF ELLIPSE IS A CIRCLE
    if EQUAL(Ellipse.a,Ellipse.b):
        return MAG(C)-Ellipse.a,MAG(C)+Ellipse.a

    #IF AN ELLIPSE
    Ellipse.C=abs(rotTrans(Ellipse.C,
                           +Ellipse.cost,
                           -Ellipse.sint,
                           AR(0,0)))

    #CLOSEST
    sE=brentq(ellipseDistanceDerivative,0,CLOSEST,
              args=(Ellipse,CLOSEST),xtol=DISTANCETOL)
    cE=CLOSEST*(1-sE**2)**0.5
    dc=MAG(ellipsePointEcc(Ellipse,cE,sE))

    #FARTHEST
    sE=brentq(ellipseDistanceDerivative,0,FARTHEST,
              args=(Ellipse,FARTHEST),xtol=DISTANCETOL)
    cE=FARTHEST*(1-sE**2)**0.5
    df=MAG(ellipsePointEcc(Ellipse,cE,sE))

    Ellipse.C=array(C)
    return dc,df

def extremePointMultiple(Figures,sgn=-1):
    ds=[]
    for Figure in Figures:
        if EQUAL(Figure.b,ZERO):continue
        ds+=[extremePoint(Figure,sgn=sgn)]
    if sgn>0:d=max(ds)
    else:d=min(ds)
    return d

def extremePointsMultiple(Figures):
    ds=[]
    for Figure in Figures:
        if EQUAL(Figure.b,ZERO):continue
        ds+=[extremePoints(Figure)]
    dc=min(array(ds)[:,0])
    df=max(array(ds)[:,1])
    return dc,df

def eIc(EL,UC):
    """
    Computes the 4 intersection points between an ellipse
    (EL) and an unitary circle (UC).
    """
    #FIND THE FIRST GUESS
    xtol=INTERTOL
    E0=fsolve(trigFunc,0,args=(EL,),xtol=xtol)
    E0=mod(E0,2*pi)
    
    #CREATE AN ARRAY OF SEARCHING POINTS
    Es=linspace(E0,E0+2*pi,10)

    #SOLVE AROUND SEARCHING POINTS
    Eos=fsolve(trigFunc,Es,args=(EL,),xtol=xtol)
    Eos=mod(Eos,2*pi)

    #SELECT ONLY THOSE POINTS FULFILLING FUNCTION 
    fEos=trigFunc(Eos,EL)
    cond=abs(fEos)<INTERFUNTOL
    Eos=Eos[cond]

    #SELECT UNIQUE ROOTS
    Es=uniqueRoots(Eos,EL,xtol=ANGLETOL)
    #print "Ellipse %s, number of points: "%EL.name,len(Es)

    #CHECK IF AN EXCESS OF POINTS HAS BEEN CALCULATED
    xtol=ANGLETOL
    while len(Es)>4:
        #print "Correcting:",array(Es)*RAD
        Es=uniqueReals(Es,xtol=xtol)
        xtol+=1
        #print "Corrected:",array(Es)*RAD

    #CREATE LIST OF POINTS
    Ps=[]
    for E in Es:
        P=toPoint(ellipsePointEcc(EL,cos(E),sin(E)))
        P.fig1=EL;P.fig2=UC
        Ps+=[P]
        
    #RETURN POINTS
    if MAG(EL.C)<1:q=+1
    else:q=-1
    Ps+=[Point(AR(123*q,123*q),EL,UC)]*(4-len(Ps))

    return Ps

def eIcNumerical(EL,UC):
    """
    Computes the 4 intersection points between an ellipse
    (EL) and an unitary circle (UC).
    """
    #FIND THE FIRST GUESS
    xtol=INTERTOL
    E0=fsolve(trigFunc,0,args=(EL,),xtol=xtol)
    E0=mod(E0,2*pi)
    
    #CREATE AN ARRAY OF SEARCHING POINTS
    Es=linspace(E0,E0+2*pi,10)

    #SOLVE AROUND SEARCHING POINTS
    Eos=fsolve(trigFunc,Es,args=(EL,),xtol=xtol)
    Eos=mod(Eos,2*pi)

    #SELECT ONLY THOSE POINTS FULFILLING FUNCTION 
    fEos=trigFunc(Eos,EL)
    cond=abs(fEos)<INTERFUNTOL
    Eos=Eos[cond]

    #SELECT UNIQUE ROOTS
    Eos=uniqueRoots(Eos,EL,xtol=ANGLETOL)

    return array(Eos)

###################################################
#PHYSICS ROUTINES
###################################################
#//////////////////////////////
#KEPLER ORBITS
#//////////////////////////////
def timeOrbit(e,n,f):
    E=2*arctan(sqrt((1-e)/(1+e))*tan(f/2))
    M=E-e*sin(E)
    t=M/n
    return t

def eccentricAnomaly(e,M):
    """
    Using a very simple method based on: 
    Mikkola 1987
    """
    Ma=abs(M)
    c=4*e+0.5
    a=3*(1-e)/c
    b=-Ma/c
    y=sqrt(0.25*b*b+a*a*a/27.0)
    x=(-0.5*b+y)**(1./3)-(0.5*b+y)**(1./3)
    s=x-0.078*x**5/(1+e)
    E=Ma+e*(3*s-4*s**3)
    return sign(M)*E

def EoK(E,e=0.0,M=0.0):
    return E-e*sin(E)-M

def eccentricAnomalyNumerical(e,M):
    E=newton(EoK,M,args=(e,M),tol=1E-15)    
    return E

def eccentricAnomalyFast(e,M):
    """
    Mikkola, 1991
    Code at: http://smallsats.org/2013/04/20/keplers-equation-iterative-and-non-iterative-solver-comparison/
    """
    if e==0:return M
    a=(1-e)*3/(4*e+0.5);
    b=-M/(4*e+0.5);
    y=(b*b/4 +a*a*a/27)**0.5;
    x=(-0.5*b+y)**(1./3)-(0.5*b+y)**(1./3);
    w=x-0.078*x**5/(1 + e);
    E=M+e*(3*w-4*x**3);

    #NEWTON CORRECTION 1
    sE=sin(E)
    cE=cos(E)

    f=(E-e*sE-M);
    fd=1-e*cE;
    f2d=e*sE;
    f3d=-e*cE;
    f4d=e*sE;
    E=E-f/fd*(1+\
                  f*f2d/(2*fd*fd)+\
                  f*f*(3*f2d*f2d-fd*f3d)/(6*fd**4)+\
                  (10*fd*f2d*f3d-15*f2d**3-fd**2*f4d)*\
                  f**3/(24*fd**6))

    #NEWTON CORRECTION 2
    f=(E-e*sE-M);
    fd=1-e*cE;
    f2d=e*sE;
    f3d=-e*cE;
    f4d=e*sE;
    E=E-f/fd*(1+\
                  f*f2d/(2*fd*fd)+\
                  f*f*(3*f2d*f2d-fd*f3d)/(6*fd**4)+\
                  (10*fd*f2d*f3d-15*f2d**3-fd**2*f4d)*\
                  f**3/(24*fd**6))
    return E

def derivedSystemProperties(S):

    #STAR
    S.Flux=S.qeff*planckPhotons(500*NANO,700*NANO,S.Tstar)*S.Rstar**2/S.Dstar**2*(4*pi*S.Ddet**2)
    
    #ROTATION MATRICES
    S.Mi=rotMat([1,0,0],-S.iorb)
    S.Mw=rotMat([0,0,1],S.wp)
    S.Mos=dot(S.Mi,S.Mw)
    
    #RINGS
    S.Rp=S.Rplanet/S.Rstar
    S.Ri=S.fi*S.Rp
    S.Re=S.fe*S.Rp

    #ORBIT
    S.Porb=2*pi*sqrt(S.ap**3/(GCONST*(S.Mstar+S.Mplanet)))
    S.norb=2*pi/S.Porb

    S.fcen=270*DEG-S.wp
    S.Ecen=2*arctan(sqrt((1-S.ep)/(1+S.ep))*tan(S.fcen/2))
    S.Mcen=S.Ecen-S.ep*sin(S.Ecen)
    
    S.tcen=S.Mcen/S.norb
    S.rcen=ellipseRadiusE(S.ap,S.ep,S.fcen)
    
    rpcen=AR3(S.rcen*cos(S.fcen),S.rcen*sin(S.fcen),0)
    S.Pcen=dot(S.Mos,rpcen)
    S.Borb=S.Pcen[1]/S.Rstar
    
    S.orbit=Orbit(S.ap/S.Rstar,S.ep,S.Porb,S.Mos)
    
    if abs(S.Borb)>1:
        print "This configuration does not lead to a transit."
        exit(1)
        
    S.C=AR(0,0)
    S.Star=Figure(AR(0,0),
                  1.0,1.0,
                  1.0,0.0,
                  'Star')
    S.Planet=Figure(S.C,
                    S.Rp,S.Rp*(1-S.fp),
                    1.0,0.0,
                    'Planet')

def updatePlanetRings(S,phir,ir):

    #//////////////////////////////////////////////////
    #ROTATION MATRIX FROM EQUATORIAL SYSTEM TO SKY
    #//////////////////////////////////////////////////
    Mpr=rotMat([0,0,1],phir)
    Mir=rotMat([1,0,0],-ir)
    Mro=dot(Mpr,Mir)
    S.Mrs=dot(S.Mi,Mro)

    #//////////////////////////////////////////////////
    #EFFECTIVE ROTATION AND ROLL
    #//////////////////////////////////////////////////
    rx=dot(S.Mrs,[1.0,0.0,0.0])
    ry=dot(S.Mrs,[0.0,1.0,0.0])
    rz=dot(S.Mrs,[0.0,0.0,1.0])
    S.ieff=arccos(abs(dot(rz,[0,0,1])))
    S.teff=-sign(rz[0])*ARCTAN(abs(rz[0]),abs(rz[1]))

    #//////////////////////////////////////////////////
    #PLANET SHAPE
    #//////////////////////////////////////////////////
    #See theorem at: 
    #http://twisee.com/all/484872495980625920
    a=S.Planet.a
    b=S.Planet.b
    S.Planet.b=b*(sin(S.ieff)**2+\
                      (a/b)**2*cos(S.ieff)**2)**0.5
    S.Planet.cost=cos(S.teff)
    S.Planet.sint=sin(S.teff)
    
    #//////////////////////////////////////////////////
    #RING BLOCK FACTOR
    #//////////////////////////////////////////////////
    if EQUAL(S.ieff,90.0):S.block=1.0
    elif EQUAL(S.tau,ZERO):S.block=0.0
    else:S.block=1-exp(-S.tau/cos(S.ieff))

    #//////////////////////////////////////////////////
    #ESTIMATED CONTACT TIMES
    #//////////////////////////////////////////////////
    df=S.Rstar*sqrt(1-S.Borb**2)/S.rcen
    f15=S.fcen-df
    r15=ellipseRadiusE(S.ap,S.ep,f15) 
    P15=dot(S.Mos,AR3(r15*cos(f15),r15*sin(f15),0))/S.Rstar
    t15=timeOrbit(S.ep,S.norb,f15)
    f35=S.fcen+df
    r35=ellipseRadiusE(S.ap,S.ep,f35) 
    t35=timeOrbit(S.ep,S.norb,f35)
    P35=dot(S.Mos,AR3(r35*cos(f35),r35*sin(f35),0))/S.Rstar
    if t15*t35<0 and t35<0:t35+=Porb

    #ESTIMATED TRANSIT TIME
    S.dtrans=t35-t15

    #TIME TO COVER THE PLANETARY RADIUS
    df=S.Rplanet/S.rcen
    f=S.fcen-df
    r=ellipseRadiusE(S.ap,S.ep,f) 
    P=dot(S.Mos,AR3(r*cos(f),r*sin(f),0))/S.Rstar
    t=timeOrbit(S.ep,S.norb,f)
    S.dtplanet=S.tcen-t

    #TIME TO COVER THE RADIUS OF THE STAR
    df=S.Rstar/S.rcen
    f=S.fcen-df
    r=ellipseRadiusE(S.ap,S.ep,f) 
    P=dot(S.Mos,AR3(r*cos(f),r*sin(f),0))/S.Rstar
    t=timeOrbit(S.ep,S.norb,f)
    S.dtstar=S.tcen-t
    
    #CREATE RING OBJECTS
    S.Ringext=Figure(S.C,
                     S.Re,S.Re*cos(S.ieff),
                     cos(S.teff),sin(S.teff),
                     'Ringext')
    S.Ringint=Figure(S.C,
                     S.Ri,S.Ri*cos(S.ieff),
                     cos(S.teff),sin(S.teff),
                     'Ringint')

    #CHECK IF GRAZING
    updatePosition(S,S.tcen)
    df=extremePointMultiple((S.Planet,
                             S.Ringext,
                             S.Ringint),
                            sgn=FARTHEST)

    if df>1:
        print "Grazing configuration."
        S.grazing=1
    else:
        S.grazing=0
    
def systemShow(S):
    pass
    print "Star primary:"
    print TAB,"Ms = %e kg"%S.Mstar
    print TAB,"Rs = %e kg"%S.Rstar
    print "Planet primary:"
    print TAB,"Mp = %e kg = %e Mstar"%(S.Mplanet,
                                       S.Mplanet/S.Mstar)
    print TAB,"Rp = %e kg = %e Rstar"%(S.Rplanet,
                                       S.Rplanet/S.Rstar)
    print "Rings primary:"
    print TAB,"fi,fe = %e,%e Rp"%(S.fi,S.fe)
    print TAB,"Inclination (orbit) = %.1f deg"%(S.ir*RAD)
    print TAB,"Roll (orbit) = %.1f deg"%(S.phir*RAD)
    print TAB,"Opacity = %.2f"%(S.tau)
    print "Orbit primary:"
    print TAB,"ap = %e km = %e AU = %e Rstar"%(S.ap,
                                               S.ap/AU,
                                               S.ap/S.Rstar)
    print TAB,"Eccentricity = %.2f"%(S.ep)
    print TAB,"Inclination (visual) = %.2f deg"%(S.iorb*RAD)
    print TAB,"Periapsis argument = %.2f deg"%(S.wp*RAD)
    print
    print "Star derivative:"
    print "Planetary derivative:"
    print TAB,"Radius (relative) = %e Rstar"%(S.Rp)
    print "Rings derivative:"
    print TAB,"Internal ring (relative) = %.2f Rstar"%(S.Ri)
    print TAB,"External ring (relative) = %.2f Rstar"%(S.Re)
    print TAB,"Apparent inclination = %.2f deg"%(S.ieff*RAD)
    print TAB,"Apparent roll = %.2f deg"%(S.teff*RAD)
    print "Orbit derivative:"
    print TAB,"Period = %e s = %e h = %e d = %e yr"%(S.Porb,S.Porb/HOUR,S.Porb/DAY,S.Porb/YEAR)
    print TAB,"Mean Angular velocity = %e rad/s = %e Rstar/s = %e Rp/s"%(S.norb,S.norb*S.rcen/S.Rstar,S.norb*S.rcen/(S.Rp*S.Rstar))
    print TAB,"Central true anomaly = %e deg"%(S.fcen*RAD)
    print TAB,"Central eccentric anomaly = %e deg"%(S.Ecen*RAD)
    print TAB,"Central mean anomaly = %e deg"%(S.Mcen*RAD)
    print TAB,"Central radius = %e km = %e AU = %e Rstar"%(S.rcen,S.rcen/AU,S.rcen/S.Rstar)
    print TAB,"Impact parameter = %e Rstar"%(S.Borb)
    print TAB,"Central time = %e s = %e Porb"%(S.tcen,S.tcen/S.Porb)
    print TAB,"Estimated transit duration = %e s = %e h"%(S.dtrans,S.dtrans/3600.0)

def updatePosition(S,t):
    S.t=t
    S.Mt=S.orbit.n*t
    S.Et=eccentricAnomalyFast(S.orbit.e,S.Mt)
    S.rp=[S.orbit.a*cos(S.Et)-S.orbit.a*S.orbit.e,
          S.orbit.b*sin(S.Et),0.0]
    S.rs=dot(S.orbit.Mos,S.rp)
    S.C=AR(S.rs[0],S.rs[1])
    S.Planet.C=S.Ringext.C=S.Ringint.C=S.C

def planckPhotonDistrib(lamb,T):
    B=2*HP*CSPEED**2/(lamb**5)/\
        (exp(HP*CSPEED/(KB*T*lamb))-1)
    J=pi*B/(HP*CSPEED/lamb)
    return J

def planckPhotons(lamb1,lamb2,T):
    N,dN=integrate(planckPhotonDistrib,
                   lamb1,lamb2,args=(T,))
    return N

def gaussianQuadrature(func,a,b,args=(0,)):
    x1=0
    x2=+(3./5)**0.5
    x3=-x2
    w1=8./9
    w2=w3=5./9
    bam2=(b-a)/2
    bap2=(b+a)/2
    try:
        integral=bam2*(w1*func(bam2*x1+bap2,*args)+
                       w2*func(bam2*x2+bap2,*args)+
                       w3*func(bam2*x3+bap2,*args))
    except:
        integral=bam2*(w1*func(bam2*x1+bap2,*args)+
                       w2*func(bam2*x2+bap2,*args)+
                       w3*func(bam2*x3+bap2,*args))
    return integral

def onlyPlanet(S):
    P=copyObject(S)
    P.Ringext.b=P.Ringint.b=0.0
    return P
