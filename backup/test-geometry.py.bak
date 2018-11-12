from geotrans import *

###################################################
#COMMON
###################################################
NP=1E4
Rp=0.0836255747894
Rp=0.085
fp=0.0
C=AR(0.9,0.4)
#C=AR(0.0,0.9995)
#C=AR(8.25560060e-01,-4.77023573e-14)

ieff=75*DEG
teff=40.0*DEG
fi=2.0
fe=2.5

Ri=fi*Rp
Re=fe*Rp

#PLANET OBLATENESS
a=Rp
b=Rp*(1-fp)
#PROJECTED
b=b*(sin(ieff)**2+\
         (a/b)**2*cos(ieff)**2)**0.5

UnitaryCircle=Figure(AR(0,0),
                     1.0,1.0,
                     1.0,0.0,
                     'BigCircle')
SmallEllipse=Figure(C,
                    Re,Re*cos(ieff),
                    cos(teff),sin(teff),
                    'SmallEllipse')
Star=Figure(AR(0,0),
            1.0,1.0,
            1.0,0.0,
            'Star')
Planet=Figure(C,
              a,b,
              cos(teff),sin(teff),
              'Planet')
Ringext=Figure(C,
               Re,Re*cos(ieff),
               cos(teff),sin(teff),
               'Ringext')
Ringint=Figure(C,
               Ri,Ri*cos(ieff),
               cos(teff),sin(teff),
               'Ringint')

###################################################
#TEST ROUTINES
###################################################
def test_area_striping(Planet,Ringext,Ringint):
    print BARL,"Testing Area Striping",RBAR

    #////////////////////////////////////////
    #SYSTEM
    #////////////////////////////////////////
    S=dict2obj(dict(Planet=Planet,
                    Ringext=Ringext,
                    Ringint=Ringint))
    c1=0.7;c2=-0.24
    
    #////////////////////////////////////////
    #EXTREMES
    #////////////////////////////////////////
    dc,df=extremePointsMultiple((Planet,Ringext))
    Sc=Figure(AR(0,0),dc,dc,1.0,0.0,'Strip low')
    Sf=Figure(AR(0,0),df,df,1.0,0.0,'Strip low')
    print "Extremes: ",dc,df

    #////////////////////////////////////////
    #STRIPING AREA
    #////////////////////////////////////////
    deltad=(df-dc)/5.0
    ds=secureArange(dc,df,deltad)
    print "Sampling points: ",ds
    Ats=areaStriping(S,ds)
    print "Areas = ",Ats
    print "Total Area via striping = ",Ats.sum()
    At=transitArea(S)[0]
    print "Total Area analytic = ",At

    #////////////////////////////////////////
    #STRIPING VIA LIMB DARKENING
    #////////////////////////////////////////
    deltad=(df-dc)/10.0
    ds=secureArange(dc,min(1,df),deltad)
    ds=limbStriping(ds,c1,c2)
    print "Limiting radii = ",ds
    dss=0.5*(ds[:-1]+ds[1:])
    print "Sampling = ",dss
    iss=limbDarkeningNormalized(dss,c1,c2)
    print "Limb darkening = ",iss
    Ats=areaStriping(S,ds)
    print "Areas = ",Ats
    its=Ats*iss
    fts=its.sum()
    print "Total Area via ld. striping = ",Ats.sum()
    print "Total flux via ld. striping = ",its.sum()

    #////////////////////////////////////////
    #PRECISE STRIPING AND LIMB DARKENING
    #////////////////////////////////////////
    dsf=[];ftsf=0
    """
    deltad=(df-dc)/30.0
    dsf=secureArange(dc,min(1,df),deltad)
    dss=0.5*(dsf[:-1]+dsf[1:])
    print "Fine Sampling = ",dss
    iss=limbDarkeningNormalized(dss,c1,c2)
    print "Fine Limb darkening = ",iss
    Ats=areaStriping(S,dsf)
    print "Areas = ",Ats
    its=Ats*iss
    print "Total Area via fine striping = ",Ats.sum()
    ftsf=its.sum()
    print "Total flux via fine striping = ",ftsf
    #"""

    #////////////////////////////////////////
    #PLOT
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    #PLOT EXTREMES
    #plotEllipse(ax,Sc,color='r',linestyle='--')
    #plotEllipse(ax,Sf,color='b',linestyle='--')

    for d in ds[:]:
        E=Figure(AR(0,0),d,d,1.0,0.0,'Strip')
        plotEllipse(ax,E,color='k',linestyle='--')
    for d in dsf[:]:
        E=Figure(AR(0,0),d,d,1.0,0.0,'Strip')
        plotEllipse(ax,E,color='b',linestyle=':')

    #PLOT FIGURES
    plotEllipse(ax,Star,color='y')
    plotEllipse(ax,Planet,color='b')
    plotEllipse(ax,Ringext,color='k')
    plotEllipse(ax,Ringint,color='r')

    rng=1.5
    xmin=Planet.C[0]-rng*Re;
    xmax=Planet.C[0]+rng*Re
    ymin=Planet.C[1]-rng*Re;
    ymax=Planet.C[1]+rng*Re
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("plots/test_area_striping.png")
    return fts,ftsf

def test_area_strip():
    print BARL,"Testing Strip Area",RBAR

    #////////////////////////////////////////
    #SYSTEM
    #////////////////////////////////////////
    S=dict2obj(dict(Planet=Planet,
                    Ringext=Ringext,
                    Ringint=Ringint))
    
    #////////////////////////////////////////
    #EXTREMES
    #////////////////////////////////////////
    dc,df=extremePointsMultiple((Planet,Ringext))
    Sc=Figure(AR(0,0),dc,dc,1.0,0.0,'Strip low')
    Sf=Figure(AR(0,0),df,df,1.0,0.0,'Strip low')
    print "Extremes: ",dc,df

    Rlow=0.7
    Rup=0.8
    Slow=Figure(AR(0,0),Rlow,Rlow,1.0,0.0,'Strip low')
    Sup=Figure(AR(0,0),Rup,Rup,1.0,0.0,'Strip low')
    print "Strip limits: ",Rlow,Rup

    #////////////////////////////////////////
    #STRIP AREA
    #////////////////////////////////////////
    At=stripArea(S,Rlow,Rup)
    print "Analytic Area = %e"%At

    #////////////////////////////////////////
    #MONTECARLO AREA
    #////////////////////////////////////////
    #"""
    mA1,dA1,xs1,ys1=montecarloArea([Sup,Slow,Ringext,Planet,Ringint],
                                   [+1,-1,+1,-1,-1],Npoints=NP,excl=2)
    #mA1,dA1,xs1,ys1=0,0,[],[]
    mA2,dA2,xs2,ys2=montecarloArea([Sup,Slow,Planet],
                                   [+1,-1,+1],Npoints=NP,excl=2)
    #mA2,dA2,xs2,ys2=0,0,[],[]
    mA=mA1+mA2
    dA=sqrt(dA1**2+dA2**2)
    xs,ys=concatenate((xs1,xs2)),concatenate((ys1,ys2))
    print "Montecarlo Area = %e +/- %e"%(mA,dA)
    if dA>0:ns=abs(At-mA)/dA
    else:ns=0
    print "Discrepance = %f sigmas"%(ns)
    #"""
    #xs,ys=[],[]

    #////////////////////////////////////////
    #PLOT
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    #PLOT EXTREMES
    plotEllipse(ax,Sc,color='r',linestyle='--')
    plotEllipse(ax,Sf,color='b',linestyle='--')

    #PLOT STRIP
    plotEllipse(ax,Slow,color='k',linewidth=3)
    plotEllipse(ax,Sup,color='k',linewidth=3)

    #PLOT MONTECARLO POINTS
    ax.plot(xs,ys,'r.',markersize=2)

    #PLOT FIGURES
    plotEllipse(ax,Star,color='y')
    plotEllipse(ax,Planet,color='b')
    plotEllipse(ax,Ringext,color='k')
    plotEllipse(ax,Ringint,color='r')

    rng=1.5
    xmin=Planet.C[0]-rng*Re;
    xmax=Planet.C[0]+rng*Re
    ymin=Planet.C[1]-rng*Re;
    ymax=Planet.C[1]+rng*Re
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("plots/test_area_strip.png")

def test_area_oblate():
    print BARL,"Testing Area Oblate Planet",RBAR
    #////////////////////////////////////////
    #FIND EXTREMES
    #////////////////////////////////////////
    dc,df=extremePoints((Planet))

    #////////////////////////////////////////
    #CALCULATE AREAS
    #////////////////////////////////////////
    S=dict2obj(dict(Planet=Planet,
                    Ringext=Ringext,
                    Ringint=Ringint))
    Es=transitAreaOblate(S)

    At=Es[0]
    Fs=Es[-1]
    Ps=Es[-2]
    print "Transit Area: ",At

    #////////////////////////////////////////
    #MONTECARLO AREAS
    #////////////////////////////////////////
    Am,dA,xs,ys=transitAreaOblateMontecarlo(Planet,NP=NP)
    print "Montecarlo Area = %f +/- %f"%(Am,dA)
    if dA>0:ns=abs(At-Am)/dA
    else:ns=0
    print "Discrepance = %f sigmas"%(ns)

    #////////////////////////////////////////
    #PLOT
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    #PLOT EXTREMES
    plotEllipse(ax,Figure(AR(0,0),dc,dc,1,0,''),
                color='r',linestyle='--')
    plotEllipse(ax,Figure(AR(0,0),df,df,1,0,''),
                color='b',linestyle='--')
    ax.plot(xs,ys,'r.',markersize=2)
    
    #PLOT FIGURES
    plotEllipse(ax,Star,color='y')
    plotEllipse(ax,Planet,color='b')

    rng=1.5
    xmin=Planet.C[0]-rng*Re;
    xmax=Planet.C[0]+rng*Re
    ymin=Planet.C[1]-rng*Re;
    ymax=Planet.C[1]+rng*Re
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("plots/test_area_oblate.png")

    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    #PLOT POINTS
    for P in Ps:
        if 'sr' in P.name or True:
            plotPoint(ax,P,label=True)

    #PLOT EXTREMES
    plotEllipse(ax,Figure(AR(0,0),dc,dc,1,0,''),
                color='r',linestyle='--')
    plotEllipse(ax,Figure(AR(0,0),df,df,1,0,''),
                color='b',linestyle='--')

    #PLOT FIGURES
    plotEllipse(ax,Star,color='y')
    plotEllipse(ax,Fs[0],color='b')
    ax.plot(xs,ys,'r.',markersize=2)

    rng=1.5
    xmin=Fs[0].C[0]-rng*Re;
    xmax=Fs[0].C[0]+rng*Re
    ymin=Fs[0].C[1]-rng*Re;
    ymax=Fs[0].C[1]+rng*Re
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("plots/test_area_oblate_eqs.png")

def test_areas():
    print BARL,"Testing Areas",RBAR
    #////////////////////////////////////////
    #FIND EXTREMES
    #////////////////////////////////////////
    dc,df=extremePointsMultiple((Planet,Ringext))
    
    #////////////////////////////////////////
    #CALCULATE AREAS
    #////////////////////////////////////////
    S=dict2obj(dict(Planet=Planet,
                    Ringext=Ringext,
                    Ringint=Ringint))

    ae=Ringext.a
    be=Ringext.b
    ai=Ringint.a
    bi=Ringint.b
    r=Planet.a
    Atap=pi*r**2
    f=1.0
    Atae=pi*ae*be-4*ae*be/(ae**2-be**2)*sqrt((r**2-be**2)*(ae**2-r**2))-2*(r-ae*sqrt((r**2-be**2)/(ae**2-be**2)))*be*sqrt((ae**2-r**2)/(ae**2-be**2))-f*(be-be*sqrt((ae**2-r**2)/(ae**2-be**2)))*ae*sqrt((r**2-be**2)/(ae**2-be**2))
    Atai=pi*ai*bi-4*ai*bi/(ai**2-bi**2)*sqrt((r**2-bi**2)*(ai**2-r**2))-2*(r-ai*sqrt((r**2-bi**2)/(ai**2-bi**2)))*bi*sqrt((ai**2-r**2)/(ai**2-bi**2))-f*(bi-bi*sqrt((ai**2-r**2)/(ai**2-bi**2)))*ai*sqrt((r**2-bi**2)/(ai**2-bi**2))
    Ata=Atap+Atae-Atai
    print "Analytic area: ",Ata


    Es=transitArea(S)

    At=Es[0]
    Fs=Es[-1]
    Ps=Es[-2]
    print "Transit Area: ",At
    print "Discrepance with analytic (%): ",abs(At-Ata)/At*100
    print TAB,"Planet Area: ",Es[1]
    print TAB,"External Ring Area: ",Es[2]
    print TAB,"Internal Ring Area: ",Es[3]
    print TAB,"Common External Ring Area: ",Es[4]
    print TAB,"Common Internal Ring Area: ",Es[5]

    #////////////////////////////////////////
    #MONTECARLO AREAS
    #////////////////////////////////////////
    Am,dA,xs,ys=transitAreaMontecarlo(Planet,Ringext,Ringint,NP=NP)
    print "Montecarlo Area = %f +/- %f"%(Am,dA)
    #PARTIAL
    A,x,x,x=montecarloArea([Star,Planet],
                           [+1,+1],Npoints=1E3)
    print TAB,"Montecarlo Planet = ",A
    A,x,x,x=montecarloArea([Star,Ringext],
                           [+1,+1],Npoints=1E3)
    print TAB,"Montecarlo External Ring = ",A
    A,x,x,x=montecarloArea([Star,Planet,Ringext],
                           [+1,+1,+1],Npoints=1E3)
    print TAB,"Montecarlo External Ring Common = ",A
    A,x,x,x=montecarloArea([Star,Planet,Ringint],
                           [+1,+1,+1],Npoints=1E3)
    print TAB,"Montecarlo Internal Ring Common = ",A

    if dA>0:ns=abs(At-Am)/dA
    else:ns=0
    print "Discrepance = %f sigmas"%(ns)

    Am,dA,xsn,ysn=transitAreaMontecarlo(Fs[0],Fs[1],Fs[2],NP=NP)

    #////////////////////////////////////////
    #PLOT
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    #PLOT EXTREMES
    plotEllipse(ax,Figure(AR(0,0),dc,dc,1,0,''),
                color='r',linestyle='--')
    plotEllipse(ax,Figure(AR(0,0),df,df,1,0,''),
                color='b',linestyle='--')
    ax.plot(xs,ys,'r.',markersize=2)
    
    #PLOT FIGURES
    plotEllipse(ax,Star,color='y')
    plotEllipse(ax,Planet,color='b')
    plotEllipse(ax,Ringext,color='k')
    plotEllipse(ax,Ringint,color='r')

    rng=1.5
    xmin=SmallEllipse.C[0]-rng*Re;
    xmax=SmallEllipse.C[0]+rng*Re
    ymin=SmallEllipse.C[1]-rng*Re;
    ymax=SmallEllipse.C[1]+rng*Re
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("plots/test_areas.png")


    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    #PLOT POINTS
    for P in Ps:
        if 'sr' in P.name or True:
            plotPoint(ax,P,label=True)

    #PLOT EXTREMES
    plotEllipse(ax,Figure(AR(0,0),dc,dc,1,0,''),
                color='r',linestyle='--')
    plotEllipse(ax,Figure(AR(0,0),df,df,1,0,''),
                color='b',linestyle='--')

    #PLOT FIGURES
    plotEllipse(ax,Star,color='y')
    plotEllipse(ax,Fs[0],color='b')
    plotEllipse(ax,Fs[1],color='k')
    plotEllipse(ax,Fs[2],color='r')
    ax.plot(xsn,ysn,'r.',markersize=2)

    rng=1.5
    xmin=Fs[0].C[0]-rng*Re;
    xmax=Fs[0].C[0]+rng*Re
    ymin=Fs[0].C[1]-rng*Re;
    ymax=Fs[0].C[1]+rng*Re
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("plots/test_areas_eqs.png")

def test_analytic_area():
    print BARL,"Testing Approximate Formula for Area",RBAR
    #////////////////////////////////////////
    #ERROR
    #////////////////////////////////////////
    imin=arccos(1/fe)
    print "Minimum inclination angle: ",imin*RAD
    S=dict2obj(dict(Planet=Planet,
                    Ringext=Ringext,
                    Ringint=Ringint))

    #////////////////////////////////////////
    #CALCULATE AREAS
    #////////////////////////////////////////
    def factor(si,ci,f,flag=1,w1=1.0,w2=1.0):
        f=1-flag*2/pi*(w1/f**2+w2/(f**2*si**2)*\
                      sqrt(1-f**2*ci**2)*sqrt(f**2-1))/((w1+w2)/2.0)
        return f
    def errorArea(ieff,Verbose=False):
        #ANALYTIC APPROXIMATE
        Z=AR(0,0)
        si=sin(ieff)
        ci=cos(ieff)
        S.Planet.C=S.Ringext.C=S.Ringint.C=Z
        S.Ringext.b=Re*ci
        S.Ringint.b=Ri*ci

        r=S.Planet.a
        ae=S.Ringext.a
        be=S.Ringext.b
        de=be*sqrt(1-r**2/ae**2)/sin(ieff)

        ai=S.Ringint.a
        bi=S.Ringint.b
        di=bi*sqrt(1-r**2/ai**2)/sin(ieff)
        
        """
        Ae=ae*be*arcsin(de/be)-r**2*arcsin(de/r)
        Ai=ai*bi*arcsin(di/bi)-r**2*arcsin(di/r)
        Ar=Ae-Ai
        """

        if cos(ieff)<1/fe:
            Fe=fe**2*cos(ieff)*2/pi*arcsin(1/(fe*sin(ieff))*sqrt(fe**2-1))-\
                2/pi*arcsin(cos(ieff)/sin(ieff)*sqrt(fe**2-1))
        else:Fe=fe**2*cos(ieff)-3/2
        if cos(ieff)<1/fi:
            Fi=fi**2*cos(ieff)*2/pi*arcsin(1/(fi*sin(ieff))*sqrt(fi**2-1))-\
                2/pi*arcsin(cos(ieff)/sin(ieff)*sqrt(fi**2-1))
        else:Fi=fi**2*cos(ieff)-3/2

        Aa=pi*Rp**2*(1+Fe-Fi)

        """
        w1=2
        w2=1
        Ae=pi*fe**2*Rp**2*cos(ieff)*factor(si,ci,fe,flag=1,w1=w1,w2=w2)
        Ai=pi*fi**2*Rp**2*cos(ieff)*factor(si,ci,fi,flag=1,w1=w1,w2=w2)
        Ar=Ae-Ai
        Aa=pi*Rp**2+Ar
        """

        #ANALYTIC EXACT
        Es=transitArea(S)
        At=Es[0]
        Fs=Es[-1]
        Ps=Es[-2]

        #ERROR
        dA=abs(At-Aa)/At*100
        
        if Verbose:
            print "Summary Exact:"
            print TAB,"Planet Area: ",Es[1]
            print TAB,"External Ring Area: ",Es[2]-Es[5]
            print TAB,"Internal Ring Area: ",Es[3]-Es[5]
            print "Transit Area: ",At

            print "Summary Approximate:"
            print TAB,"Planet Area: ",pi*Rp**2
            print TAB,"External Ring Area: ",Ae
            print TAB,"Internal Ring Area: ",Ai
            print "Transit Area: ",At

            print "Error area (%): ",dA
        return At,Aa,dA

    #////////////////////////////////////////
    #CALCULATE AREAS
    #////////////////////////////////////////
    ieffs=linspace(1*DEG,89*DEG,100)
    dAs=[errorArea(ieff)[2] for ieff in ieffs]

    #////////////////////////////////////////
    #PLOT
    #////////////////////////////////////////
    fig=plt.figure(figsize=(6,8))
    ax=fig.gca()
    ax.plot(ieffs*RAD,dAs,'-')
    ax.set_xlabel(r"$i_{\rm eff}$")
    ax.set_xlabel(r"Area Relative Error (%)")
    ax.grid()
    fig.savefig("plots/test_analytic_area_error.png")

def test_points():
    print BARL,"Testing Intersection Points",RBAR
    #////////////////////////////////////////
    #FIND EXTREMES
    #////////////////////////////////////////
    dc,df=extremePointsMultiple((Planet,Ringext))

    #////////////////////////////////////////
    #FIND INTERSECTION POINTS
    #////////////////////////////////////////
    #PLANET AND STAR
    Ps=[]
    Psp1,Psp2=cIc(Planet,Star)
    Psp1.name='Psp1';Psp2.name='Psp2'
    Ps+=[Psp1,Psp2]

    #PLANET AND RINGEXT
    Ppre1,Ppre2,Ppre3,Ppre4=cIe(Ringext,Planet)
    Ppre1.name='Ppre1';Ppre2.name='Ppre2';
    Ppre3.name='Ppre3';Ppre4.name='Ppre4';
    Ps+=[Ppre1,Ppre2,Ppre3,Ppre4]
    
    #PLANET AND RINGINT
    Ppri1,Ppri2,Ppri3,Ppri4=cIe(Ringint,Planet)
    Ppri1.name='Ppri1';Ppri2.name='Ppri2';
    Ppri3.name='Ppri3';Ppri4.name='Ppri4';
    Ps+=[Ppri1,Ppri2,Ppri3,Ppri4]

    #STAR AND RINGEXT
    Psre1,Psre2,Psre3,Psre4=eIc(Ringext,Star)
    Psre1.name='Psre1';Psre2.name='Psre2';
    Psre3.name='Psre3';Psre4.name='Psre4';
    Ps+=[Psre1,Psre2,Psre3,Psre4]
    
    #STAR AND RINGINT
    Psri1,Psri2,Psri3,Psri4=eIc(Ringint,Star)
    Psri1.name='Psri1';Psri2.name='Psri2';
    Psri3.name='Psri3';Psri4.name='Psri4';
    Ps+=[Psri1,Psri2,Psri3,Psri4]

    #////////////////////////////////////////
    #PLOT ELLIPSES
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    #PLOT POINTS
    for P in Ps:
        if 's' in P.name or False:
            plotPoint(ax,P,label=True)

    #PLOT EXTREMES
    plotEllipse(ax,Figure(AR(0,0),dc,dc,1,0,''),
                color='r',linestyle='--')
    plotEllipse(ax,Figure(AR(0,0),df,df,1,0,''),
                color='b',linestyle='--')

    #PLOT FIGURES
    plotEllipse(ax,Star,color='y')
    plotEllipse(ax,Planet,color='b')
    plotEllipse(ax,Ringext,color='k')
    plotEllipse(ax,Ringint,color='r')

    rng=1.5
    xmin=SmallEllipse.C[0]-rng*Re;
    xmax=SmallEllipse.C[0]+rng*Re
    ymin=SmallEllipse.C[1]-rng*Re;
    ymax=SmallEllipse.C[1]+rng*Re
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("plots/test_points.png")
    
def test_inclusion():
    print BARL,"Testing Extremes",RBAR
    #////////////////////////////////////////
    #FIGURE
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    #////////////////////////////////////////
    #TEST DISTANCE
    #////////////////////////////////////////
    dc,sE,cE=extremePoint2(SmallEllipse,sgn=CLOSEST)

    print sE,cE

    a=SmallEllipse.a
    b=SmallEllipse.b
    x=C[0]
    y=C[1]
    c2=a**2-b**2

    #NO INCLINATION
    u=b*y/(c2-a*x)
    print "Approximation no inclination: ",u
    x=sqrt(1-y**2*(1-b**2/(a)))

    #x=1-0.5*y**2*(1-b**2/a)
    #dc,df=extremePoints(SmallEllipse)
    #x=sqrt(1-y**2*(1-b**2/(a*dc))**2)
    #dc,df=extremePoints(SmallEllipse)
    #x=sqrt(1-y**2*(1-b**2/(a*dc))**2)
    #x=sqrt(1-y**2*(1-b**2/(a*x)))
    #x=sqrt(1-y**2*(1-b**2/(a*x)))
    
    print dc,sE,cE,u,x

    SmallEllipse.C[0]=x+a

    #WITH INCLINATION
    #"""
    x=C[0]
    y=C[1]
    #Approx 0: E<<1
    sEt0=b*(y*cos(teff)-x*sin(teff))/\
        (c2-a*(x*cos(teff)+y*sin(teff)))
    #Approx 1: x=1
    sEt1=b*(y*cos(teff)-sin(teff))/\
        (c2-a*(cos(teff)+y*sin(teff)))
    #Approx 2: x=1
    sEt1=b*(y*cos(teff)-sin(teff))/\
        (c2-a*(cos(teff)+y*sin(teff)))
    print "Approximation inclination 0: ",sEt0
    print "Approximation inclination 1: ",sEt1

    x=sqrt(1-(y-a*sin(teff)-b**2*(y*cos(teff)+sin(teff)/a))**2)
    x=sqrt(1-a**2*(sin(teff)-y/a)**2*(1-b**2/a))
    SmallEllipse.C[0]=x+a*cos(teff)
    #"""

    dc,df=extremePoints(SmallEllipse)
    print "Extremes: ",dc,df,abs(dc-1)*100,"%"

    #////////////////////////////////////////
    #PLOT ELLIPSES
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()
    
    plotEllipse(ax,Figure(AR(0,0),dc,dc,1,0,''),
                color='r',linestyle='--')
    plotEllipse(ax,Figure(AR(0,0),df,df,1,0,''),
                color='b',linestyle='--')
    plotEllipse(ax,UnitaryCircle)
    plotEllipse(ax,SmallEllipse)
    
    rng=1.5
    xmin=SmallEllipse.C[0]-rng*Re;
    xmax=SmallEllipse.C[0]+rng*Re
    ymin=SmallEllipse.C[1]-rng*Re;
    ymax=SmallEllipse.C[1]+rng*Re
    #xmin=-1.5;xmax=1.5
    #ymin=-1.5;ymax=1.5
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("plots/test_inclusion.png")

def test_analytic_contact():
    print BARL,"Testing Extremes",RBAR

    #////////////////////////////////////////
    #MAKE A MAP
    #////////////////////////////////////////
    ieffs=linspace(0*DEG,90*DEG,30)
    teffs=linspace(0.0*DEG,40*DEG,30)
    #ieffs=array([50.0*DEG])
    #teffs=array([50.0*DEG])
    IS,TS=meshgrid(ieffs,teffs)
    DIS=zeros_like(IS)

    a=SmallEllipse.a
    B=0.4
    i=0
    for ieff in ieffs:
        print "Testing ieff = ",ieff*RAD
        SmallEllipse.b=Re*cos(ieff)
        j=0
        for teff in teffs:
            SmallEllipse.cost=cos(teff)
            SmallEllipse.sint=sin(teff)
            b=SmallEllipse.b
            x=a*cos(teff)+sqrt(1-a**2*(sin(teff)-B/a)**2*\
                                   (1-b**2/a))
            SmallEllipse.C[0]=x
            SmallEllipse.C[1]=B
            dc,df=extremePoints(SmallEllipse)
            print dc,df
            DIS[i,j]=abs(1-dc)*100
            #print abs(1-dc)*100;exit(0)
            j+=1
        i+=1

    #////////////////////////////////////////
    #PLOT MAP
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,6))
    ax=fig.gca()
    dmin=DIS.min()
    dmax=DIS.max()
    levels=linspace(dmin,dmax,100)
    c=ax.contourf(IS*RAD,TS*RAD,transpose(DIS),
                  levels=levels)
    ax.set_xlabel(r"$i$",fontsize=20)
    ax.set_ylabel(r"$\theta$",fontsize=20)
    #ax.set_title("b = %.2f, min = %.4f%%, max = %.2f%%"%(B,dmin,dmax))
    cbar=fig.colorbar(c)
    cbar.ax.set_ylabel(r"Analytical Contact Formula Error (%%$R_{\star})$",fontsize=14)
    yts=cbar.ax.get_yticks()
    yl=[]
    for yt in yts:
        yl+=["%.1f%%"%((dmin+yt*(dmax-dmin))*1.00)]
    cbar.ax.set_yticklabels(yl)
    fig.savefig("plots/test_analytic_contact_contour.png")
    
    #////////////////////////////////////////
    #TEST DISTANCE
    #////////////////////////////////////////
    ieff=80*DEG
    teff=40.0*DEG

    SmallEllipse.b=Re*cos(ieff)
    SmallEllipse.cost=cos(teff)
    SmallEllipse.sint=sin(teff)
    a=SmallEllipse.a
    b=SmallEllipse.b
    B=0.2
    x=a*cos(teff)+sqrt(1-a**2*(sin(teff)-B/a)**2*\
                           (1-b**2/a))
    SmallEllipse.C[0]=x
    SmallEllipse.C[1]=B
    dc,df=extremePoints(SmallEllipse)
    #print "Extremes: ",dc,df,abs(dc-1)*100,"%"

    #////////////////////////////////////////
    #PLOT ELLIPSES
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()
    
    plotEllipse(ax,Figure(AR(0,0),dc,dc,1,0,''),
                color='r',linestyle='--')
    plotEllipse(ax,Figure(AR(0,0),df,df,1,0,''),
                color='b',linestyle='--')
    plotEllipse(ax,UnitaryCircle)
    plotEllipse(ax,SmallEllipse)
    
    rng=1.5
    xmin=SmallEllipse.C[0]-rng*Re;
    xmax=SmallEllipse.C[0]+rng*Re
    ymin=SmallEllipse.C[1]-rng*Re;
    ymax=SmallEllipse.C[1]+rng*Re
    #xmin=-1.5;xmax=1.5
    #ymin=-1.5;ymax=1.5
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("plots/test_analytic_contact.png")

def test_intersect():
    print BARL,"Testing Intersection",RBAR
    #////////////////////////////////////////
    #TEST INTERSECTION
    #////////////////////////////////////////
    Ps=eIcNumerical(SmallEllipse,UnitaryCircle)
    Ps=eIcNumerical(SmallEllipse,UnitaryCircle)

    #////////////////////////////////////////
    #PLOT CHARACTERIZING FUNCTION
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()
    Es=linspace(0,2*pi)
    ax.plot(Es*RAD,trigFunc(Es,SmallEllipse),
            label='Function',linewidth=2)
    ax.plot(Es*RAD,trigFuncD(Es,SmallEllipse),
            label='Derivarive')
    ax.plot(Es*RAD,trigFuncD2(Es,SmallEllipse),
            label='Curvature')
    ax.axhline(0,color='k',linewidth=2)
    for E in Ps:
        ax.axvline(E*RAD,linestyle='--')
    ax.legend(loc='best')
    ax.grid()
    fig.savefig("plots/test_intersect_charact.png")

    #////////////////////////////////////////
    #PLOT ELLIPSES
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    for E in Ps:
        P=toPoint(ellipsePointEcc(SmallEllipse,cos(E),sin(E)))
        plotPoint(ax,P)
        
    plotEllipse(ax,UnitaryCircle)
    plotEllipse(ax,SmallEllipse)
    
    rng=1.5
    ax.set_xlim((SmallEllipse.C[0]-rng*Re,SmallEllipse.C[0]+rng*Re))
    ax.set_ylim((SmallEllipse.C[1]-rng*Re,SmallEllipse.C[1]+rng*Re))
    ax.grid()
    fig.savefig("plots/test_intersect_conf.png")

###################################################
#ROUTINE TO TEST
###################################################
#test_intersect()
#test_inclusion()
#test_points()
#test_areas()
#test_area_oblate()
#test_area_strip()
#it,itf=test_area_striping(Planet,Ringext,Ringint)
"""
its=[]
itfs=[]
for x in linspace(0.0,1.0,10):
    Planet.C=Ringext.C=Ringint.C=AR(x,0)
    it,itf=test_area_striping(Planet,Ringext,Ringint)
    its+=[it]
    itfs+=[itf]
#"""
test_analytic_contact()
#test_analytic_area()
