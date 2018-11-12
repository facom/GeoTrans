from geotrans import *
from system import *

###################################################
#TEST ROUTINES
###################################################
def test_transit_curve_limb():
    print BARL,"Testing Transit Limb Darkening Curve",RBAR
    
    #////////////////////////////////////////
    #TRANSIT CURVE
    #////////////////////////////////////////
    #SPHERICAL PLANET
    SphericalPlanet=copyObject(System)
    SphericalPlanet.fp=0.0
    derivedSystemProperties(SphericalPlanet)
    updatePlanetRings(SphericalPlanet)
    SphericalPlanet.Ringext.b=SphericalPlanet.Ringint.b=0.0

    #OBLATE PLANET
    OblatePlanet=copyObject(System)
    OblatePlanet.fp=0.5
    derivedSystemProperties(OblatePlanet)
    updatePlanetRings(OblatePlanet)
    OblatePlanet.Ringext.b=OblatePlanet.Ringint.b=0.0

    #CONTACT TIMES
    tcsp=contactTimes(SphericalPlanet)
    tcso=contactTimes(OblatePlanet)
    tcsr=contactTimes(System)
    
    #EXTREMA
    tmin=tcsr[1]-System.dtplanet
    tmax=tcsr[4]+System.dtplanet

    #LIGHT CURVE RANGE
    maxd=(System.Rp**2+System.Re**2*cos(System.ieff)-System.Ri**2*cos(System.ieff))
    targs=dict(verticalalignment='top',
               horizontalalignment='center',
               bbox=dict(fc='w',ec='none'))

    #NUMBER OF POINTS
    N=100
    n=N/5
    
    #RINGED PLANET AREA
    Art=ringedPlanetArea(System)
    Apt=FIGUREAREA(SphericalPlanet.Planet)
    Apo=FIGUREAREA(OblatePlanet.Planet)
    
    print "Area of the planet spherical (star area units) = ",Apt
    print "Area of the planet oblate (star area units) = ",Apo
    print "Area of the ringed planet (star area units) = ",Art

    #CURVE
    ts=linspace(tmin,tmax,N)
    Ats=[]
    Aps=[]
    Aos=[]
    Flts=[]
    Flps=[]
    Flos=[]
    i=1
    for t in ts:
        At=transitAreaTime(t,System)
        Ap=transitAreaTime(t,SphericalPlanet)
        Ao=transitAreaOblateTime(t,OblatePlanet)
        Ats+=[1-At/pi]
        Aps+=[1-Ap/pi]
        Aos+=[1-Ao/pi]

        Flt=fluxLimbTime(t,Art,System)
        Flts+=[Flt]
        Flp=fluxLimbTime(t,Apt,SphericalPlanet)
        Flps+=[Flp]
        Flo=fluxLimbTime(t,Apt,OblatePlanet,
                         areas=areaOblateStriping)
        Flos+=[Flo]

        if (i%n)==0:
            print "i = %d, t = %e, At = %e, Ap = %e, Ao = %e"%\
                (i,(t-System.tcen)/HOUR,At,Ap,Ao)
        i+=1

    #////////////////////////////////////////
    #PLOT
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()
    ax2=fig.add_axes([0.4,0.6,0.2,0.2])

    ax.plot((ts-System.tcen)/HOUR,Ats,'b-')
    ax.plot((ts-System.tcen)/HOUR,Flts,'ko-',
            linewidth=2,markersize=3)
    ax.plot((ts-System.tcen)/HOUR,Aps,'g-')
    ax.plot((ts-System.tcen)/HOUR,Flps,'ko-',
            linewidth=2,markersize=3)
    ax.plot((ts-System.tcen)/HOUR,Aos,'g-')
    ax.plot((ts-System.tcen)/HOUR,Flos,'ko-',
            linewidth=2,markersize=3)

    ax.axvline((tcsr[1]-System.tcen)/HOUR,
               color='b',linestyle='-')
    ax.axvline((tcsr[2]-System.tcen)/HOUR,
               color='b',linestyle='--')
    ax.axvline((tcsr[3]-System.tcen)/HOUR,
               color='b',linestyle='--')
    ax.axvline((tcsr[4]-System.tcen)/HOUR,
               color='b',linestyle='-')

    ax.axvline((tcsp[1]-SphericalPlanet.tcen)/HOUR,
               color='g',linestyle='-')
    ax.axvline((tcsp[2]-SphericalPlanet.tcen)/HOUR,
               color='g',linestyle='--')
    ax.axvline((tcsp[3]-SphericalPlanet.tcen)/HOUR,
               color='g',linestyle='--')
    ax.axvline((tcsp[4]-SphericalPlanet.tcen)/HOUR,
               color='g',linestyle='-')

    ax.axvline(-4.370816e+00,
               color='k',linestyle=':')

    """
    ax.axvline((tcsp[1]-OblatePlanet.tcen)/HOUR,
               color='r',linestyle='-')
    ax.axvline((tcsp[2]-OblatePlanet.tcen)/HOUR,
               color='r',linestyle='--')
    ax.axvline((tcsp[3]-OblatePlanet.tcen)/HOUR,
               color='r',linestyle='--')
    ax.axvline((tcsp[4]-OblatePlanet.tcen)/HOUR,
               color='r',linestyle='-')
               """
    OblatePlanet.Planet.C=System.Planet.C=System.Ringext.C=System.Ringint.C=AR(0,0)
    plotEllipse(ax2,OblatePlanet.Planet,color='b')
    plotEllipse(ax2,System.Ringext,color='k')
    plotEllipse(ax2,System.Ringint,color='k')
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_xlim((-System.Re,System.Re))
    ax2.set_ylim((-System.Re,System.Re))
    
    ymin=1.0-10*maxd/9
    ymax=1.0+maxd/10
    tmin=(tmin-System.tcen)/HOUR
    tmax=(tmax-System.tcen)/HOUR
    ax.set_xlim((tmin,tmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("plots/test_transit_curve_limb.png")

def test_limb_darkening():
    print BARL,"Testing Limb Darkening",RBAR
    
    #////////////////////////////////////////
    #PLOT
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()
    rhos=linspace(0.0,1.0,100)
    ls=limbDarkeningPhysical(rhos,0.7,-0.23)
    ax.plot(rhos,ls,label='Sun')
    ls=limbDarkeningPhysical(rhos,0.64,-0.055)
    ax.plot(rhos,ls,label='Brown et al. (2001)')

    ax.set_xlim((0,1))
    ax.set_title("Limb darkening",fontsize=16)
    ax.set_xlabel(r"$\rho/R_\star$",fontsize=14)
    ax.set_ylabel(r"$I(\rho)/I(0)$",fontsize=14)
    ax.legend(loc='best')
    fig.savefig("plots/test_limb_darkening.png")

    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()
    rhos=linspace(0.0,1.0,100)
    ls=limbDarkeningNormalized(rhos,0.7,-0.23)
    ax.plot(rhos,ls,label='Sun')
    ls=limbDarkeningNormalized(rhos,0.64,-0.055)
    ax.plot(rhos,ls,label='Brown et al. (2001)')

    ax.set_xlim((0,1))
    ax.set_title("Limb darkening",fontsize=16)
    ax.set_xlabel(r"$\rho/R_\star$",fontsize=14)
    ax.set_ylabel(r"$i(\rho)=I(\rho)/\int_0^1 2\pi\rho I(\rho) d\rho$",fontsize=10)
    ax.legend(loc='best')
    fig.savefig("plots/test_limb_darkening_normalized.png")

def test_transit_curve():
    print BARL,"Testing Transit Curve",RBAR
    
    #////////////////////////////////////////
    #TRANSIT CURVE
    #////////////////////////////////////////
    #COPY
    NotRinged=copyObject(System)
    NotRinged.Ringext.b=NotRinged.Ringint.b=0.0
    
    #CONTACT TIMES
    tcsp=contactTimes(NotRinged)
    tcsr=contactTimes(System)

    #EXTREMA
    tmin=tcsr[1]-System.dtplanet
    tmax=tcsr[4]+System.dtplanet
    
    #LIGHT CURVE RANGE
    maxd=(System.Rp**2+System.Re**2*cos(System.ieff)-System.Ri**2*cos(System.ieff))
    targs=dict(verticalalignment='top',
               horizontalalignment='center',
               bbox=dict(fc='w',ec='none'))

    #NUMBER OF POINTS
    N=100
    n=N/5
    
    #RINGED PLANET AREA
    Art=ringedPlanetArea(System)
    Apt=pi*System.Rp**2
    
    print "Area of the planet alone (star area units) = ",Apt
    print "Area of the ringed planet (star area units) = ",Art

    #CURVE
    ts=linspace(tmin,tmax,N)
    Ats=[]
    Aps=[]
    i=1
    for t in ts:
        At=transitAreaTime(t,System)
        Ap=transitAreaTime(t,NotRinged)
        Ats+=[1-At/pi]
        Aps+=[1-Ap/pi]
        if (i%n)==0:
            print "i = %d, t = %e, At = %e"%\
                (i,(t-System.tcen)/HOUR,At)
        i+=1

    #////////////////////////////////////////
    #PLOT
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()
    ax2=fig.add_axes([0.4,0.6,0.2,0.2])

    ax.plot((ts-System.tcen)/HOUR,Ats,'b-')
    ax.plot((ts-System.tcen)/HOUR,Aps,'g-')

    ax.axvline((tcsr[1]-System.tcen)/HOUR,
               color='b',linestyle='-')
    ax.axvline((tcsr[2]-System.tcen)/HOUR,
               color='b',linestyle='--')
    ax.axvline((tcsr[3]-System.tcen)/HOUR,
               color='b',linestyle='--')
    ax.axvline((tcsr[4]-System.tcen)/HOUR,
               color='b',linestyle='-')

    ax.axvline((tcsp[1]-NotRinged.tcen)/HOUR,
               color='g',linestyle='-')
    ax.axvline((tcsp[2]-NotRinged.tcen)/HOUR,
               color='g',linestyle='--')
    ax.axvline((tcsp[3]-NotRinged.tcen)/HOUR,
               color='g',linestyle='--')
    ax.axvline((tcsp[4]-NotRinged.tcen)/HOUR,
               color='g',linestyle='-')

    System.Planet.C=System.Ringext.C=System.Ringint.C=AR(0,0)
    plotEllipse(ax2,System.Planet,color='b')
    plotEllipse(ax2,System.Ringext,color='k')
    plotEllipse(ax2,System.Ringint,color='k')
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_xlim((-System.Re,System.Re))
    ax2.set_ylim((-System.Re,System.Re))
    
    ymin=1.0-10*maxd/9
    ymax=1.0+maxd/10
    tmin=(tmin-System.tcen)/HOUR
    tmax=(tmax-System.tcen)/HOUR
    ax.set_xlim((tmin,tmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("plots/test_transit_curve.png")

def test_contacts():
    print BARL,"Testing Contact Times",RBAR
    #////////////////////////////////////////
    #CONTACT TIMES
    #////////////////////////////////////////
    tcs=contactTimes(System)
    print "Transit duration = %.3f h"%((tcs[-1]-tcs[1])/HOUR)

    #////////////////////////////////////////
    #PLOT
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()
    
    #PLOT FIGURES
    cl=['b','r','r','b']
    ls=['-','--','--','-']
    i=1
    for t in tcs[1:]:
        updatePosition(System,t)
        args=dict(color=cl[i-1],linestyle=ls[i-1])
        plotEllipse(ax,System.Star,patch=True,color='y')
        plotEllipse(ax,System.Planet,**args)
        plotEllipse(ax,System.Ringext,**args)
        plotEllipse(ax,System.Ringint,**args)
        ax.text(System.Planet.C[0],
                System.Planet.C[1]+System.Re,
                r"$t_%d=%.2f$"%(i,(t-System.tcen)/HOUR),
                horizontalalignment='center')
        i+=1
    
    rng=1.5
    xmin=-rng;xmax=rng
    ymin=-rng;ymax=rng
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("plots/test_contacts.png")
    
    
def test_system():
    print BARL,"Testing System Properties",RBAR
    #systemShow(System)

    #////////////////////////////////////////
    #GETTING POSITION AT A GIVEN TIME
    #////////////////////////////////////////
    t=System.tcen-System.dtrans/2+System.dtplanet
    updatePosition(System,t)

    #////////////////////////////////////////
    #PLOT
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    plotEllipse(ax,System.Star,patch=True,color='y')
    plotEllipse(ax,System.Planet,color='b')
    plotEllipse(ax,System.Ringext,color='b')
    plotEllipse(ax,System.Ringint,color='r')

    rng=1.5
    xmin=System.Planet.C[0]-rng*System.Re;
    xmax=System.Planet.C[0]+rng*System.Re
    ymin=System.Planet.C[1]-rng*System.Re;
    ymax=System.Planet.C[1]+rng*System.Re
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("plots/test_system.png")

    #////////////////////////////////////////
    #ORBIT
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    plotEllipse(ax,System.Star,patch=True,color='y')
    plotEllipse(ax,System.Planet,color='b')

    #PLOT ORBIT
    fs=linspace(System.fcen-1*DEG,System.fcen+1*DEG,100)
    rs=ellipseRadiusE(System.ap/System.Rstar,System.ep,fs)
    xs=rs*cos(fs)
    ys=rs*sin(fs)
    zs=zeros_like(xs)
    rs=array([dot(System.Mos,AR3(x,y,z)) for x,y,z in zip(xs,ys,zs)])
    ax.plot(rs[:,0],rs[:,1],'b-')

    fs=linspace(pi+System.fcen-1*DEG,pi+System.fcen+1*DEG,100)
    rs=ellipseRadiusE(System.ap/System.Rstar,System.ep,fs)
    xs=rs*cos(fs)
    ys=rs*sin(fs)
    zs=zeros_like(xs)
    rs=array([dot(System.Mos,AR3(x,y,z)) for x,y,z in zip(xs,ys,zs)])
    ax.plot(rs[:,0],rs[:,1],'b--',zorder=-10)

    #PLOT RINGS
    fs=linspace(0,2*pi,100)
    xs=System.Re*cos(fs)
    ys=System.Re*sin(fs)
    zs=zeros_like(xs)
    rs=array([dot(System.Mrs,AR3(x,y,z))+System.rs for x,y,z in zip(xs,ys,zs)])
    for i in xrange(len(rs)):
        c='k'
        if (rs[i,2]-System.rs[2])<0:c='r'
        ax.plot([rs[i,0]],[rs[i,1]],'o',markeredgecolor='none',color=c,markersize=2)

    #PLOT AXIS
    rx=dot(System.Mrs,[1.0,0.0,0.0])
    ry=dot(System.Mrs,[0.0,1.0,0.0])
    rz=dot(System.Mrs,[0.0,0.0,1.0])

    rp=3*System.Rp*rz+System.rs
    ax.plot([System.rs[0],rp[0]],
            [System.rs[1],rp[1]],'k-',linewidth=3)
    ax.text(rp[0],rp[1],'z')
    
    rp=3*System.Rp*ry+System.rs
    ax.plot([System.rs[0],rp[0]],
            [System.rs[1],rp[1]],'k-',linewidth=3)
    ax.text(rp[0],rp[1],'y')
    
    rp=3*System.Rp*rx+System.rs
    ax.plot([System.rs[0],rp[0]],
            [System.rs[1],rp[1]],'k-',linewidth=3)
    ax.text(rp[0],rp[1],'x')

    xmin=-1.5;xmax=1.5
    ymin=-1.5;ymax=1.5
    rng=1.5
    xmin=System.Planet.C[0]-rng*System.Re;
    xmax=System.Planet.C[0]+rng*System.Re
    ymin=System.Planet.C[1]-rng*System.Re;
    ymax=System.Planet.C[1]+rng*System.Re
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("plots/test_system_orbit.png")

###################################################
#ROUTINE TO TEST
###################################################
#test_system()
#test_contacts()
#test_transit_curve()
#test_limb_darkening()
test_transit_curve_limb()
