from geotrans import *
from system import *

#########################################
#SYSTEM
#########################################
Ringed=System
NotRinged=copyObject(Ringed)
NotRinged.Ringext.b=NotRinged.Ringint.b=0.0

VAR=1
FIX=0

SHOW=1
HIDE=0

DEF=0
MIN=1
MAX=2
FUNC=3
SCAL=4
STAT=5

IDENT=lambda x:x
POW10=lambda x:10**x

PARAMETERS=dict(
    
    #PLANETARY RADIUS: SATURN/FSTAR, JUPITER/MSTAR
    Rplanet=[RSAT,
             RSAT,
             RJUP,
             IDENT,
             "S.Rstar",
             FIX],
    
    #STELLAR MASS
    Mstar=[1.0*MSUN,
           0.6*MSUN,
           1.2*MSUN,
           IDENT,
           "MSUN",
           FIX],
    
    #ORBITAL SEMI MAJOR AXIS
    ap=[1.0*AU,
        1.0*AU,
        3.0*AU,
        IDENT,
        "AU",
        FIX],

    #ECCENTRICITY
    ep=[0.0,
        0.0,
        0.5,
        IDENT,
        "1",
        FIX],

    #ORBITAL INCLINATION
    iorb=[0.0,
          0.0,
          1.0,
          arccos,
          "DEG",
          FIX],

    #ECCENTRICITY
    wp=[0.0*DEG,
        0.0*DEG,
        360.0*DEG,
        IDENT,
        "DEG",
        FIX],
    
    #RING RADIUS 
    fe=[2.35,
        2.00,
        4.00,
        IDENT,
        "1",
        FIX],
    
    fi=[1.58,
        1.50,
        2.00,
        IDENT,
        "1",
        FIX],

    #ORBITAL INCLINATION
    ir=[1.0,
        0.0,
        1.0,
        arccos,
        "DEG",
        FIX],
    
    #ROLL
    phir=[0.0*DEG,
          #-90.0*DEG,
          0.0*DEG,
          #+90.0*DEG,
          360.0*DEG,
          IDENT,
          "DEG",
          FIX],

    #OPACITY
    tau=[log10(4.0),
         log10(1.0),
         log10(4.0),
         POW10,
         "1",
         FIX],

    )
PARKEYS=sorted(PARAMETERS.keys())

PROPERTIES=dict(
    ieff=[0,
          0,
          0,
          IDENT,
          "DEG",
          SHOW],
    teff=[0,
          0,
          0,
          IDENT,
          "DEG",
          SHOW],
    r=[0,
       0,
       0,
       IDENT,
       "1",
       SHOW],

    p=[0,
       0,
       0,
       IDENT,
       "1",
       SHOW],

    PR=[0,
        0,
        0,
        IDENT,
        "1",
        SHOW],

    logPR=[0,
           0,
           0,
           IDENT,
           "1",
           SHOW],
    )
PROPKEYS=sorted(PROPERTIES.keys())

def transitPosition(Rp,fe,i,t,B,direction=+1,sign=+1,qcorrected=False):
    """
    direction = +1 (out of disk), -1 (toward disk)
    sign: -1 (before contact), +1 (after contact)
    
    Example:
      Contact 1: direction=-1, sign=-1
      Contact 2: direction=-1, sign=+1
      Contact 3: direction=+1, sign=-1
      Contact 4: direction=+1, sign=+1
    """
    a=fe*Rp
    b=fe*Rp*cos(i)

    if qcorrected:
        if cos(i)>0.6:
            xp=direction*sqrt((1+direction*sign*a)**2-B**2)
            return xp

    a=fe*Rp
    b=fe*Rp*cos(i)

    xp=direction*sqrt(1-a**2*(sin(t)-sign*B/a)**2*\
                          (1-b**2/a))+\
                          sign*a*cos(t)

    #COMPARE WITH THE NOT-RINGED CASE
    xpP=direction*sqrt((1+direction*sign*Rp)**2-B**2)
    if sign<0:
        if xpP<xp:xp=xpP
    else:
        if xpP>xp:xp=xpP
    return xp

def testTransitDepth():
    print BARL,"Test Transit Depth",RBAR

    #========================================
    #FIX RINGED PROPERTIES BY HAND
    #========================================
    #MANUAL i,t
    i=60.0*DEG;t=30*DEG
    
    Ringed.ieff=i;Ringed.teff=t

    Ringed.Ringext.b=Ringed.Ringext.a*cos(i)
    Ringed.Ringext.cost=cos(t);Ringed.Ringext.sint=sin(t)

    Ringed.Ringint.b=Ringed.Ringint.a*cos(i)
    Ringed.Ringint.cost=cos(t);Ringed.Ringint.sint=sin(t)

    Ringed.block=blockFactor(Ringed.tau,i)

    #========================================
    #NOT RINGED TRANSIT DEPTH
    #========================================
    Anr=pi*NotRinged.Rp**2
    print TAB,"Transit area (not ringed): %.17e"%Anr

    #========================================
    #ANALYTICAL TRANSIT DEPTH
    #========================================
    Aarg=analyticalTransitArea(Ringed.Rp,Ringed.block,Ringed.fi,Ringed.fe,Ringed.ieff)
    print TAB,"Analytical Transit area (ringed): %.17e"%Aarg
    
    #========================================
    #RINGED TRANSIT DEPTH
    #========================================
    Arg=ringedPlanetArea(Ringed)
    print TAB,"Transit area (ringed): %.17e"%Arg
    r=sqrt(Arg/Anr)
    print TAB,"Ratio of depths: %.17e"%(Arg/Anr)
    print TAB,"Ratio of radii: %.17e"%(r)

    #========================================
    #MONTECARLO AREA
    #========================================
    NP=10000
    #"""
    Am,dA,xs,ys=transitAreaMontecarlo(Ringed.Planet,
                                      Ringed.Ringext,
                                      Ringed.Ringint,
                                      NP=NP)
    print TAB,"Montecarlo area: %.17e +/- %.1e"%(Am,dA)
    #"""

    #========================================
    #PLOT
    #========================================
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    ax.plot(xs,ys,'ro',markersize=1)
    plotEllipse(ax,Ringed.Star,color='y')
    plotEllipse(ax,Ringed.Planet,color='b')
    plotEllipse(ax,Ringed.Ringext,color='k')
    plotEllipse(ax,Ringed.Ringint,color='r')
    
    rng=1.5
    Re=Ringed.Ringext.a
    #Re=1.0
    xmin=Ringed.Planet.C[0]-rng*Re;
    xmax=Ringed.Planet.C[0]+rng*Re
    ymin=Ringed.Planet.C[1]-rng*Re;
    ymax=Ringed.Planet.C[1]+rng*Re
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("figures/TestAreas.png")

def testTransitDuration():

    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    #========================================
    #FIX RINGED PROPERTIES BY HAND
    #========================================
    #MANUAL i,t
    i=80.0*DEG;t=60*DEG
    
    Ringed.ieff=i;Ringed.teff=t

    Ringed.Ringext.b=Ringed.Ringext.a*cos(i)
    Ringed.Ringext.cost=cos(t);Ringed.Ringext.sint=sin(t)

    Ringed.Ringint.b=Ringed.Ringint.a*cos(i)
    Ringed.Ringint.cost=cos(t);Ringed.Ringint.sint=sin(t)

    Ringed.block=blockFactor(Ringed.tau,i)

    print BARL,"Test Transit Duration",RBAR

    print "Orientation parameters:"
    print TAB,"i = %.2f deg"%(i*RAD)
    print TAB,"t = %.2f deg"%(t*RAD)

    #========================================
    #PROPERTIES
    #========================================
    RingedC=copyObject(Ringed)
    ap=Ringed.ap/Ringed.Rstar
    P=Ringed.Porb/HOUR
    ip=Ringed.iorb
    
    #========================================
    #NOT RINGED TRANSIT DURATION (NUMERICAL)
    #========================================
    tcsp=contactTimes(NotRinged)
    tT=(tcsp[-1]-tcsp[1])/HOUR
    tF=(tcsp[-2]-tcsp[2])/HOUR

    """
    updatePosition(NotRinged,tcsp[1])
    plotEllipse(ax,NotRinged.Planet,color='b',
                linestyle='-',linewidth=1)
    updatePosition(NotRinged,tcsp[-1])
    plotEllipse(ax,NotRinged.Planet,color='b',
                linestyle='-',linewidth=1)
    #"""

    print TAB,"Transit duration numerical (not ringed):"
    print 2*TAB,"Full: %.17e"%tT
    print 2*TAB,"Total: %.17e"%tF

    #========================================
    #NOT RINGED TRANSIT DURATION (ANALYTICAL)
    #========================================
    print TAB,"Transit duration analytical (not ringed):"

    xp=sqrt((1+NotRinged.Rp)**2-NotRinged.Borb**2)
    xm=sqrt((1-NotRinged.Rp)**2-NotRinged.Borb**2)

    """
    CP1=Figure(AR(-xp,Ringed.Borb),NotRinged.Rp,NotRinged.Rp,
              1.0,0.0,'Contact 1')
    CP4=Figure(AR(+xp,Ringed.Borb),NotRinged.Rp,NotRinged.Rp,
              1.0,0.0,'Contact 4')
    plotEllipse(ax,CP1,color='k',linestyle=':',linewidth=2)
    plotEllipse(ax,CP4,color='k',linestyle=':',linewidth=2)
    #"""

    tT=P*arcsin(2*xp/(ap*sin(ip)))/(2*pi)
    tF=P*arcsin(2*xm/(ap*sin(ip)))/(2*pi)

    print 2*TAB,"Full: %.17e"%tT
    print 2*TAB,"Total: %.17e"%tF

    xp1=-xp
    xp2=-xm
    xp3=xm
    xp4=xp
    print 2*TAB,"Contact point:"
    print 3*TAB,"xp1 = %.17e"%xp1
    print 3*TAB,"xp2 = %.17e"%xp2
    print 3*TAB,"xp3 = %.17e"%xp3
    print 3*TAB,"xp4 = %.17e"%xp4
    
    #========================================
    #RINGED TRANSIT DURATION (NUMERICAL)
    #========================================
    lw=1
    print TAB,"Transit duration numerical (ringed):"

    tcsp=contactTimes(RingedC)
    tT=(tcsp[-1]-tcsp[1])/HOUR
    tF=(tcsp[-2]-tcsp[2])/HOUR

    updatePosition(RingedC,tcsp[1])
    xp1=RingedC.C[0]
    plotEllipse(ax,RingedC.Ringext,color='r',
                linestyle='-',linewidth=lw)
    plotEllipse(ax,RingedC.Planet,color='r',
                linestyle='-',linewidth=lw)

    updatePosition(RingedC,tcsp[2])
    xp2=RingedC.C[0]
    plotEllipse(ax,RingedC.Ringext,color='r',
                linestyle='-',linewidth=lw)
    plotEllipse(ax,RingedC.Planet,color='r',
                linestyle='-',linewidth=lw)

    updatePosition(RingedC,tcsp[3])
    xp3=RingedC.C[0]
    plotEllipse(ax,RingedC.Ringext,color='r',
                linestyle='-',linewidth=lw)
    plotEllipse(ax,RingedC.Planet,color='r',
                linestyle='-',linewidth=lw)

    updatePosition(RingedC,tcsp[4])
    xp4=RingedC.C[0]
    plotEllipse(ax,RingedC.Ringext,color='r',
                linestyle='-',linewidth=lw)
    plotEllipse(ax,RingedC.Planet,color='r',
                linestyle='-',linewidth=lw)

    print 2*TAB,"Full: %.17e"%tT
    print 2*TAB,"Total: %.17e"%tF

    print 2*TAB,"Contact point:"
    print 3*TAB,"xp1 = %.17e"%xp1
    print 3*TAB,"xp2 = %.17e"%xp2
    print 3*TAB,"xp3 = %.17e"%xp3
    print 3*TAB,"xp4 = %.17e"%xp4

    #========================================
    #RINGED TRANSIT ANALYTICAL 
    #========================================
    lw=2
    print TAB,"Transit duration analytical (ringed):"

    a=RingedC.Ringext.a
    b=RingedC.Ringext.b
    B=RingedC.Borb

    #C1
    xpa1=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                         direction=-1,sign=-1)
    #C2
    xpa2=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                         direction=-1,sign=+1)
    #C3
    xpa3=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                         direction=+1,sign=-1)
    #C4
    xpa4=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                         direction=+1,sign=+1)

    taT=P*arcsin((xpa4-xpa1)/(ap*sin(ip)))/(2*pi)
    taF=P*arcsin((xpa3-xpa2)/(ap*sin(ip)))/(2*pi)

    print 2*TAB,"Full: %.17e"%taT
    print 2*TAB,"Total: %.17e"%taF

    print 2*TAB,"Contact points:"
    print 3*TAB,"xpa1 = %.17e"%xpa1
    print 3*TAB,"xpa2 = %.17e"%xpa2
    print 3*TAB,"xpa3 = %.17e"%xpa3
    print 3*TAB,"xpa4 = %.17e"%xpa4

    RingedC.Ringext.C[0]=xpa1
    RingedC.Planet.C[0]=xpa1
    plotEllipse(ax,RingedC.Planet,color='g',
                linestyle='-',linewidth=lw)
    plotEllipse(ax,RingedC.Ringext,color='g',
                linestyle='-',linewidth=lw)

    RingedC.Ringext.C[0]=xpa2
    RingedC.Planet.C[0]=xpa2
    plotEllipse(ax,RingedC.Planet,color='g',
                linestyle='-',linewidth=lw)
    plotEllipse(ax,RingedC.Ringext,color='g',
                linestyle='-',linewidth=lw)

    RingedC.Ringext.C[0]=xpa3
    RingedC.Planet.C[0]=xpa3
    plotEllipse(ax,RingedC.Planet,color='g',
                linestyle='-',linewidth=lw)
    plotEllipse(ax,RingedC.Ringext,color='g',
                linestyle='-',linewidth=lw)

    RingedC.Ringext.C[0]=xpa4
    RingedC.Planet.C[0]=xpa4
    plotEllipse(ax,RingedC.Planet,color='g',
                linestyle='-',linewidth=lw)
    plotEllipse(ax,RingedC.Ringext,color='g',
                linestyle='-',linewidth=lw)

    #========================================
    #ERROR IN CONTACT POSITIONS
    #========================================
    dxp1=abs(xp1-xpa1)/abs(xpa1)*100
    dxp2=abs(xp2-xpa2)/abs(xpa2)*100
    dxp3=abs(xp3-xpa3)/abs(xpa3)*100
    dxp4=abs(xp4-xpa4)/abs(xpa4)*100
    print 1*TAB,"Error in Contact points:"
    print 2*TAB,"dxp1 = %.4f %%"%dxp1
    print 2*TAB,"dxp2 = %.4f %%"%dxp2
    print 2*TAB,"dxp3 = %.4f %%"%dxp3
    print 2*TAB,"dxp4 = %.4f %%"%dxp4

    dtT=abs(tT-taT)*60
    dtF=abs(tF-taF)*60
    print 1*TAB,"Error in Contact points:"
    print 2*TAB,"dtT = %.2f min"%dtT
    print 2*TAB,"dtF = %.2f min"%dtF
    
    #========================================
    #PLOT
    #========================================
    plotEllipse(ax,Ringed.Star,color='y')
    plotEllipse(ax,Ringed.Planet,color='b')
    plotEllipse(ax,Ringed.Ringext,color='k')
    plotEllipse(ax,Ringed.Ringint,color='r')
    
    rng=1.5
    Re=Ringed.Ringext.a
    Re=1.0
    xmin=Ringed.Planet.C[0]-rng*Re;
    xmax=Ringed.Planet.C[0]+rng*Re
    ymin=Ringed.Planet.C[1]-rng*Re;
    ymax=Ringed.Planet.C[1]+rng*Re
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("figures/TestDuration.png")

def errorTransitPositions():

    print BARL,"Testing Extremes",RBAR

    #////////////////////////////////////////
    #DATA
    #////////////////////////////////////////
    Rp=0.0836255747894
    fi=2.0
    fe=2.5
    Ri=fi*Rp
    Re=fe*Rp
    ieff=0.0*DEG
    teff=0.0*DEG

    C=AR(0.9,0.6)
    SmallEllipse=Figure(C,
                        Re,Re*cos(ieff),
                        cos(teff),sin(teff),
                        'SmallEllipse')

    UnitaryCircle=Figure(AR(0,0),
                         1.0,1.0,
                         1.0,0.0,
                         'BigCircle')

    #////////////////////////////////////////
    #MAKE A MAP
    #////////////////////////////////////////
    ieffs=linspace(0*DEG,90*DEG,30)
    teffs=linspace(0.0*DEG,40*DEG,30)
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
            #print dc,df
            DIS[i,j]=abs(1-dc)*cos(teff)
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
    ax.set_title("$b = %.2f$"%(B),position=(0.5,1.02),
                 fontsize=16)
    cbar=fig.colorbar(c)
    cbar.ax.set_ylabel(r"$\Delta x_{\pm}$(%)",fontsize=14)
    yts=cbar.ax.get_yticks()
    yl=[]
    for yt in yts:
        yl+=["%.1f%%"%((dmin+yt*(dmax-dmin))*100.0)]
    cbar.ax.set_yticklabels(yl)
    fig.savefig("figures/ErrorTransitPositionsAnalytic.png")

def errorTransitTimes():

    qload=False
    #qload=True
    #qcorrected=True
    qcorrected=False
    print BARL,"Testing Analytical Times",RBAR

    #////////////////////////////////////////
    #SYSTEM
    #////////////////////////////////////////
    RingedC=copyObject(Ringed)
    ap=Ringed.ap/Ringed.Rstar
    P=Ringed.Porb/HOUR
    ip=Ringed.iorb
    B=RingedC.Borb

    #////////////////////////////////////////
    #MAKE A MAP
    #////////////////////////////////////////
    cieffmin=0.01
    cieffmax=1.0
    Ncieffs=30
    teffmin=0.0*DEG
    teffmax=90.0*DEG
    Nteffs=30

    if not qload:
        #NOT-RINGED
        tcsp=contactTimes(NotRinged)
        tPT=(tcsp[-1]-tcsp[1])/HOUR
        tPF=(tcsp[-2]-tcsp[2])/HOUR

        cieffs=linspace(cieffmin,cieffmax,Ncieffs)
        teffs=linspace(teffmin,teffmax,Nteffs)

        IS,TS=meshgrid(cieffs,teffs)
        DIST=zeros_like(IS)   
        DISF=zeros_like(IS)
        DISTP=zeros_like(IS)   
        DISFP=zeros_like(IS)
        DISR=zeros_like(IS)
 
        ii=0
        for ci in cieffs:
            i=arccos(ci)
            print "Testing ieff = ",i*RAD
            jj=0
            for t in teffs:
                RingedC.Ringext.b=RingedC.Ringext.a*cos(i)
                RingedC.Ringext.cost=cos(t)
                RingedC.Ringext.sint=sin(t)

                #NUMERICAL
                tcsp=contactTimes(RingedC)
                tT=(tcsp[-1]-tcsp[1])/HOUR
                tF=(tcsp[-2]-tcsp[2])/HOUR
                rt=(tT**2-tF**2)**1.5
                
                #ANALYTICAL
                xpa1=transitPosition(RingedC.Rp,RingedC.fe,
                                     i,t,B,
                                     direction=-1,sign=-1,
                                     qcorrected=qcorrected)
                xpa2=transitPosition(RingedC.Rp,RingedC.fe,
                                     i,t,B,
                                     direction=-1,sign=+1,
                                     qcorrected=qcorrected)
                xpa3=transitPosition(RingedC.Rp,RingedC.fe,
                                     i,t,B,
                                     direction=+1,sign=-1,
                                     qcorrected=qcorrected)
                xpa4=transitPosition(RingedC.Rp,RingedC.fe,
                                     i,t,B,
                                     direction=+1,sign=+1,
                                     qcorrected=qcorrected)

                taT=P*arcsin((xpa4-xpa1)/(ap*sin(ip)))/(2*pi)
                taF=P*arcsin((xpa3-xpa2)/(ap*sin(ip)))/(2*pi)
                rat=(taT**2-taF**2)**1.5

                dT=(taT-tT)
                dF=(taF-tF)
                dr=(rt-rat)
                dTp=(tT-tPT)

                """
                print TAB,"Testing teff = ",t*RAD
                print 3*TAB,"Numerical times: ",tT,tF
                print 3*TAB,"Analytical times: ",taT,taF
                print 3*TAB,"Errors: ",dT,dF
                #"""

                DIST[ii,jj]=abs(dT)/tT
                DISF[ii,jj]=abs(dF)/tF
                DISR[ii,jj]=dT/dTp
                if abs(dTp*60)<10:DISR[ii,jj]=0.0
                DISTP[ii,jj]=abs(dTp)/tT
                DISFP[ii,jj]=dTp*60
                
                jj+=1
            ii+=1

        savetxt("IS.dat",IS)
        savetxt("TS.dat",TS)
        savetxt("DIST.dat",DIST)
        savetxt("DISF.dat",DISF)
        savetxt("DISR.dat",DISR)
        savetxt("DISTP.dat",DISTP)
        savetxt("DISFP.dat",DISFP)

    IS=loadtxt("IS.dat")
    TS=loadtxt("TS.dat")
    DIST=loadtxt("DIST.dat")
    DISF=loadtxt("DISF.dat")
    DISR=loadtxt("DISR.dat")
    DISTP=loadtxt("DISTP.dat")
    DISFP=loadtxt("DISFP.dat")
    
    cmap=plt.get_cmap("bwr")
    cmap=plt.get_cmap("spectral")
    cmap=plt.get_cmap("PiYG_r")
    cmap=plt.get_cmap("rainbow")
    #////////////////////////////////////////
    #PLOT DIST
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,6))
    ax=fig.add_axes([0.1,0.1,0.8,0.8])

    dmin=DIST.min()
    dmax=DIST.max()

    levels=linspace(dmin,dmax,100)
    c=ax.contourf(IS,TS*RAD,transpose(DIST),
                  levels=levels,cmap=cmap)

    ax.set_xlabel(r"$\cos\,i$",fontsize=14)
    ax.set_ylabel(r"$\theta$",fontsize=14)
    #ax.set_title("Transit Total Duration (Corrected Formula)"%(B),position=(0.5,1.02),fontsize=14)
    ax.set_title("Transit Total Duration"%(B),position=(0.5,1.02),fontsize=14)

    cbar=fig.colorbar(c)
    cbar.ax.set_ylabel(r"$\Delta t_{\rm T}/t_{\rm T}$ (%)",fontsize=14)
    yts=cbar.ax.get_yticks()
    yl=[]
    levels=[]
    for yt in yts:
        yv=(dmin+yt*(dmax-dmin))
        yl+=["%.1f"%(yv*1E2)]
        levels+=[yv]
    cbar.ax.set_yticklabels(yl)
    ax.contour(IS,TS*RAD,transpose(DIST),
               levels=levels,
               colors=['k'],linestyles=[':'])

    plotPlanets(ax,RingedC,
                xmin=cieffmin,
                scalex=(cieffmax-cieffmin),
                ymin=teffmin*RAD,
                scaley=(teffmax-teffmin)*RAD)

    ax.set_xlim(cieffmin,cieffmax)
    ax.set_ylim(teffmin*RAD,teffmax*RAD)

    fig.savefig("figures/AccuracyTimes-TT.png")
    
    #////////////////////////////////////////
    #PLOT DISTF
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,6))
    ax=fig.gca()

    dmin=DISF.min()
    dmax=DISF.max()

    levels=linspace(dmin,dmax,100)
    c=ax.contourf(IS,TS*RAD,transpose(DISF),
                  levels=levels,cmap=cmap)

    ax.set_xlabel(r"$\cos\,i$",fontsize=20)
    ax.set_ylabel(r"$\theta$",fontsize=20)
    ax.set_title("Duration of Full Transit"%(B),position=(0.5,1.02),
                 fontsize=14)

    cbar=fig.colorbar(c)
    cbar.ax.set_ylabel(r"$\Delta t_{\rm F}/t_{\rm F}$ (%)",fontsize=14)
    yts=cbar.ax.get_yticks()
    yl=[]
    levels=[]
    for yt in yts:
        yv=(dmin+yt*(dmax-dmin))
        yl+=["%.1f"%((dmin+yt*(dmax-dmin))*1E2)]
        levels+=[yv]
    cbar.ax.set_yticklabels(yl)
    ax.contour(IS,TS*RAD,transpose(DISF),
               levels=levels,
               colors=['k'],linestyles=[':'])

    plotPlanets(ax,RingedC,
                xmin=cieffmin,
                scalex=(cieffmax-cieffmin),
                ymin=teffmin*RAD,
                scaley=(teffmax-teffmin)*RAD)
    ax.set_xlim(cieffmin,cieffmax)
    ax.set_ylim(teffmin*RAD,teffmax*RAD)

    fig.savefig("figures/AccuracyTimes-TF.png")

    #////////////////////////////////////////
    #PLOT DISTP
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,6))
    ax=fig.gca()

    dmin=DISTP.min()
    dmax=DISTP.max()

    levels=linspace(dmin,dmax,100)
    c=ax.contourf(IS,TS*RAD,transpose(DISTP),
                  levels=levels,cmap=cmap)

    ax.set_xlabel(r"$\cos\,i$",fontsize=20)
    ax.set_ylabel(r"$\theta$",fontsize=20)
    ax.set_title("Duration of Full Transit"%(B),position=(0.5,1.02),
                 fontsize=14)

    cbar=fig.colorbar(c)
    cbar.ax.set_ylabel(r"$|t_{\rm T}-t_{\rm T,P}|/t_{\rm T}$ (%)",fontsize=14)
    yts=cbar.ax.get_yticks()
    yl=[]
    levels=[]
    for yt in yts:
        yv=(dmin+yt*(dmax-dmin))
        yl+=["%.1f"%((dmin+yt*(dmax-dmin))*1E2)]
        levels+=[yv]
    cbar.ax.set_yticklabels(yl)
    ax.contour(IS,TS*RAD,transpose(DISTP),
               levels=levels,
               colors=['k'],linestyles=[':'])

    plotPlanets(ax,RingedC,
                xmin=cieffmin,
                scalex=(cieffmax-cieffmin),
                ymin=teffmin*RAD,
                scaley=(teffmax-teffmin)*RAD)
    ax.set_xlim(cieffmin,cieffmax)
    ax.set_ylim(teffmin*RAD,teffmax*RAD)

    fig.savefig("figures/AccuracyTimes-TP.png")

def contourTransitDepths():

    qload=False
    #qload=True
    print BARL,"Testing Analytical Times",RBAR

    #////////////////////////////////////////
    #SYSTEM
    #////////////////////////////////////////
    RingedC=copyObject(Ringed)
    ap=Ringed.ap/Ringed.Rstar
    P=Ringed.Porb/HOUR
    ip=Ringed.iorb
    B=RingedC.Borb
    Rp=RingedC.Rp

    #////////////////////////////////////////
    #MAKE A MAP
    #////////////////////////////////////////
    cieffmin=0.01
    cieffmax=1.0
    Ncieffs=30
    teffmin=0.0*DEG
    teffmax=90.0*DEG
    Nteffs=30

    if not qload:
        #NOT-RINGED
        Ap=pi*Rp**2

        cieffs=linspace(cieffmin,cieffmax,Ncieffs)
        teffs=linspace(teffmin,teffmax,Nteffs)

        IS,TS=meshgrid(cieffs,teffs)
        TDEPTH=zeros_like(IS)   
        TDEPTHN=zeros_like(IS)   
 
        ii=0
        for ci in cieffs:
            i=arccos(ci)
            print "Testing ieff = ",i*RAD
            jj=0
            for t in teffs:
                RingedC.ieff=i
                RingedC.Ringext.b=RingedC.Ringext.a*cos(i)
                RingedC.Ringext.cost=cos(t)
                RingedC.Ringext.sint=sin(t)

                RingedC.Ringint.b=RingedC.Ringint.a*cos(i)
                RingedC.Ringint.cost=cos(t)
                RingedC.Ringint.sint=sin(t)

                #NUMERICAL
                An=ringedPlanetArea(RingedC)

                #ANALYTICAL
                Aa=analyticalTransitArea(RingedC.Rp,RingedC.block,
                                         RingedC.fi,RingedC.fe,
                                         i)
                
                TDEPTH[ii,jj]=sqrt(Aa/Ap)
                TDEPTHN[ii,jj]=sqrt(An/Ap)

                jj+=1
            ii+=1

        savetxt("ISD.dat",IS)
        savetxt("TSD.dat",TS)
        savetxt("TDEPTH.dat",TDEPTH)
        savetxt("TDEPTHN.dat",TDEPTHN)

    IS=loadtxt("ISD.dat")
    TS=loadtxt("TSD.dat")
    TDEPTH=loadtxt("TDEPTH.dat")
    TDEPTHN=loadtxt("TDEPTHN.dat")
    
    cmap=plt.get_cmap("rainbow")
    #////////////////////////////////////////
    #CONTOUR
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,6))
    ax=fig.gca()

    dmin=TDEPTH.min()
    dmax=TDEPTH.max()

    levels=linspace(dmin,dmax,100)
    c=ax.contourf(IS,TS*RAD,transpose(TDEPTH),
                  levels=levels,cmap=cmap)

    ax.set_xlabel(r"$\cos\,i$",fontsize=20)
    ax.set_ylabel(r"$\theta$",fontsize=20)
    ax.set_title("Duration of Full Transit"%(B),position=(0.5,1.02),
                 fontsize=14)

    cbar=fig.colorbar(c)
    cbar.ax.set_ylabel(r"$R_{\rm p,obs}/R_{\rm p}$",fontsize=14)
    yts=cbar.ax.get_yticks()
    yl=[]
    levels=[]
    for yt in yts:
        yv=(dmin+yt*(dmax-dmin))
        yl+=["%.1f"%((dmin+yt*(dmax-dmin))*1E0)]
        levels+=[yv]
    cbar.ax.set_yticklabels(yl)
    ax.contour(IS,TS*RAD,transpose(TDEPTH),
               levels=levels,
               colors=['k'],linestyles=[':'])

    plotPlanets(ax,RingedC,
                xmin=cieffmin,
                scalex=(cieffmax-cieffmin),
                ymin=teffmin*RAD,
                scaley=(teffmax-teffmin)*RAD)
    ax.set_xlim(cieffmin,cieffmax)
    ax.set_ylim(teffmin*RAD,teffmax*RAD)

    fig.savefig("figures/TransitDepthContour.png")

def curveTransitDepths():

    print BARL,"Curve Transit Depth",RBAR

    #////////////////////////////////////////
    #SYSTEM
    #////////////////////////////////////////
    RingedC=copyObject(Ringed)
    ap=Ringed.ap/Ringed.Rstar
    P=Ringed.Porb/HOUR
    ip=Ringed.iorb
    B=RingedC.Borb
    Rp=RingedC.Rp
    Ap=pi*Rp**2

    #////////////////////////////////////////
    #MAKE A MAP
    #////////////////////////////////////////
    cieffmin=0.01
    cieffmax=1.0
    Ncieffs=100
    cieffs=linspace(cieffmin,cieffmax,Ncieffs)
    t=0.0
    rs=[]
    for ci in cieffs:
        i=arccos(ci)
        RingedC.ieff=i
        RingedC.Ringext.b=RingedC.Ringext.a*cos(i)
        RingedC.Ringext.cost=cos(t)
        RingedC.Ringext.sint=sin(t)

        RingedC.Ringint.b=RingedC.Ringint.a*cos(i)
        RingedC.Ringint.cost=cos(t)
        RingedC.Ringint.sint=sin(t)
        
        #NUMERICAL
        An=ringedPlanetArea(RingedC)
        
        #ANALYTICAL
        Aa=analyticalTransitArea(RingedC.Rp,RingedC.block,
                                 RingedC.fi,RingedC.fe,
                                 i)
        rs+=[sqrt(Aa/Ap)]
       
    #////////////////////////////////////////
    #CURVE
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    """
    tau=0.4
    Aas=[analyticalTransitArea(RingedC.Rp,blockFactor(tau,arccos(cieff)),
                               RingedC.fi,RingedC.fe,
                               arccos(cieff)) for cieff in cieffs]

    rs=sqrt(array(Aas)/Ap)
    ax.plot(cieffs,rs,'-',label=r"$\tau=%.1f$"%tau)

    tau=1.0
    Aas=[analyticalTransitArea(RingedC.Rp,blockFactor(tau,arccos(cieff)),
                               RingedC.fi,RingedC.fe,
                               arccos(cieff)) for cieff in cieffs]
    rs=sqrt(array(Aas)/Ap)
    ax.plot(cieffs,rs,'-',label=r"$\tau=%.1f$"%tau)

    tau=4.0
    Aas=[analyticalTransitArea(RingedC.Rp,blockFactor(tau,arccos(cieff)),
                               RingedC.fi,RingedC.fe,
                               arccos(cieff)) for cieff in cieffs]
    rs=sqrt(array(Aas)/Ap)
    ax.plot(cieffs,rs,'-',label=r"$\tau=%.1f$"%tau)
    """
    Ntau=60
    i=0
    cmap=cm.spectral
    for tau in logspace(log10(0.4),log10(4.0),Ntau):
        Aas=[analyticalTransitArea(RingedC.Rp,blockFactor(tau,arccos(cieff)),
                                   RingedC.fi,RingedC.fe,
                                   arccos(cieff)) for cieff in cieffs]
        rs=sqrt(array(Aas)/Ap)
        ax.plot(cieffs,rs,'-',label=r"$\tau=%.1f$"%tau,
                color=cmap((1.0*i)/Ntau))

        if (i%10)==0:
            ax.text(1.0,rs[-1],r"$\tau$=%.1f"%tau,
                    horizontalalignment='left',
                    verticalalignment='center',fontsize=10)
        i+=1
    ax.text(1.0,rs[-1],r"$\tau$=%.1f"%tau,
            horizontalalignment='left',
            verticalalignment='center',fontsize=10)

    ax.set_xlabel(r"$\cos\,i$",fontsize=20)
    ax.set_ylabel(r"$p_{\rm obs}/p$",fontsize=20)

    
    psize=0.02
    plotPlanets(ax,RingedC,Nx=6,Ny=2,ymin=0,yoffset=0.1,
                fh=psize/RingedC.Rp,fv=psize/RingedC.Rp)
    ax.set_xlim(cieffmin,cieffmax)
    ax.set_ylim(1.0,sqrt(RingedC.Ringext.a**2-RingedC.Ringint.a**2+Rp**2)/Rp)

    yt=ax.get_yticks()
    for y in yt[1:-1]:
        delta=y**2
        ax.text(0.03,y,"%.1f"%delta,
                horizontalalignment='left',
                verticalalignment='center')
    ax.text(0.07,0.5,r"$\delta_{\rm Ringed}/\delta_{\rm Not-ringed}$",
            horizontalalignment='left',
            verticalalignment='center',
            rotation=90,
            transform=ax.transAxes,fontsize=20,
            bbox=dict(fc='w',ec='none'))

    #ax.legend(loc="best")
    ax.grid(which="both")

    fig.savefig("figures/TransitDepthCurve.png")

def transitDepthPosterior():

    verbose=True
    #verbose=False
    try:
        qcalc=int(argv[1])
    except:
        qcalc=1
    
    #########################################
    #INPUT PARAMETERS
    #########################################
    """
    Transit depth input variables:
    Rp/Rstar, tau, fe, fi, i
    
    Simplification:

    fi can be fixed since transit depth depends on fe^2-fi^2. 

    Thus there are pair of values (fe,fi) and (fe',fi') such that:
    fe^2-fi^2 = fe'^2-fi'^2

    In summary: Rp, tau, fe, i

    Format of parameter array: [min, max, nominal, variable?]
    """
    S=System

    Nplanets=2000
    Nsamples=5

    Nbins=30

    Ntotal=Nplanets*Nsamples

    PARAMETERS["iorb"][DEF]=0.0
    PARAMETERS["ir"][STAT]=VAR
    PARAMETERS["phir"][STAT]=VAR

    PROPERTIES['p'][STAT]=0
    PROPERTIES['PR'][STAT]=0
    PROPERTIES['logPR'][STAT]=0
    
    #########################################
    #SAMPLE GENERATION
    #########################################
    if qcalc:

        i=1
        header=""
        header+="%-16s"%("#+0:id")
        for parkey in PARKEYS:
            parameter=PARAMETERS[parkey]
            scale=parameter[SCAL]
            if parameter[STAT]:
                header+="%-16s"%("+%d:%s(%s)"%(i,parkey,scale))
                i+=1

        for propkey in PROPKEYS:
            properti=PROPERTIES[propkey]
            scale=properti[SCAL]
            if properti[STAT]:
                header+="%-16s"%("+%d:%s(%s)"%(i,propkey,scale))
                i+=1

        header+="\n"
        npar=i-1

        line=[]
        for i in xrange(Nplanets*Nsamples):
            if (i%Nplanets)==0:
                print "Sample %d: %d planets generated..."%(i/Nplanets,i)
            data=[i]

            #========================================
            #RANDOM INPUT PARAMETERS
            #========================================
            for parkey in PARKEYS:
                parameter=PARAMETERS[parkey]
                func=parameter[FUNC]
                #GENERATE VALUE
                if parameter[STAT]:
                    val=func(randomVal(parameter[MIN],parameter[MAX]))
                    exec("data+=[val/%s]"%parameter[SCAL])
                else:
                    val=func(parameter[DEF])
                    
                #PREPARE SYSTEM
                exec("S.%s=val"%parkey)
                if verbose and 0:
                    print "Parameter %s:"%parkey
                    print TAB,"Default value = %e"%(func(parameter[DEF]))
                    print TAB,"Range = %e-%e"%(func(parameter[MIN]),
                                               func(parameter[MAX]))
                    print TAB,"Variable? = %d"%(parameter[STAT])
                    print TAB,"Adopted value = %e"%(val)
                    
            #========================================
            #TWEAK PARAMETERS
            #========================================
            """
            imax=arctan(S.ap/(S.Rstar-2*S.fe*S.Rplanet))
            cimax=cos(imax)
            S.iorb=arccos(randomVal(-cimax,cimax))
            data[0]=S.iorb/DEG
            """
            """
            S.ir=randomFisher(kappa=8,nsample=1)[0]
            if S.ir>pi/2:S.ir=pi/2-S.ir
            data[1]=S.ir*RAD
            #"""

            #========================================
            #UPDATE SYSTEM
            #========================================
            derivedSystemProperties(S)
            updatePlanetRings(S)

            #========================================
            #COMPUTE
            #========================================
            Ap=pi*S.Rp**2
            Aa=analyticalTransitAreaSystem(S)
            S.r=sqrt(Aa/Ap)

            S.p=S.PR=S.logPR=0
            #========================================
            #SAVE DERIVATIVE QUANTITIES
            #========================================
            for propkey in PROPKEYS:
                properti=PROPERTIES[propkey]
                func=properti[FUNC]
                #GENERATE VALUE
                if properti[STAT]:
                    exec("val=S.%s"%propkey)
                    exec("data+=[val/%s]"%properti[SCAL])

            line+=[data]
            
        savetxtheader("posterior-TransitDepth.dat",
                      header,line,fmt="%+.8e")
    
    data=loadtxt("posterior-TransitDepth.dat")
    rs=data[:,4]

    #########################################
    #STATISTICS
    #########################################
    xs,hs,dhs=histPosterior(rs,Nsamples,nbins=Nbins,
                               normed=True)

    #########################################
    #HISTOGRAM
    #########################################
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()
    error=True

    xms=histPlot(ax,xs,hs,dhs,error=error,color='b',alpha=0.1)
    hms=softArraySG(hs,frac=2,nP=3)
    histog=transpose(vstack((xms,hms)))
    savetxt("p-Posterior.dat",histog)

    ax.plot(xms,hms,'b-',linewidth=2,zorder=10)
    ax.set_xlabel(r"$p_{\rm obs}/p$",fontsize=20)
    ax.set_ylabel("Probability Density",fontsize=18)
    """
    ax.set_title(r"Posterior Distribution of Observed Radius",
                 position=(0.5,1.02),fontsize=16)
                 """
    #ax.set_yticks([])
    #ax.set_yticklabels([])
    yt=ax.get_yticks()
    yls=[]
    for y in yt:
        if y==0:yls+=["0"]
        else:yls+=[""]
    ax.set_yticklabels(yls)

    histogc=loadtxt("p-posterior-Concentrated.dat")
    datac=loadtxt("posterior-TransitDepth-Concentrated.dat")
    ax.plot(histogc[:,0],histogc[:,1],'r--',linewidth=2,zorder=10)

    #########################################
    #PLOT ORIENTATIONS
    #########################################
    R=0.11
    D=(R+0.05)
    X=1.0-D
    X=0.5
    Y=1-D
    tilt=20.0*DEG
    plotEllipse(ax,Figure(AR(X,Y),R,R,1,0,''),
                color='k',linestyle='-',
                transform=ax.transAxes)
    plotEllipse(ax,Figure(AR(X,Y),R,R*sin(tilt),1,0,''),
                color='k',linestyle='--',
                transform=ax.transAxes)
    ax.text(X,Y+R+0.01,"Ring Axis Orientation",fontsize=12,
            horizontalalignment='center',verticalalignment='bottom',
            transform=ax.transAxes)

    msize=10
    style=dict(markersize=2,markeredgecolor='none',
               color='b',marker='o',
               transform=ax.transAxes)

    nobjs=len(data)
    ies=arange(nobjs)
    for i in ies[::int(nobjs/200)]:
        theta=data[i,1]*DEG
        phi=data[i,2]*DEG
        x=sin(theta)*cos(phi)
        y=cos(theta)
        z=-sin(theta)*sin(phi)
        xp=x
        yp=y*cos(tilt)+z*sin(tilt)
        ax.plot([R*xp+X],[R*yp+Y],**style)
        ax.plot([-R*xp+X],[-R*yp+Y],**style)

    nobjs=len(datac)
    ies=arange(nobjs)
    style['color']='r'
    for i in ies[::int(nobjs/200)]:
        theta=datac[i,1]*DEG
        phi=datac[i,2]*DEG
        x=sin(theta)*cos(phi)
        y=cos(theta)
        z=-sin(theta)*sin(phi)
        xp=x
        yp=y*cos(tilt)+z*sin(tilt)
        ax.plot([R*xp+X],[R*yp+Y],**style)
        ax.plot([-R*xp+X],[-R*yp+Y],**style)

    ax.set_xlim((1.0,xs[-1]))
    hmin,hmax=ax.get_ylim()
    xmin,xmax=ax.get_xlim()
    hmax=max(histog.max(),histogc.max(),(hs+dhs).max())

    ax.set_ylim((0.0,hmax))
    fig.savefig("figures/posterior-TransitDepth.png")

def contourPhotoRing():

    S=System
    try:
        qcalc=int(argv[1])
    except:
        qcalc=1

    print BARL,"Calculating Photo-ring effect contours",RBAR

    #////////////////////////////////////////
    #SYSTEM
    #////////////////////////////////////////
    RingedC=copyObject(Ringed)
    ap=Ringed.ap/Ringed.Rstar
    P=Ringed.Porb/HOUR
    ip=Ringed.iorb
    B=RingedC.Borb
    Rp=RingedC.Rp
    rho_true=S.Mstar/(4*pi/3*S.Rstar**3)

    #////////////////////////////////////////
    #MAKE A MAP
    #////////////////////////////////////////
    cieffmin=0.01
    cieffmax=1.0
    Ncieffs=30
    teffmin=0.0*DEG
    teffmax=90.0*DEG
    Nteffs=30
    qcorrected=False

    if qcalc:
        cieffs=linspace(cieffmin,cieffmax,Ncieffs)
        teffs=linspace(teffmin,teffmax,Nteffs)
        #cieffs=[0.5];teffs=[45*DEG]
        
        IS,TS=meshgrid(cieffs,teffs)
        PR=zeros_like(IS)   
        PRN=zeros_like(IS)   
 
        ii=0
        for ci in cieffs:
            i=arccos(ci)
            print "Testing ieff (cos ieff) = ",i*RAD,ci
            jj=0
            for t in teffs:
                print TAB,"Testing teff = ",t*RAD
                RingedC.ieff=i

                RingedC.block=blockFactor(RingedC.tau,i)

                RingedC.Ringext.b=RingedC.Ringext.a*cos(i)
                RingedC.Ringext.cost=cos(t)
                RingedC.Ringext.sint=sin(t)

                RingedC.Ringint.b=RingedC.Ringint.a*cos(i)
                RingedC.Ringint.cost=cos(t)
                RingedC.Ringint.sint=sin(t)
                
                B=RingedC.Borb

                #NUMERICAL
                tcsp=contactTimes(RingedC)
                tT=(tcsp[-1]-tcsp[1])/HOUR
                tF=(tcsp[-2]-tcsp[2])/HOUR

                p=ringedPlanetArea(RingedC)/pi
                pRn=rhoObserved_Seager(p,
                                       RingedC.Rstar,
                                       tT,tF,
                                       P)/rho_true

                #ANALYTICAL
                xpa1=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                                     direction=-1,sign=-1,
                                     qcorrected=qcorrected)
                xpa2=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                                     direction=-1,sign=+1,
                                     qcorrected=qcorrected)
                xpa3=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                                     direction=+1,sign=-1,
                                     qcorrected=qcorrected)
                xpa4=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                                     direction=+1,sign=+1,
                                     qcorrected=qcorrected)
                tT=P*arcsin((xpa4-xpa1)/(ap*sin(ip)))/(2*pi)
                tF=P*arcsin((xpa3-xpa2)/(ap*sin(ip)))/(2*pi)

                p=analyticalTransitArea(RingedC.Rp,RingedC.block,
                                        RingedC.fi,RingedC.fe,
                                        i)/pi
                pR=rhoObserved_Seager(p,
                                      RingedC.Rstar,
                                      tT,tF,
                                      P)/rho_true

                print 2*TAB,"PR (Numerical) = %.6e, PR (Analytical) = %.6e"%(log10(pRn),log10(pR))

                PR[ii,jj]=log10(pR)
                PRN[ii,jj]=log10(pRn)

                jj+=1
            ii+=1

        savetxt("ISP.dat",IS)
        savetxt("TSP.dat",TS)
        savetxt("PR.dat",PR)
        savetxt("PRN.dat",PRN)

    IS=loadtxt("ISP.dat")
    TS=loadtxt("TSP.dat")
    PR=loadtxt("PR.dat")
    PRN=loadtxt("PRN.dat")
    data=loadtxt("posterior-PhotoRing-Concentrated.dat") #EXPERIMENT
    
    cmap=plt.get_cmap("rainbow")

    #////////////////////////////////////////
    #CONTOUR
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,6))
    ax=fig.gca()

    dmin=PRN.min()
    dmax=PRN.max()
    
    levels=linspace(dmin,dmax,100)
    c=ax.contourf(IS,TS*RAD,transpose(PRN),
                  levels=levels,cmap=cmap)

    freq=50
    ax.plot(cos(data[::freq,4]*DEG),abs(data[::freq,7]),
            '+',color=cm.gray(0.2),zorder=5) #EXP

    ax.set_xlabel(r"$\cos\,i$",fontsize=20)
    ax.set_ylabel(r"$\theta$",fontsize=20)
    ax.set_title("Photo-Ring Effect"%(B),position=(0.5,1.02),
                 fontsize=14)

    cbar=fig.colorbar(c)
    cbar.ax.set_ylabel(r"$\log_{10}(\rho_{\rm obs}/\rho_\star)$",fontsize=20)
    yts=cbar.ax.get_yticks()
    yl=[]
    levels=[]
    for yt in yts:
        yv=(dmin+yt*(dmax-dmin))
        yl+=["%.3f"%((dmin+yt*(dmax-dmin))*1E0)]
        levels+=[yv]
    cbar.ax.set_yticklabels(yl)
    c=ax.contour(IS,TS*RAD,transpose(PRN),
               levels=levels,
               colors=['k'],linestyles=[':'])
    ax.clabel(c,inline=1,fontsize=10,zorder=50)
    levelsC=levels

    ax.contour(IS,TS*RAD,transpose(PRN),
               levels=[0.0],
               colors=['k'],linestyles=['-'],linewidths=['2'])
    
    plotPlanets(ax,RingedC,
                xmin=cieffmin,
                scalex=(cieffmax-cieffmin),
                ymin=teffmin*RAD,
                scaley=(teffmax-teffmin)*RAD)
    ax.set_xlim(cieffmin,cieffmax)
    ax.set_ylim(teffmin*RAD,teffmax*RAD)

    fig.savefig("figures/PhotoRingContour.png")

    #////////////////////////////////////////
    #CONTOUR
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,6))
    ax=fig.gca()

    dmin=PR.min()
    dmax=PR.max()
    
    levels=linspace(dmin,dmax,100)
    c=ax.contourf(IS,TS*RAD,transpose(PR),
                  levels=levels,cmap=cmap)

    ax.set_xlabel(r"$\cos\,i$",fontsize=20)
    ax.set_ylabel(r"$\theta$",fontsize=20)
    ax.set_title("Photo-Ring Effect"%(B),position=(0.5,1.02),
                 fontsize=14)

    cbar=fig.colorbar(c)
    cbar.ax.set_ylabel(r"$\log(\rho_{\rm obs}/\rho_\star)$",fontsize=14)
    yts=cbar.ax.get_yticks()
    yl=[]
    levels=[]
    for yt in yts:
        yv=(dmin+yt*(dmax-dmin))
        yl+=["%.3f"%((dmin+yt*(dmax-dmin))*1E0)]
        levels+=[yv]
    cbar.ax.set_yticklabels(yl)
    c=ax.contour(IS,TS*RAD,transpose(PR),
                 levels=levelsC,
                 colors=['k'],linestyles=[':'])
    ax.clabel(c,inline=1,fontsize=10)

    ax.contour(IS,TS*RAD,transpose(PR),
               levels=[0.0],
               colors=['k'],linestyles=['-'],linewidths=['2'])
    
    plotPlanets(ax,RingedC,
                xmin=cieffmin,
                scalex=(cieffmax-cieffmin),
                ymin=teffmin*RAD,
                scaley=(teffmax-teffmin)*RAD)
    ax.set_xlim(cieffmin,cieffmax)
    ax.set_ylim(teffmin*RAD,teffmax*RAD)

    fig.savefig("figures/PhotoRingContour-Analytical.png")

def photoRingPosterior():

    verbose=True
    #verbose=False
    try:
        qcalc=int(argv[1])
    except:
        qcalc=1
    
    #########################################
    #INPUT PARAMETERS
    #########################################
    """
    Transit depth input variables:
    Rp/Rstar, tau, fe, fi, i
    
    Simplification:

    fi can be fixed since transit depth depends on fe^2-fi^2. 

    Thus there are pair of values (fe,fi) and (fe',fi') such that:
    fe^2-fi^2 = fe'^2-fi'^2

    In summary: Rp, tau, fe, i

    Format of parameter array: [min, max, nominal, variable?]
    """
    S=System

    NplanetsPR=2000
    NsamplesPR=5

    Nbins=30

    Ntotal=NplanetsPR*NsamplesPR

    PARAMETERS["iorb"][DEF]=0.0
    PARAMETERS["ir"][STAT]=VAR
    PARAMETERS["phir"][STAT]=VAR
    PROPERTIES['r'][STAT]=SHOW

    #########################################
    #SAMPLE GENERATION
    #########################################
    Ap=pi*S.Rp**2
    if qcalc:

        i=1
        header=""
        header+="%-16s"%("#+0:id")
        for parkey in PARKEYS:
            parameter=PARAMETERS[parkey]
            scale=parameter[SCAL]
            if parameter[STAT]:
                header+="%-16s"%("+%d:%s(%s)"%(i,parkey,scale))
                i+=1

        for propkey in PROPKEYS:
            properti=PROPERTIES[propkey]
            scale=properti[SCAL]
            if properti[STAT]:
                header+="%-16s"%("+%d:%s(%s)"%(i,propkey,scale))
                i+=1

        header+="\n"
        npar=i-1

        line=[]
        for i in xrange(NplanetsPR*NsamplesPR):
            if (i%NplanetsPR)==0:
                print "Sample %d: %d planets generated..."%(i/NplanetsPR,i)
            data=[i]

            #========================================
            #RANDOM INPUT PARAMETERS
            #========================================
            for parkey in PARKEYS:
                parameter=PARAMETERS[parkey]
                func=parameter[FUNC]
                #GENERATE VALUE
                if parameter[STAT]:
                    val=func(randomVal(parameter[MIN],parameter[MAX]))
                    exec("data+=[val/%s]"%parameter[SCAL])
                else:
                    val=func(parameter[DEF])
                    
                #PREPARE SYSTEM
                exec("S.%s=val"%parkey)
                if verbose and 0:
                    print "Parameter %s:"%parkey
                    print TAB,"Default value = %e"%(func(parameter[DEF]))
                    print TAB,"Range = %e-%e"%(func(parameter[MIN]),
                                               func(parameter[MAX]))
                    print TAB,"Variable? = %d"%(parameter[STAT])
                    print TAB,"Adopted value = %e"%(val)
                    
            #========================================
            #TWEAK PARAMETERS
            #========================================
            """
            imax=arctan(S.ap/(S.Rstar-2*S.fe*S.Rplanet))
            cimax=cos(imax)
            S.iorb=arccos(randomVal(-cimax,cimax))
            data[0]=S.iorb/DEG
            """
            #"""
            S.ir=randomFisher(kappa=8,nsample=1)[0]
            if S.ir>pi/2:S.ir=pi/2-S.ir
            data[1]=S.ir*RAD
            #"""
            
            #========================================
            #UPDATE SYSTEM
            #========================================
            derivedSystemProperties(S)
            updatePlanetRings(S)

            #========================================
            #COMPUTE
            #========================================
            rho_true=S.Mstar/(4*pi/3*S.Rstar**3)

            tcsp=contactTimes(S)
            tT=(tcsp[-1]-tcsp[1])/HOUR
            tF=(tcsp[-2]-tcsp[2])/HOUR
            Aa=S.p=analyticalTransitArea(S.Rp,S.block,
                                         S.fi,S.fe,
                                         S.ieff)/pi
            S.PR=rhoObserved_Seager(S.p,S.Rstar,
                                    tT,tF,S.Porb/HOUR)/rho_true
            
            S.logPR=log10(S.PR)
            S.r=sqrt(Aa/Ap)

            #========================================
            #SAVE DERIVATIVE QUANTITIES
            #========================================
            for propkey in PROPKEYS:
                properti=PROPERTIES[propkey]
                func=properti[FUNC]
                #GENERATE VALUE
                if properti[STAT]:
                    exec("val=S.%s"%propkey)
                    exec("data+=[val/%s]"%properti[SCAL])

            line+=[data]
            
        savetxtheader("posterior-PhotoRing.dat",
                      header,line,fmt="%+.8e")
    
    data=loadtxt("posterior-PhotoRing.dat")
    PRs=data[:,5]
    #exit(0)

    #########################################
    #STATISTICS
    #########################################
    xs,hs,dhs=histPosterior(PRs,NsamplesPR,nbins=Nbins,
                               normed=True)

    #########################################
    #HISTOGRAM
    #########################################
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()
    error=True

    xms=histPlot(ax,xs,hs,dhs,error=error,color='b',alpha=0.1)
    hms=hs
    hms=softArraySG(hs,frac=6,nP=2)
    histog=transpose(vstack((xms,hms)))
    savetxt("pr-Posterior.dat",histog)
    
    hfuncs=interpolant(xms,hms,kind='cubic')
    xvec=linspace(xms[0],xms[-1],1000)
    hvec=hfuncs(xvec)
    ax.plot(xvec,hvec,'b-',linewidth=2,zorder=10)

    ax.set_xlabel(r"$\log_{10}(\rho_{\star,\mathrm{transit}}/\rho_{\star,\mathrm{indep}})$",fontsize=20)
    ax.set_ylabel("Probability Density",fontsize=18)
    """
    ax.set_title(r"Posterior Distribution of Photo-Ring Effect",
                 position=(0.5,1.02),fontsize=16)
    """
    ax.axvline(0.0,color='k',linestyle='--',linewidth=2)
    #ax.set_yticks([])
    #ax.set_yticklabels([])
    yt=ax.get_yticks()
    yls=[]
    for y in yt:
        if y==0:yls+=["0"]
        else:yls+=[""]
    ax.set_yticklabels(yls)

    #"""
    histogc=loadtxt("pr-posterior-Concentrated.dat")
    datac=loadtxt("posterior-PhotoRing-Concentrated.dat")
    ax.plot(histogc[:,0],histogc[:,1],'r--',linewidth=2,zorder=10)
    #"""
    #########################################
    #PLOT ORIENTATIONS
    #########################################
    R=0.11
    D=(R+0.05)
    X=D
    #X=0.5
    Y=1-D
    tilt=20.0*DEG
    plotEllipse(ax,Figure(AR(X,Y),R,R,1,0,''),
                color='k',linestyle='-',
                transform=ax.transAxes)
    plotEllipse(ax,Figure(AR(X,Y),R,R*sin(tilt),1,0,''),
                color='k',linestyle='--',
                transform=ax.transAxes)
    ax.text(X,Y+R+0.01,"Ring Axis Orientation",fontsize=12,
            horizontalalignment='center',verticalalignment='bottom',
            transform=ax.transAxes)

    msize=10
    style=dict(markersize=2,markeredgecolor='none',
               color='b',marker='o',
               transform=ax.transAxes)

    #"""
    nobjs=len(data)
    ies=arange(nobjs)
    for i in ies[::int(nobjs/200)]:
        theta=data[i,1]*DEG
        phi=data[i,2]*DEG
        x=sin(theta)*cos(phi)
        y=cos(theta)
        z=-sin(theta)*sin(phi)
        xp=x
        yp=y*cos(tilt)+z*sin(tilt)
        ax.plot([R*xp+X],[R*yp+Y],**style)
        ax.plot([-R*xp+X],[-R*yp+Y],**style)

    nobjs=len(datac)
    ies=arange(nobjs)
    style['color']='r'
    for i in ies[::int(nobjs/200)]:
        theta=datac[i,1]*DEG
        phi=datac[i,2]*DEG
        x=sin(theta)*cos(phi)
        y=cos(theta)
        z=-sin(theta)*sin(phi)
        xp=x
        yp=y*cos(tilt)+z*sin(tilt)
        ax.plot([R*xp+X],[R*yp+Y],**style)
        ax.plot([-R*xp+X],[-R*yp+Y],**style)
    #"""

    PRmin=10**xs[0]
    PRmax=10**xs[-1]
    PRs=concatenate((arange(0.1,1.0,0.1),
                     arange(1.1,2.0,0.1)
                     ))
    """
    xt=[];xl=[]
    for PR in PRs:
        xt+=[log10(PR)]
        xl+=["%.1f"%PR]
    ax.set_xticks(xt)
    ax.set_xticklabels(xl)
    """
    ax.set_xlim((xs[0],xs[-1]))
    xmin,xmax=ax.get_xlim()
    hmin,hmax=ax.get_ylim()
    xt=ax.get_xticks()
    for x in xt:
        if x<xmin or x>xmax:continue
        PR=1/10**x
        ax.text(x,hmax,"%.1f"%PR,
                horizontalalignment='center',verticalalignment='bottom',
                transform=offSet(0,5))
        #"""
    ax.text(0.5,1.05,r"$\rho_{\star,\mathrm{indep}}/\rho_{\star,\mathrm{transit}}$",
            horizontalalignment='center',verticalalignment='bottom',
            fontsize=20,transform=ax.transAxes)
            #"""
    
    ax.set_ylim((0.0,hmax))
    fig.savefig("figures/posterior-PhotoRing-relabeled.pdf")

    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()
    
    data=loadtxt("posterior-PhotoRing-Concentrated.dat")
    PRs=data[:,5]
    ieff=data[:,4]
    teff=data[:,7]
    rs=sqrt(data[:,6]/S.Rp**2)

    ax.plot(PRs,rs,'bo',markersize=2,markeredgecolor='none')

    cond=abs(teff)<10.0
    ax.plot(PRs[cond],rs[cond],'ro',markersize=2,markeredgecolor='none')

    ax.set_xlabel(r'$\log_{10}(\rho_{\rm obs}/\rho_\star)$',fontsize=16)
    ax.set_ylabel(r'$p_{\rm obs}/p$',fontsize=16)

    fig.savefig("figures/scatter-PR-TransitDepth.png")

def photoRingPosteriorTiming():

    verbose=True
    #verbose=False
    try:
        qcalc=int(argv[1])
    except:
        qcalc=1
    
    #########################################
    #INPUT PARAMETERS
    #########################################
    """
    Transit depth input variables:
    Rp/Rstar, tau, fe, fi, i
    
    Simplification:

    fi can be fixed since transit depth depends on fe^2-fi^2. 

    Thus there are pair of values (fe,fi) and (fe',fi') such that:
    fe^2-fi^2 = fe'^2-fi'^2

    In summary: Rp, tau, fe, i

    Format of parameter array: [min, max, nominal, variable?]
    """
    S=System

    NplanetsPR=2000
    NsamplesPR=5

    Nbins=30

    Ntotal=NplanetsPR*NsamplesPR

    PARAMETERS["iorb"][DEF]=0.0
    PARAMETERS["ir"][STAT]=VAR
    PARAMETERS["phir"][STAT]=VAR
    PROPERTIES['r'][STAT]=SHOW

    #########################################
    #SAMPLE GENERATION
    #########################################
    Ap=pi*S.Rp**2
    if qcalc:

        i=1
        header=""
        header+="%-16s"%("#+0:id")
        for parkey in PARKEYS:
            parameter=PARAMETERS[parkey]
            scale=parameter[SCAL]
            if parameter[STAT]:
                header+="%-16s"%("+%d:%s(%s)"%(i,parkey,scale))
                i+=1

        for propkey in PROPKEYS:
            properti=PROPERTIES[propkey]
            scale=properti[SCAL]
            if properti[STAT]:
                header+="%-16s"%("+%d:%s(%s)"%(i,propkey,scale))
                i+=1

        header+="\n"
        npar=i-1

        line=[]
        for i in xrange(NplanetsPR*NsamplesPR):
            if (i%NplanetsPR)==0:
                print "Sample %d: %d planets generated..."%(i/NplanetsPR,i)
            data=[i]

            #========================================
            #RANDOM INPUT PARAMETERS
            #========================================
            for parkey in PARKEYS:
                parameter=PARAMETERS[parkey]
                func=parameter[FUNC]
                #GENERATE VALUE
                if parameter[STAT]:
                    val=func(randomVal(parameter[MIN],parameter[MAX]))
                    exec("data+=[val/%s]"%parameter[SCAL])
                else:
                    val=func(parameter[DEF])
                    
                #PREPARE SYSTEM
                exec("S.%s=val"%parkey)
                if verbose and 0:
                    print "Parameter %s:"%parkey
                    print TAB,"Default value = %e"%(func(parameter[DEF]))
                    print TAB,"Range = %e-%e"%(func(parameter[MIN]),
                                               func(parameter[MAX]))
                    print TAB,"Variable? = %d"%(parameter[STAT])
                    print TAB,"Adopted value = %e"%(val)
                    
            #========================================
            #TWEAK PARAMETERS
            #========================================
            """
            imax=arctan(S.ap/(S.Rstar-2*S.fe*S.Rplanet))
            cimax=cos(imax)
            S.iorb=arccos(randomVal(-cimax,cimax))
            data[0]=S.iorb/DEG
            """
            #"""
            S.ir=randomFisher(kappa=8,nsample=1)[0]
            if S.ir>pi/2:S.ir=pi/2-S.ir
            data[1]=S.ir*RAD
            #"""
            
            #========================================
            #UPDATE SYSTEM
            #========================================
            derivedSystemProperties(S)
            updatePlanetRings(S)
            #print S.ieff*RAD,S.teff*RAD,S.Borb

            #========================================
            #COMPUTE
            #========================================
            rho_true=S.Mstar/(4*pi/3*S.Rstar**3)

            #"""
            tcsp=contactTimes(S)
            tT=(tcsp[-1]-tcsp[1])/HOUR
            tF=(tcsp[-2]-tcsp[2])/HOUR
            #print "Numerical: ",tT,tF
            #"""

            """
            xpa1=transitPosition(S.Rp,S.fe,S.ieff,S.teff,S.Borb,
                         direction=-1,sign=-1)
            xpa2=transitPosition(S.Rp,S.fe,S.ieff,S.teff,S.Borb,
                                 direction=-1,sign=+1)
            xpa3=transitPosition(S.Rp,S.fe,S.ieff,S.teff,S.Borb,
                                 direction=+1,sign=-1)
            xpa4=transitPosition(S.Rp,S.fe,S.ieff,S.teff,S.Borb,
                                 direction=+1,sign=+1)
            tT=S.Porb*arcsin((xpa4-xpa1)/\
                                 (S.ap/S.Rstar*sin(S.iorb)))/(2*pi)/HOUR
            tF=S.Porb*arcsin((xpa3-xpa2)\
                                 /(S.ap/S.Rstar*sin(S.iorb)))/(2*pi)/HOUR
            #print "Analytical: ",tT,tF
            #"""
            #raw_input()

            Aa=S.p=analyticalTransitArea(S.Rp,S.block,
                                      S.fi,S.fe,
                                      S.ieff)/pi
            S.PR=rhoObserved_Seager(S.p,S.Rstar,
                                    tT,tF,S.Porb/HOUR)/rho_true
            
            S.logPR=log10(S.PR)
            S.r=sqrt(Aa/Ap)

            #========================================
            #SAVE DERIVATIVE QUANTITIES
            #========================================
            for propkey in PROPKEYS:
                properti=PROPERTIES[propkey]
                func=properti[FUNC]
                #GENERATE VALUE
                if properti[STAT]:
                    exec("val=S.%s"%propkey)
                    exec("data+=[val/%s]"%properti[SCAL])

            line+=[data]
            
        savetxtheader("posterior-PhotoRing-timing.dat",
                      header,line,fmt="%+.8e")
    
    data=loadtxt("posterior-PhotoRing-timing.dat")
    PRs=data[:,5]
    #return

    #########################################
    #STATISTICS
    #########################################
    xs,hs,dhs=histPosterior(PRs,NsamplesPR,nbins=Nbins,
                               normed=True)

    #########################################
    #HISTOGRAM
    #########################################
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()
    error=True

    xms=histPlot(ax,xs,hs,dhs,error=error,color='w',alpha=0.1)
    hms=hs
    hms=softArraySG(hs,frac=6,nP=2)
    histog=transpose(vstack((xms,hms)))
    savetxt("pr-Posterior.dat",histog)
    
    hfuncs=interpolant(xms,hms,kind='cubic')
    xvec=linspace(xms[0],xms[-1],1000)
    hvec=hfuncs(xvec)
    ax.plot(xvec,hvec,'b-',linewidth=2,zorder=10,label='Analytical')

    ax.set_xlabel(r"$\log_{10}(\rho_{\rm obs}/\rho_\star)$",fontsize=20)
    ax.set_ylabel("Probability Density",fontsize=18)
    """
    ax.set_title(r"Posterior Distribution of Photo-Ring Effect",
                 position=(0.5,1.02),fontsize=16)
    """
    ax.axvline(0.0,color='k',linestyle='--',linewidth=2)
    ax.set_yticks([])
    ax.set_yticklabels([])

    #"""
    histogc=loadtxt("pr-posterior-Concentrated-timing.dat")
    datac=loadtxt("posterior-PhotoRing-Concentrated-timing.dat")
    ax.plot(histogc[:,0],histogc[:,1],'r-',label='Numerical',linewidth=2,zorder=10)
    #"""
    #########################################
    #PLOT ORIENTATIONS
    #########################################
    """
    R=0.11
    D=(R+0.05)
    X=D
    #X=0.5
    Y=1-D
    tilt=20.0*DEG
    plotEllipse(ax,Figure(AR(X,Y),R,R,1,0,''),
                color='k',linestyle='-',
                transform=ax.transAxes)
    plotEllipse(ax,Figure(AR(X,Y),R,R*sin(tilt),1,0,''),
                color='k',linestyle='--',
                transform=ax.transAxes)
    ax.text(X,Y+R+0.01,"Ring Axis Orientation",fontsize=12,
            horizontalalignment='center',verticalalignment='bottom',
            transform=ax.transAxes)
    #"""
    msize=10
    style=dict(markersize=2,markeredgecolor='none',
               color='b',marker='o',
               transform=ax.transAxes)

    """
    nobjs=len(data)
    ies=arange(nobjs)
    for i in ies[::int(nobjs/200)]:
        theta=data[i,1]*DEG
        phi=data[i,2]*DEG
        x=sin(theta)*cos(phi)
        y=cos(theta)
        z=-sin(theta)*sin(phi)
        xp=x
        yp=y*cos(tilt)+z*sin(tilt)
        ax.plot([R*xp+X],[R*yp+Y],**style)
        ax.plot([-R*xp+X],[-R*yp+Y],**style)

    nobjs=len(datac)
    ies=arange(nobjs)
    style['color']='r'
    for i in ies[::int(nobjs/200)]:
        theta=datac[i,1]*DEG
        phi=datac[i,2]*DEG
        x=sin(theta)*cos(phi)
        y=cos(theta)
        z=-sin(theta)*sin(phi)
        xp=x
        yp=y*cos(tilt)+z*sin(tilt)
        ax.plot([R*xp+X],[R*yp+Y],**style)
        ax.plot([-R*xp+X],[-R*yp+Y],**style)
    #"""

    PRmin=10**xs[0]
    PRmax=10**xs[-1]
    PRs=concatenate((arange(0.1,1.0,0.1),
                     arange(1.1,2.0,0.1)
                     ))
    """
    xt=[];xl=[]
    for PR in PRs:
        xt+=[log10(PR)]
        xl+=["%.1f"%PR]
    ax.set_xticks(xt)
    ax.set_xticklabels(xl)
    """
    ax.set_xlim((xs[0],xs[-1]))
    xmin,xmax=ax.get_xlim()
    hmin,hmax=ax.get_ylim()
    xt=ax.get_xticks()
    for x in xt:
        if x<xmin or x>xmax:continue
        PR=1/10**x
        ax.text(x,hmax,"%.1f"%PR,
                horizontalalignment='center',verticalalignment='bottom',
                transform=offSet(0,5))
    ax.text(0.5,1.05,r"$\rho_\star/\rho_{\rm obs}$",
            horizontalalignment='center',verticalalignment='bottom',
            fontsize=20,transform=ax.transAxes)

    ax.legend(loc="upper right")
    ax.set_ylim((0.0,hmax))
    fig.savefig("figures/posterior-PhotoRing-timing.png")

def testPhotoRing():

    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    S=Ringed
    t=Ringed.teff
    i=Ringed.ieff

    print BARL,"Test Transit Duration",RBAR


    #========================================
    #FIX RINGED PROPERTIES BY HAND
    #========================================
    #MANUAL i,t
    i=80.0*DEG;t=60*DEG
    
    Ringed.ieff=i;Ringed.teff=t

    Ringed.Ringext.b=Ringed.Ringext.a*cos(i)
    Ringed.Ringext.cost=cos(t);Ringed.Ringext.sint=sin(t)

    Ringed.Ringint.b=Ringed.Ringint.a*cos(i)
    Ringed.Ringint.cost=cos(t);Ringed.Ringint.sint=sin(t)

    Ringed.block=blockFactor(Ringed.tau,i)

    print "Orientation parameters:"
    print TAB,"i = %.2f deg"%(i*RAD)
    print TAB,"t = %.2f deg"%(t*RAD)

    #========================================
    #PROPERTIES
    #========================================
    RingedC=copyObject(Ringed)
    ap=Ringed.ap/Ringed.Rstar
    P=Ringed.Porb/HOUR
    ip=Ringed.iorb
    rho_true=S.Mstar/(4*pi/3*S.Rstar**3)

    print "True density = ",rho_true
    
    #========================================
    #NOT RINGED DENSITY
    #========================================
    xp=sqrt((1+NotRinged.Rp)**2-NotRinged.Borb**2)
    xm=sqrt((1-NotRinged.Rp)**2-NotRinged.Borb**2)
    tT=P*arcsin(2*xp/(ap*sin(ip)))/(2*pi)
    tF=P*arcsin(2*xm/(ap*sin(ip)))/(2*pi)
    p=S.Rp**2

    print "t_T, t_F = ",tT,tF
    rho_obs=rhoObserved_Seager(p,S.Rstar,
                               tT,tF,S.Porb/HOUR)

    print "Observed density (not ringed, Seager) = ",rho_obs
    rho_obs=rhoObserved_Kipping(p,S.Rstar,
                               tT,tF,S.Porb/HOUR)
    print "Observed density (not ringed, Kipping) = ",rho_obs

    #========================================
    #FIX RINGEDC PROPERTIES BY HAND
    #========================================
    #MANUAL i,t
    
    #GOOD SPOT
    i=89.43*DEG;t=0.0*DEG

    print "Orientation parameters (MANUAL):"
    print TAB,"b = %.2f"%(RingedC.Borb)
    print TAB,"i = %.2f deg"%(i*RAD)
    print TAB,"t = %.2f deg"%(t*RAD)
    
    RingedC.Ringext.b=RingedC.Ringext.a*cos(i)
    RingedC.Ringint.b=RingedC.Ringint.a*cos(i)
    RingedC.Ringext.cost=cos(t);RingedC.Ringext.sint=sin(t)
    RingedC.Ringint.cost=cos(t);RingedC.Ringint.sint=sin(t)
    RingedC.block=blockFactor(RingedC.tau,i)

    #========================================
    #RINGED TRANSIT DURATION (NUMERICAL)
    #========================================
    print BARL,"NUMERICAL",RBAR

    p=ringedPlanetArea(RingedC)/pi

    tcsp=contactTimes(RingedC)
    tT=(tcsp[-1]-tcsp[1])/HOUR
    tF=(tcsp[-2]-tcsp[2])/HOUR

    rho_obs=rhoObserved_Seager(p,S.Rstar,
                               tT,tF,S.Porb/HOUR)
    pR=log10(rho_obs/rho_true)

    print "p = ",p
    print "t_T, t_F = ",tT,tF
    print "Observed density (ringed, numerical, Seager) = ",rho_obs
    print "Photo-ring effect (numerical) = ",pR

    #========================================
    #RINGED TRANSIT ANALYTICAL 
    #========================================
    print BARL,"ANALYTICAL",RBAR

    p=analyticalTransitArea(RingedC.Rp,RingedC.block,
                            RingedC.fi,RingedC.fe,
                            i)/pi
    
    a=RingedC.Ringext.a
    b=RingedC.Ringext.b
    B=RingedC.Borb

    xpa1=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                         direction=-1,sign=-1,qcorrected=True)
    xpa2=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                         direction=-1,sign=+1,qcorrected=True)
    xpa3=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                         direction=+1,sign=-1,qcorrected=True)
    xpa4=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                         direction=+1,sign=+1,qcorrected=True)

    taT=P*arcsin((xpa4-xpa1)/(ap*sin(ip)))/(2*pi)
    taF=P*arcsin((xpa3-xpa2)/(ap*sin(ip)))/(2*pi)

    rho_obsa=rhoObserved_Seager(p,S.Rstar,
                                taT,taF,S.Porb/HOUR)
    pRa=log10(rho_obsa/rho_true)

    print "p = ",p
    print "t_T, t_F = ",taT,taF
    print "Observed density (ringed, analytical, Seager) = ",rho_obsa
    print "Photo-ring effect (analytical) = ",pRa

    #========================================
    #ERRORS
    #========================================
    drho=abs(rho_obs-rho_obsa)/rho_true*100
    print "Error in observed density = %.2f%%"%(drho)
    dT=abs(tT-taT)/tT*100
    print "Error in total transit time = %.2f%%"%(dT)
    dF=abs(tF-taF)/tF*100
    print "Error in full transit = %.2f%%"%(dF)
    
    #========================================
    #PLOT
    #========================================
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    plotEllipse(ax,RingedC.Star,color='y')
    plotEllipse(ax,RingedC.Planet,color='b')
    plotEllipse(ax,RingedC.Ringext,color='k')
    plotEllipse(ax,RingedC.Ringint,color='r')
    
    rng=1.5
    Re=RingedC.Ringext.a
    Re=1.0
    xmin=RingedC.Planet.C[0]-rng*Re;
    xmax=RingedC.Planet.C[0]+rng*Re
    ymin=RingedC.Planet.C[1]-rng*Re;
    ymax=RingedC.Planet.C[1]+rng*Re
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((ymin,ymax))
    ax.grid()
    fig.savefig("figures/TestPhotoRing.png")

def testFisherDistribution():

    #AVERAGE AS A FUNCTION OF KAPPA
    tstd=[]
    ks=logspace(log10(0.1),log10(100),1000)
    for k in ks:
        thetas=randomFisher(kappa=k,nsample=10000)
        tstd+=[std(thetas)*RAD]
        
    fig=plt.figure()
    ax=fig.gca()
    ax.plot(ks,tstd,'b+')
    ax.set_xscale("log")
    ax.grid(which='both')
    ax.set_xlabel(r"$\kappa$",fontsize=20)
    ax.set_ylabel(r"$\sigma_{\theta}$ (degrees)",fontsize=20)
    ax.set_title("Fisher Distribution",position=(0.5,1.02))
    fig.savefig("figures/FisherStd.png")
    
    thetas=randomFisher(kappa=10.0,nsample=1000)
    mus=cos(thetas)
    fig=plt.figure()
    ax=fig.gca()
    ax.hist(mus,bins=20)
    fig.savefig("figures/Fisher.png")

def contourPhotoRing2():

    S=System
    try:
        qcalc=int(argv[1])
    except:
        qcalc=1

    print BARL,"Calculating Photo-ring effect contours",RBAR

    #////////////////////////////////////////
    #SYSTEM
    #////////////////////////////////////////
    RingedC=copyObject(Ringed)
    ap=Ringed.ap/Ringed.Rstar
    P=Ringed.Porb/HOUR
    ip=Ringed.iorb
    B=RingedC.Borb
    print "Impact parameter:",B
    print "Correction:",log10((1-B**2)**-0.75)
    print "Maximum:",log10(Ringed.fe**-1.5*(1-B**2)**-0.75)
    raw_input()
    Rp=RingedC.Rp
    rho_true=S.Mstar/(4*pi/3*S.Rstar**3)

    #////////////////////////////////////////
    #MAKE A MAP
    #////////////////////////////////////////
    cieffmin=0.01
    cieffmax=1.0
    Ncieffs=30
    teffmin=0.0*DEG
    teffmax=90.0*DEG
    Nteffs=30
    qcorrected=False

    if qcalc:
        cieffs=linspace(cieffmin,cieffmax,Ncieffs)
        teffs=linspace(teffmin,teffmax,Nteffs)
        #cieffs=[0.5];teffs=[45*DEG]
        
        IS,TS=meshgrid(cieffs,teffs)
        PR=zeros_like(IS)   
        PRN=zeros_like(IS)   
 
        ii=0
        for ci in cieffs:
            i=arccos(ci)
            print "Testing ieff (cos ieff) = ",i*RAD,ci
            jj=0
            for t in teffs:
                print TAB,"Testing teff = ",t*RAD
                RingedC.ieff=i

                RingedC.block=blockFactor(RingedC.tau,i)

                RingedC.Ringext.b=RingedC.Ringext.a*cos(i)
                RingedC.Ringext.cost=cos(t)
                RingedC.Ringext.sint=sin(t)

                RingedC.Ringint.b=RingedC.Ringint.a*cos(i)
                RingedC.Ringint.cost=cos(t)
                RingedC.Ringint.sint=sin(t)
                
                B=RingedC.Borb

                #NUMERICAL
                tcsp=contactTimes(RingedC)
                tT=(tcsp[-1]-tcsp[1])/HOUR
                tF=(tcsp[-2]-tcsp[2])/HOUR

                p=ringedPlanetArea(RingedC)/pi
                pRn=rhoObserved_Seager(p,
                                       RingedC.Rstar,
                                       tT,tF,
                                       P)/rho_true

                #ANALYTICAL
                xpa1=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                                     direction=-1,sign=-1,
                                     qcorrected=qcorrected)
                xpa2=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                                     direction=-1,sign=+1,
                                     qcorrected=qcorrected)
                xpa3=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                                     direction=+1,sign=-1,
                                     qcorrected=qcorrected)
                xpa4=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                                     direction=+1,sign=+1,
                                     qcorrected=qcorrected)
                tT=P*arcsin((xpa4-xpa1)/(ap*sin(ip)))/(2*pi)
                tF=P*arcsin((xpa3-xpa2)/(ap*sin(ip)))/(2*pi)

                p=analyticalTransitArea(RingedC.Rp,RingedC.block,
                                        RingedC.fi,RingedC.fe,
                                        i)/pi
                pR=rhoObserved_Seager(p,
                                      RingedC.Rstar,
                                      tT,tF,
                                      P)/rho_true

                print 2*TAB,"PR (Numerical) = %.6e, PR (Analytical) = %.6e"%(log10(pRn),log10(pR))

                PR[ii,jj]=log10(pR)
                PRN[ii,jj]=log10(pRn)

                jj+=1
            ii+=1

        savetxt("ISP.dat",IS)
        savetxt("TSP.dat",TS)
        savetxt("PR.dat",PR)
        savetxt("PRN.dat",PRN)

    IS=loadtxt("ISP.dat")
    TS=loadtxt("TSP.dat")
    PR=loadtxt("PR.dat")
    PRN=loadtxt("PRN.dat")
    data=loadtxt("posterior-PhotoRing-Concentrated.dat") #EXPERIMENT
    
    cmap=plt.get_cmap("rainbow")

    #////////////////////////////////////////
    #CONTOUR
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,6))
    ax=fig.gca()

    dmin=PRN.min()
    dmax=PRN.max()
    print "Minimum = ",dmin

    levels=linspace(dmin,dmax,100)
    c=ax.contourf(IS,TS*RAD,transpose(PRN),
                  levels=levels,cmap=cmap)

    """
    freq=50
    ax.plot(cos(data[::freq,4]*DEG),abs(data[::freq,7]),
            '+',color=cm.gray(0.2),zorder=5) #EXP
            """

    ax.set_xlabel(r"$\cos\,i$",fontsize=20)
    ax.set_ylabel(r"$\theta$",fontsize=20)
    ax.set_title("Photo-Ring Effect, b = %.2f"%(abs(B)),position=(0.5,1.02),
                 fontsize=14)

    cbar=fig.colorbar(c)
    cbar.ax.set_ylabel(r"$\log_{10}(\rho_{\rm obs}/\rho_\star)$",fontsize=20)
    yts=cbar.ax.get_yticks()
    yl=[]
    levels=[]
    for yt in yts:
        yv=(dmin+yt*(dmax-dmin))
        yl+=["%.3f"%((dmin+yt*(dmax-dmin))*1E0)]
        levels+=[yv]
    cbar.ax.set_yticklabels(yl)
    c=ax.contour(IS,TS*RAD,transpose(PRN),
               levels=levels,
               colors=['k'],linestyles=[':'])
    ax.clabel(c,inline=1,fontsize=10,zorder=50)
    levelsC=levels

    ax.contour(IS,TS*RAD,transpose(PRN),
               levels=[0.0],
               colors=['k'],linestyles=['-'],linewidths=['2'])
    
    plotPlanets(ax,RingedC,
                xmin=cieffmin,
                scalex=(cieffmax-cieffmin),
                ymin=teffmin*RAD,
                scaley=(teffmax-teffmin)*RAD)
    ax.set_xlim(cieffmin,cieffmax)
    ax.set_ylim(teffmin*RAD,teffmax*RAD)

    fig.savefig("figures/PhotoRingContour-Tests-%.2f.png"%B)

    #////////////////////////////////////////
    #CONTOUR
    #////////////////////////////////////////
    fig=plt.figure(figsize=(8,6))
    ax=fig.gca()

    dmin=PR.min()
    dmax=PR.max()
    
    levels=linspace(dmin,dmax,100)
    c=ax.contourf(IS,TS*RAD,transpose(PR),
                  levels=levels,cmap=cmap)

    ax.set_xlabel(r"$\cos\,i$",fontsize=20)
    ax.set_ylabel(r"$\theta$",fontsize=20)
    ax.set_title("Photo-Ring Effect, b = %.2f"%(abs(B)),position=(0.5,1.02),
                 fontsize=14)

    cbar=fig.colorbar(c)
    cbar.ax.set_ylabel(r"$\log(\rho_{\rm obs}/\rho_\star)$",fontsize=14)
    yts=cbar.ax.get_yticks()
    yl=[]
    levels=[]
    for yt in yts:
        yv=(dmin+yt*(dmax-dmin))
        yl+=["%.3f"%((dmin+yt*(dmax-dmin))*1E0)]
        levels+=[yv]
    cbar.ax.set_yticklabels(yl)
    c=ax.contour(IS,TS*RAD,transpose(PR),
                 levels=levelsC,
                 colors=['k'],linestyles=[':'])
    ax.clabel(c,inline=1,fontsize=10)

    ax.contour(IS,TS*RAD,transpose(PR),
               levels=[0.0],
               colors=['k'],linestyles=['-'],linewidths=['2'])
    
    plotPlanets(ax,RingedC,
                xmin=cieffmin,
                scalex=(cieffmax-cieffmin),
                ymin=teffmin*RAD,
                scaley=(teffmax-teffmin)*RAD)
    ax.set_xlim(cieffmin,cieffmax)
    ax.set_ylim(teffmin*RAD,teffmax*RAD)

    fig.savefig("figures/PhotoRingContour-Analytical-Tests-%.2f.png"%B)

#testTransitDepth()
#testTransitDuration()
#errorTransitPositions()
#errorTransitTimes()
#contourTransitDepths()
#curveTransitDepths()
#transitDepthPosterior()
#testPhotoRing()
#contourPhotoRing()
photoRingPosterior()
#photoRingPosteriorTiming()
#testFisherDistribution()
#contourPhotoRing2()
