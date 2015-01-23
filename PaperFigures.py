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
          -90.0*DEG,
          +90.0*DEG,
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

    print "Fixed values:"
    print TAB,"i = %.2f deg"%(Ringed.ieff*RAD)
    print TAB,"t = %.2f deg"%(Ringed.teff*RAD)

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

    t=Ringed.teff
    i=Ringed.ieff

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
    #FIX RINGEDC PROPERTIES BY HAND
    #========================================
    #MANUAL i,t
    i=89.427032*DEG;t=0*DEG
    RingedC.Ringext.b=RingedC.Ringext.a*cos(i)
    RingedC.Ringext.cost=cos(t);RingedC.Ringext.sint=sin(t)

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


    ax.set_xlabel(r"$\cos\,i$",fontsize=16)
    ax.set_ylabel(r"$R_{\rm p,obs}/R_{\rm p}$",fontsize=16)

    
    psize=0.02
    plotPlanets(ax,RingedC,Nx=6,Ny=2,ymin=0,yoffset=0.1,
                fh=psize/RingedC.Rp,fv=psize/RingedC.Rp)
    ax.set_xlim(cieffmin,cieffmax)
    ax.set_ylim(1.0,sqrt(RingedC.Ringext.a**2-RingedC.Ringint.a**2+Rp**2)/Rp)

    ax.legend(loc="best")
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

    Nplanets=2
    Nsamples=5

    Nbins=30

    Ntotal=Nplanets*Nsamples

    PARAMETERS["iorb"][DEF]=0.0
    #PARAMETERS["ir"][DEF]=0.98
    #PARAMETERS["phir"][DEF]=-40.0*DEG

    PARAMETERS["ir"][STAT]=VAR
    PARAMETERS["phir"][STAT]=VAR
    
    #########################################
    #SAMPLE GENERATION
    #########################################
    if qcalc:

        i=1
        header=""
        header+="%-17s\t"%("#0:id")
        for parkey in PARKEYS:
            parameter=PARAMETERS[parkey]
            scale=parameter[SCAL]
            if parameter[STAT]:
                header+="%-17s\t"%("%d:%s(%s)"%(i,parkey,scale))
                i+=1

        for propkey in PROPKEYS:
            properti=PROPERTIES[propkey]
            scale=properti[SCAL]
            if properti[STAT]:
                header+="%-18s\t"%("%d:%s(%s)"%(i,propkey,scale))
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

            """
            #Z-AXIS: 
            #print dot(S.Mrs,[0,0,1])

            fig=plt.figure(figsize=(8,8))
            ax=fig.gca()
            
            plotEllipse(ax,S.Star,color='y')
            plotEllipse(ax,S.Planet,color='b')
            plotEllipse(ax,S.Ringext,color='k')
            plotEllipse(ax,S.Ringint,color='r')
            
            rng=1.5
            Re=1.0
            xmin=S.Planet.C[0]-rng*Re;
            xmax=S.Planet.C[0]+rng*Re
            ymin=S.Planet.C[1]-rng*Re;
            ymax=S.Planet.C[1]+rng*Re
            ax.set_xlim((xmin,xmax))
            ax.set_ylim((ymin,ymax))
            ax.grid()
            fig.savefig("tmp/systems/system-%d.png"%i)
            """
            
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
                      header,line,fmt="%+.17e")
    
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

    xms=histPlot(ax,xs,hs,dhs,error=error,color='r',alpha=0.1)
    hms=softArraySG(hs,frac=2,nP=3)
    ax.plot(xms,hms,'b-',linewidth=2,zorder=10)
    ax.set_xlabel(r"$R_{\rm p,obs}/R_{\rm p}$",fontsize=14)
    ax.set_ylabel("Frequency",fontsize=12)
    ax.set_title(r"Observed Radius Posterior Distribution",
                 position=(0.5,1.02))
    ax.axvline(0.0,linewidth=2)
    ax.set_xlim((xs[0],xs[-1]))
    hmin,hmax=ax.get_ylim()
    ax.set_ylim((0.0,hmax))
    fig.savefig("figures/posterior-TransitDepth.png")

def contourPhotoRing():

    S=System
    try:
        qcalc=int(argv[1])
    except:
        qcalc=1

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
    dp=S.Rp**2
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

    if qcalc:
        #NOT-RINGED
        Ap=pi*Rp**2

        cieffs=linspace(cieffmin,cieffmax,Ncieffs)
        teffs=linspace(teffmin,teffmax,Nteffs)

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
                pRn=rhoObserved_Seager(RingedC.Rp,
                                       RingedC.Rstar,
                                       tT,tF,
                                       P)/rho_true

                #ANALYTICAL
                xpa1=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                                     direction=-1,sign=-1,qcorrected=True)
                xpa2=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                                     direction=-1,sign=+1,qcorrected=True)
                xpa3=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                                     direction=+1,sign=-1,qcorrected=True)
                xpa4=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                                     direction=+1,sign=+1,qcorrected=True)
                tT=P*arcsin((xpa4-xpa1)/(ap*sin(ip)))/(2*pi)
                tF=P*arcsin((xpa3-xpa2)/(ap*sin(ip)))/(2*pi)
                pR=rhoObserved_Seager(RingedC.Rp,
                                      RingedC.Rstar,
                                      tT,tF,
                                      P)/rho_true

                print 2*TAB,"PR (Numerical) = %.6e, PR (Analytical) = %.6e"%(pRn,pR)

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

    ax.set_xlabel(r"$\cos\,i$",fontsize=20)
    ax.set_ylabel(r"$\theta$",fontsize=20)
    ax.set_title("Photo-Ring Effect (Numerical calculation)"%(B),position=(0.5,1.02),
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
    c=ax.contour(IS,TS*RAD,transpose(PRN),
               levels=levels,
               colors=['k'],linestyles=[':'])
    ax.clabel(c,inline=1,fontsize=10)
    
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
    ax.contour(IS,TS*RAD,transpose(PR),
               levels=levels,
               colors=['k'],linestyles=[':'])

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

    Nplanets=2
    Nsamples=5

    Nbins=30

    Ntotal=Nplanets*Nsamples

    PARAMETERS["iorb"][DEF]=0.0
    #PARAMETERS["ir"][DEF]=0.98
    #PARAMETERS["phir"][DEF]=-40.0*DEG

    PARAMETERS["ir"][STAT]=VAR
    PARAMETERS["phir"][STAT]=VAR
    
    PROPERTIES['r'][STAT]=HIDE

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
            S.PR=rhoObserved_Seager(S.Rp,S.Rstar,
                                    tT,tF,S.Porb/HOUR)/rho_true
            S.logPR=log10(S.PR)

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

    xms=histPlot(ax,xs,hs,dhs,error=error,color='r',alpha=0.1)
    hms=softArraySG(hs,frac=2,nP=3)
    ax.plot(xms,hms,'b-',linewidth=2,zorder=10)
    ax.set_xlabel(r"$R_{\rm p,obs}/R_{\rm p}$",fontsize=14)
    ax.set_ylabel("Frequency",fontsize=12)
    ax.set_title(r"Observed Radius Posterior Distribution",
                 position=(0.5,1.02))
    ax.axvline(0.0,linewidth=2)
    ax.set_xlim((xs[0],xs[-1]))
    hmin,hmax=ax.get_ylim()
    ax.set_ylim((0.0,hmax))
    fig.savefig("figures/posterior-PhotoRing.png")

def testPhotoRing():

    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()

    S=Ringed
    t=Ringed.teff
    i=Ringed.ieff

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
    rho_true=S.Mstar/(4*pi/3*S.Rstar**3)
    dp=S.Rp**2
    print "True density = ",rho_true
    
    #========================================
    #NOT RINGED DENSITY
    #========================================
    xp=sqrt((1+NotRinged.Rp)**2-NotRinged.Borb**2)
    xm=sqrt((1-NotRinged.Rp)**2-NotRinged.Borb**2)
    tT=P*arcsin(2*xp/(ap*sin(ip)))/(2*pi)
    tF=P*arcsin(2*xm/(ap*sin(ip)))/(2*pi)
    
    aR_obs=2*S.Porb/pi*\
        dp**0.25/((tT*HOUR)**2-(tF*HOUR)**2)**0.5
    b_obs=(((1-sqrt(dp))**2-(tF/tT)**2*(1+sqrt(dp))**2)/\
               (1-(tF/tT)**2))**0.5
    rho_obs=3*pi*aR_obs**3/(GCONST*S.Porb**2) #Kipping, 2013

    rho_obs=rhoObserved_Seager(S.Rp,S.Rstar,
                               tT,tF,S.Porb/HOUR)
    print "Observed density (not ringed, Seager) = ",rho_obs
    rho_obs=rhoObserved_Kipping(S.Rp,S.Rstar,
                               tT,tF,S.Porb/HOUR)
    print "Observed density (not ringed, Kipping) = ",rho_obs

    #========================================
    #FIX RINGEDC PROPERTIES BY HAND
    #========================================
    #MANUAL i,t
    
    #GOOD SPOT
    i=53*DEG;t=80*DEG

    #GOOD SPOT
    i=10*DEG;t=80*DEG

    i=89.427033*DEG;t=0*DEG
    
    RingedC.Ringext.b=RingedC.Ringext.a*cos(i)
    RingedC.Ringint.b=RingedC.Ringint.a*cos(i)
    RingedC.Ringext.cost=cos(t);RingedC.Ringext.sint=sin(t)
    RingedC.Ringint.cost=cos(t);RingedC.Ringint.sint=sin(t)

    #========================================
    #RINGED TRANSIT DURATION (NUMERICAL)
    #========================================
    lw=1

    tcsp=contactTimes(RingedC)
    tT=(tcsp[-1]-tcsp[1])/HOUR
    tF=(tcsp[-2]-tcsp[2])/HOUR

    rho_obs=rhoObserved_Kipping(S.Rp,S.Rstar,
                               tT,tF,S.Porb/HOUR)
    pR=rho_obs/rho_true
    print "Observed density (ringed, numerical, Kipping) = ",rho_obs
    print "Photo-ring effect (numerical) = ",pR

    #========================================
    #RINGED TRANSIT ANALYTICAL 
    #========================================
    a=RingedC.Ringext.a
    b=RingedC.Ringext.b
    B=RingedC.Borb

    #C1
    xpa1=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                         direction=-1,sign=-1,qcorrected=True)
    #C2
    xpa2=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                         direction=-1,sign=+1,qcorrected=True)
    #C3
    xpa3=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                         direction=+1,sign=-1,qcorrected=True)
    #C4
    xpa4=transitPosition(RingedC.Rp,RingedC.fe,i,t,B,
                         direction=+1,sign=+1,qcorrected=True)
    taT=P*arcsin((xpa4-xpa1)/(ap*sin(ip)))/(2*pi)
    taF=P*arcsin((xpa3-xpa2)/(ap*sin(ip)))/(2*pi)

    rho_obsa=rhoObserved_Kipping(S.Rp,S.Rstar,
                                 taT,taF,S.Porb/HOUR)
    pRa=rho_obsa/rho_true
    print "Observed density (ringed, analytical, Kipping) = ",rho_obsa
    print "Photo-ring effect (analytical) = ",pRa

    #========================================
    #ERRORS
    #========================================
    drho=abs(rho_obs-rho_obsa)/rho_obs*100
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
