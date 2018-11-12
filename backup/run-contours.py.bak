from geotrans import *
from system import *
print(BARL,"Creating Transit Parameters Contours",RBAR)

#########################################
#NOT RINGED PLANET
#########################################
S=System
NotRinged=copyObject(S)
NotRinged.Ringext.b=NotRinged.Ringint.b=0.0
print("ieff,teff = ",S.ieff*RAD,S.teff*RAD)


#########################################
#CONTACT TIMESCALES
#########################################
print("Single planet parameters:")
print(TAB,"Calculating contact times...")
tcsp=contactTimes(NotRinged)
t14p=(tcsp[-1]-tcsp[1])/HOUR
t23p=(tcsp[-2]-tcsp[2])/HOUR
dp=S.Rp**2
print(TAB,"T_14 = %e h = %e s"%(t14p,t14p*HOUR))
print(TAB,"Transit depth = %e"%(dp))

#THEORETICAL
aR_obs=2*S.Porb/pi*\
    dp**0.25/((t14p*HOUR)**2-(t23p*HOUR)**2)**0.5
b_obs=(((1-sqrt(dp))**2-(t23p/t14p)**2*(1+sqrt(dp))**2)/\
    (1-(t23p/t14p)**2))**0.5
rho_obs=3*pi*aR_obs**3/(GCONST*S.Porb**2) #Kipping, 2013

#TRUE
rho_true=S.Mstar/(4*pi/3*S.Rstar**3)
aR_true=S.ap/S.Rstar
b_true=S.Borb

print(TAB,"rho_star(obs) = %e"%rho_obs)
print(TAB,"rho_star(true) = %e"%rho_true)

print(TAB,"a/Rstar (obs) = %e"%aR_obs)
print(TAB,"a/Rstar (true) = %e"%aR_true)

print(TAB,"b (obs) = %e"%b_obs)
print(TAB,"b (true) = %e"%b_true)

print(TAB,"Planet transit duration: T14 = %.3fh, T23 = %.3fh"%(t14p,t23p))

tcsr=contactTimes(S)
t14ref=(tcsr[-1]-tcsr[1])/HOUR
t23ref=(tcsr[-2]-tcsr[2])/HOUR

print(TAB,"Contact timescales: T14 = %.3fh, T23 = %.3fh"%(t14ref,t23ref))

############################################
#PLOT TRANSIT DEPTH AS A FUNCTION OF PR.INCL
############################################
print()
print("Calculating effect on transit depth...")
ieffs=linspace(0.0,90.0,100)*DEG
drs=[]
for ieff in ieffs:
    S.Ringext.b=S.Re*cos(ieff)
    S.Ringint.b=S.Ri*cos(ieff)
    S.Ringext.ct=S.Ringint.ct=1.0
    S.Ringext.st=S.Ringint.st=0.0
    Ar=ringedPlanetArea(S)
    dpr=Ar/pi
    drs+=[sqrt(dpr/dp)]

fig=plt.figure(figsize=(8,6))
ax=fig.gca()
ax.plot(ieffs*RAD,drs,'b-')
ax.set_xlabel(r"Projected inclination, $i_{\rm eff}$ (degrees)")
ax.set_ylabel(r"$\sqrt{d_{\rm Ring}/d_{\rm Planet}}$")
ax.set_xlim((0,90))
ax.set_ylim((min(drs),max(drs)))
fig.savefig("plots/depth_projieff.png")

#########################################
#CONTOURS
#########################################
try:qsave=int(argv[1])
except:qsave=1

print()
print("Calculating contours...")
print()
if qsave:
    #ieffs=linspace(0.0,90.0,30)*DEG
    cieffs=linspace(0.0,1.0,30)
    teffs=linspace(0.0,90.0,30)*DEG
    #ieffs=linspace(30.0,90.0,10)*DEG
    #teffs=linspace(0.0,90.0,10)*DEG
    #ieffs=array([30,60])*DEG
    #teffs=array([30,60])*DEG
    #ieffs=array([30])*DEG;teffs=array([30])*DEG
    #ieffs=array([90])*DEG;teffs=array([45])*DEG
    cieffs=array([0.5]);teffs=array([45])*DEG
    
    IS,TS=meshgrid(cieffs,teffs*RAD)
    DST14=zeros_like(IS)
    DST23=zeros_like(IS)
    DSRHO=zeros_like(IS)
    
    i=0
    VERBCONT=1
    for cieff in cieffs:
        ieff=arccos(cieff)
        print("Projected Inclination: ieff = %.2f"%(ieff*RAD))
        j=0

        S.Ringext.b=S.Re*cieff
        S.Ringint.b=S.Ri*cieff
        S.Ringext.cost=S.Ringint.cost=1.0
        S.Ringext.sint=S.Ringint.sint=0.0
        Ar=ringedPlanetArea(S)
        dpr=Ar/pi

        for teff in teffs:
            if VERBCONT:print(TAB,"Projected Roll: t = %.2f (ieff=%.2f)"%(teff*RAD,ieff*RAD))
            S.Ringext.cost=S.Ringint.cost=cos(teff)
            S.Ringext.sint=S.Ringint.sint=sin(teff)
            
            ts=contactTimes(S)
            t14=(ts[-1]-ts[1])/HOUR
            t23=(ts[-2]-ts[2])/HOUR
            if VERBCONT:print(2*TAB,"Transit times = ",t14,t23)
            
            aR_obsR=2*S.Porb/pi*\
                dpr**0.25/((t14*HOUR)**2-(t23*HOUR)**2)**0.5
            rho_obsR=3*pi*aR_obsR**3/(GCONST*S.Porb**2)
            if VERBCONT:print(2*TAB,"a/Rstar (obs,ring) = ",aR_obsR)
            if VERBCONT:print(2*TAB,"rho(obs,ring) = ",rho_obsR)

            DST14[i,j]=t14/t14p-1
            DST23[i,j]=t23/t23p-1
            if VERBCONT:print(2*TAB,"Comparison: dt14 = %e. dt23 = %e"%(DST14[i,j],DST23[i,j]))

            DSRHO[i,j]=log10(rho_obsR/rho_obs)
            if VERBCONT:print(2*TAB,"Astrodensity: rhoR/rhoP = %e"%(DSRHO[i,j]))            
            if VERBCONT:exit(0)

            j+=1
        i+=1

    savetxt("contours-IS.dat",IS)
    savetxt("contours-TS.dat",TS)
    savetxt("contours-DST14.dat",DST14)
    savetxt("contours-DST23.dat",DST23)
    savetxt("contours-DSRHO.dat",DSRHO)
    
IS=loadtxt("contours-IS.dat")
TS=loadtxt("contours-TS.dat")
DST14=transpose(loadtxt("contours-DST14.dat"))
DST23=transpose(loadtxt("contours-DST23.dat"))
DSRHO=transpose(loadtxt("contours-DSRHO.dat"))

#########################################
#PLOT
#########################################
Nlevels=100
cmap=plt.get_cmap("Spectral")

def plotPlanets(ax):
    #ieffs=array([0,20,45,65,90])*DEG
    cieffs=array([0,0.25,0.50,0.75,1.0])
    teffs=array([0,20,45,65,90])*DEG
    for cieff in cieffs:
        ieff=arccos(cieff)
        for teff in teffs:
            fh=0.03/S.Rp
            fv=fh
            C=AR(cieff,teff*RAD/90)
            Planet=Figure(C,fh*S.Rp,fv*S.Rp,1.0,0.0,'Planet')
            Ringe=Figure(C,fh*S.Re,fv*S.Re*cos(ieff),cos(teff),sin(teff),'Ringext')
            Ringi=Figure(C,fh*S.Ri,fv*S.Ri*cos(ieff),cos(teff),sin(teff),'Ringint')
            plotEllipse(ax,Planet,patch=True,zorder=10,color='b',transform=ax.transAxes)
            plotEllipse(ax,Ringe,zorder=10,color='k',transform=ax.transAxes)
            plotEllipse(ax,Ringi,zorder=10,color='k',transform=ax.transAxes)

#//////////////////////////////
#TRANSIT DURATION
#//////////////////////////////
fig=plt.figure(figsize=(16,12))
ax=fig.gca()

dmin=DST14.min()*t14p
dmax=DST14.max()*t14p
levels=linspace(dmin,dmax,Nlevels)

c=ax.contourf(IS,TS,DST14*t14p,levels=levels,cmap=cmap)
fig.colorbar(c)

levels=[15*MINUTE/(t14ref*HOUR)]
levels=[15*MINUTE/HOUR]
c=ax.contour(IS,TS,DST14*t14p,levels=levels,color='k',zorder=10,linewidth=4)

ax.set_title(r"$t_{14,\rm{Ring}}-t_{14,\rm{Planet}}$ (hours)",position=(0.5,1.01),fontsize=14)
ax.set_xlabel(r"Projected Inclination, $i_{\rm eff}$ (degrees)",fontsize=14)
ax.set_ylabel(r"Projected Roll, $\theta_{\rm R}$ (degrees)",fontsize=14)

plotPlanets(ax)
ax.grid()
ax.set_xlim((0,1))
ax.set_ylim((0,90))
fig.savefig("plots/ring-contour-ts14.png")

#//////////////////////////////
#CENTRAL PART
#//////////////////////////////
fig=plt.figure(figsize=(16,12))
ax=fig.gca()

dmin=DST23.min()*t23p
dmax=DST23.max()*t23p
levels=linspace(dmin,dmax,Nlevels)

c=ax.contourf(IS,TS,DST23*t23p,levels=levels,cmap=cmap)
fig.colorbar(c)

levels=[-15*MINUTE/(t23ref*HOUR)]
levels=[-15*MINUTE/HOUR]
c=ax.contour(IS,TS,DST23*t23p,levels=levels,color='k',zorder=10,linewidth=4)

ax.set_title(r"$t_{23,\rm{Ring}}-t_{23,\rm{Planet}}$ (hours)",position=(0.5,1.01),fontsize=14)
ax.set_xlabel(r"Projected Inclination, $i_{\rm eff}$ (degrees)",fontsize=14)
ax.set_ylabel(r"Projected Roll, $\theta_{\rm R}$ (degrees)",fontsize=14)

plotPlanets(ax)
ax.grid()
ax.set_xlim((0,1))
ax.set_ylim((0,90))
fig.savefig("plots/ring-contour-ts23.png")

#//////////////////////////////
#ASTRODENSITY
#//////////////////////////////
fig=plt.figure(figsize=(16,12))
ax=fig.gca()

dmin=DSRHO.min()
dmax=DSRHO.max()
levels=linspace(dmin,dmax,Nlevels)

c=ax.contourf(IS,TS,DSRHO,levels=levels,cmap=cmap)
cbar=fig.colorbar(c)
cbar.ax.set_ylabel(r"$log(\rho_{\rm obs}/\rho_{\rm true})$",fontsize=20)
xt=cbar.ax.get_xticks()
xl=[]
for x in xt:
    xl+=["%+.3f"%x]
cbar.ax.set_xticklabels(xl)

levels=[0]
c=ax.contour(IS,TS,DSRHO,levels=levels,cmap=cmap,
             color='k',zorder=10,linewidth=4)

ax.set_xlabel(r"Projected Inclination, $cos(i_{\rm eff})$ (degrees)",fontsize=20)
ax.set_ylabel(r"Projected tilt $\theta_{\rm R}$ (degrees)",fontsize=20)

plotPlanets(ax)
ax.grid()
ax.set_xlim((0,1))
ax.set_ylim((0,90))
fig.savefig("plots/ring-contour-rhos.png")
