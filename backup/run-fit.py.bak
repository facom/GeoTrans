from geotrans import *
from system import *
print BARL,"Fit light curve",RBAR
S=System
Planet=onlyPlanet(System)
FVERBOSE=1

#########################################
#FITTING ROUTINE
#########################################
def residualTransitPlanet(params,
                          times,signal,dsignal,
                          S):
    Rp=params['Rp'].value
    fp=params['fp'].value
    stp=params['stp'].value
    ctp=(1-stp**2)**0.5
    S.Planet.a=Rp
    S.Planet.b=Rp
    if EQUAL(fp,ZERO):
        transitFunc=areaStriping
        qoblate=False
    else:
        if FVERBOSE:print "Testing: Rp = %e, fp = %e, stp = %e"%(Rp,fp,stp)
        S.Planet.b*=(1-fp)
        S.Planet.sint=stp
        S.Planet.cost=ctp
        transitFunc=areaOblateStriping
        qoblate=True

    Ar=FIGUREAREA(S.Planet)
    tcs=contactTimes(S)
    transit=[]
    for t in times:
        a=fluxLimbTime(t+S.tcen,Ar,S,areas=transitFunc)
        transit+=[a]
    transit=array(transit)
    residual=transit-signal
    if qoblate:
        chisq=(residual**2).sum()
        if FVERBOSE:print "Chi-square = %e"%chisq
    return residual

#########################################
#OPTIONS
#########################################
try:qfit_obl=int(argv[1])
except:qfit_obl=1
try:qfit_sph=int(argv[2])
except:qfit_sph=1

#########################################
#LOAD SIGNAL DATA
#########################################
data=loadtxt("lightcurve.dat")
times=data[:,0]
signal=data[:,1]
dsignal=data[:,2]
data=loadtxt("transit.dat")
ttimes=data[:,0]
tsignal=data[:,1]
tzeroes=zeros_like(tsignal)

#########################################
#RING RESIDUAL
#########################################
transitRinged=interpolant(ttimes,tsignal)
residual_rng=signal-transitRinged(times)

#########################################
#PARAMETERS
#########################################
params=Parameters()
#EQUATORIAL RADIUS
params.add("Rp",value=0.17813250,
           min=S.Planet.a/2,
           max=1.0)
#OBLATENESS
params.add("fp",value=0.63893422,
           min=0.0,
           max=1.0)
#SINE OF PROJECTED ROLL
params.add("stp",value=0.0,
           min=-1.0,
           max=+1.0,
           vary=False)

#########################################
#FIT OBLATE
#########################################
print "Fitting oblate..."
#params['Rp'].value=0.1715
params['Rp'].value=0.11
if qfit_obl:
    result=minimize(residualTransitPlanet,
                    params,
                    args=(times,signal,dsignal,Planet))
report_fit(params)
Rp_obl=params['Rp'].value
fp_obl=params['fp'].value
stp_obl=params['stp'].value
residual_obl=residualTransitPlanet(params,times,
                                   signal,dsignal,
                                   Planet)
bestfit_obl=residualTransitPlanet(params,ttimes,
                                  tzeroes,tzeroes,
                                  Planet)
#########################################
#FIT SPHERICAL
#########################################
print "Fitting spherical..."
params['Rp'].value=0.103033
params['fp'].value=0.0
params['fp'].vary=False
params['stp'].value=0.0
params['stp'].vary=False
if qfit_sph:
    result=minimize(residualTransitPlanet,
                    params,
                    args=(times,signal,dsignal,Planet))
report_fit(params)
Rp_sph=params['Rp'].value
residual_sph=residualTransitPlanet(params,times,
                                  signal,dsignal,
                                  Planet)
bestfit_sph=residualTransitPlanet(params,ttimes,
                                  tzeroes,tzeroes,
                                  Planet)

#########################################
#PLOT SIGNAL
#########################################
fig=plt.figure(figsize=(8,8))
ax=fig.gca()

#DATA POINTS
ax.plot(times/HOUR,signal,'ko',markersize=2)
ax.plot(ttimes/HOUR,tsignal,'r-',label='Rings, $R_p/R_\star$ = %.3f'%S.Rp)
ax.plot(ttimes/HOUR,bestfit_sph,'b-',label='Spherical, $R_p/R_\star$ = %.3f'%Rp_sph)
ax.plot(ttimes/HOUR,bestfit_obl,'c-',label='Oblate, $R_p/R_\star$ = %.3f, $f_p$ = %.3f'%(Rp_obl,fp_obl))

#//////////////////////////////
#DECORATION SIGNAL
#//////////////////////////////
ymin=min(signal);ymax=max(signal)
tmin=min(times/HOUR);tmax=max(times/HOUR)
ax.set_xlim((tmin,tmax))
ax.set_ylim((ymin,ymax))
ax.set_xlabel(r"$t-t_{\rm cen}$ (hours)",fontsize=14)
ax.set_xlabel(r"Relative Flux",fontsize=14)
ax.legend(loc='best',prop=dict(size=10))
ax.grid()
fig.savefig("plots/fit_bestfit.png")

#########################################
#PLOT RESIDUALS
#########################################
fig=plt.figure(figsize=(8,8))

l=0.1;b=0.08;w=0.85;h=0.28;s=0.02
ax_rng=fig.add_axes([l,b,w,h]);b+=h+s
ax_sph=fig.add_axes([l,b,w,h]);b+=h+s
ax_obl=fig.add_axes([l,b,w,h]);b+=h+s

#DATA POINTS
ax_rng.plot(times/HOUR,residual_rng,'ro',markersize=1,label='Residual Spherical')
ax_rng.errorbar(times/HOUR,residual_rng,yerr=dsignal,linestyle='none',color='r',alpha=0.3)

ax_sph.plot(times/HOUR,residual_sph,'bo',markersize=1,label='Residual Spherical')
ax_sph.errorbar(times/HOUR,residual_sph,yerr=dsignal,linestyle='none',color='b',alpha=0.3)

ax_obl.plot(times/HOUR,residual_obl,'co',markersize=1,label='Residual Oblate')
ax_obl.errorbar(times/HOUR,residual_obl,yerr=dsignal,linestyle='none',color='c',alpha=0.3)

#FLUX RANGE
ymin=1.0;ymax=-1.0
for res in residual_rng,residual_sph,residual_obl:
    ymin=min(ymin,min(res))
    ymax=max(ymax,max(res))

for ax in ax_rng,ax_sph,ax_obl:
    #LABELS
    yt=ax.get_yticks()
    yl=[]
    for y in yt:
        yl+=["%.0f"%(y*1E6)]
    ax.set_yticklabels(yl)
    ax.set_ylim((ymin,ymax))
    ax.set_xlim((tmin,tmax))
    ax.axhline(0,color='k')
    ax.set_ylabel('Residuals (ppm)')
    ax.grid()

for ax in ax_sph,ax_obl:ax.set_xticks([])

ax_rng.set_xlabel(r"$t-t_{\rm cen}$ (hours)",fontsize=14)
fig.savefig("plots/fit_residual.png")
