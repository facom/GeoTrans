from geotrans import *
from system import *
print BARL,"Distribution of PR",RBAR

#########################################
#SYSTEM
#########################################
S=System
NotRinged=copyObject(S)
NotRinged.Ringext.b=NotRinged.Ringint.b=0.0

#########################################
#PARAMETERS
#########################################
Nplanets=500
Nsamples=5
#Nbins=Nplanets/10
Nbins=30
Ntotal=Nplanets*Nsamples
verbose=True
verbose=False
qcalc=True
qcalc=False #Uncomment if just plot

fac=1/1.5
fac=1.0
fac=1.5
S.fe=2.5*fac
S.fi=1.5*fac

#########################################
#INPUT VARIABLES
#########################################
imax=arctan(S.ap/(S.Rstar-2*S.fe*S.Rplanet))
cimax=cos(imax)

#########################################
#COMPARISON VALUES
#########################################
tcsp=contactTimes(NotRinged)
t14p=(tcsp[-1]-tcsp[1])/HOUR
t23p=(tcsp[-2]-tcsp[2])/HOUR
dp=S.Rp**2
aR_obs=2*S.Porb/pi*\
    dp**0.25/((t14p*HOUR)**2-(t23p*HOUR)**2)**0.5
b_obs=(((1-sqrt(dp))**2-(t23p/t14p)**2*(1+sqrt(dp))**2)/\
    (1-(t23p/t14p)**2))**0.5
rho_true=3*pi*aR_obs**3/(GCONST*S.Porb**2) #Kipping, 2013
print "True density = ",rho_true

#########################################
#SAMPLE GENERATION
#########################################
if qcalc:
    PRs=[]
    cirs=[]
    phirs=[]
    ciorbs=[]
    cieffs=[]
    teffs=[]
    Borbs=[]
    for i in xrange(Nplanets*Nsamples):
        if (i%Nplanets)==0:
            print "Sample %d: %d planets generated..."%(i/Nplanets,i)
        #RANDOM INPUT PARAMETERS
        cir=randomVal(0.0,1.0)
        ir=arccos(cir)
        ciorb=randomVal(-cimax,cimax)
        iorb=arccos(ciorb)
        phir=randomVal(0.0,360.0)*DEG
        
        if verbose:
            print "Testing values:"
            print TAB,"ir = ",ir*RAD
            print TAB,"phir = ",phir*RAD
            print TAB,"iorb = ",iorb*RAD
        
        #CONVERTING INTO OBSERVED
        S.ir=ir
        S.phir=phir
        S.iorb=iorb
        derivedSystemProperties(S)
        updatePlanetRings(S)
        
        if verbose:
            print "Derived Properties:"
            print TAB,"ieff = ",S.ieff*RAD
            print TAB,"teff = ",S.teff*RAD
            print TAB,"Borb = ",S.Borb
            
        #STORING
        cirs+=[cir]
        ciorbs+=[ciorb]
        phirs+=[phir]
        teffs+=[S.teff]
        cieffs+=[cos(S.ieff)]
        Borbs+=[S.Borb]
        
        #PHOTO-RING EFFECT
        Ar=ringedPlanetArea(S)
        dpr=Ar/pi
        ts=contactTimes(S)
        t14=(ts[-1]-ts[1])/HOUR
        t23=(ts[-2]-ts[2])/HOUR
        aR_obsR=2*S.Porb/pi*\
            dpr**0.25/((t14*HOUR)**2-(t23*HOUR)**2)**0.5
        rho_obs=3*pi*aR_obsR**3/(GCONST*S.Porb**2)
        PR=log10(rho_obs/rho_true)
        
        PRs+=[PR]
    
    cirs=array(cirs)
    phirs=array(phirs)*RAD
    ciorbs=array(ciorbs)
    
    teffs=array(teffs)*RAD
    cieffs=array(cieffs)
    Borbs=array(Borbs)
    
    PRs=array(PRs)

    data=transpose(vstack((cirs,phirs,ciorbs,
                           cieffs,teffs,Borbs,
                           PRs)))
    savetxt("PR-sample-fe_%.1f.dat"%S.fe,data)

data=loadtxt("PR-sample-fe_%.1f.dat"%S.fe)
cirs=data[:,0]
phirs=data[:,1]
ciorbs=data[:,2]
cieffs=data[:,3]
teffs=data[:,4]
Borbs=data[:,5]
PRs=data[:,6]

#########################################
#HISTOGRAM AND ERRORS
#########################################
xs,hs,dhs=histPosterior(PRs,Nsamples,nbins=Nbins,
                        normed=True)

#########################################
#HISTOGRAM PR-EFFECT
#########################################
fig=plt.figure()
ax=fig.gca()
error=True
histPlot(ax,xs,hs,dhs,error=error,color='r')
ax.set_xlabel(r"$\log(\rho^{\rm PR}_{\rm obs}/\rho_{\rm true})$")
ax.set_ylabel("Frequency")
ax.set_title(r"$R_{\rm Ring}/R_{\rm P}=%.2f$"%(S.fe),position=(0.5,1.02))
ax.axvline(0.0,linewidth=2)
ax.set_xlim((xs[0],xs[-1]))
fig.savefig("plots/PR-histogram--fe_%.1f.png"%S.fe)
