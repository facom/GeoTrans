#!/usr/bin/env python
# coding: utf-8

# # GeoTrans
# ## Syntethic Exoring Light-curve
# 
# Using this script you will be able to generate a synthetic light-curve of a planet with an exoring.

# ### Import modules and load magics

# In[12]:


from geotrans import *
get_ipython().magic('matplotlib notebook')
get_ipython().magic('reload_ext autoreload')
get_ipython().magic('autoreload 2')
et=elapsed()


# ### Configuration

# In[2]:


if QIPY:
    CONF.SYSTEM="kepler421"


# ### Import system information

# In[3]:


exec(f"from data.{CONF.SYSTEM} import *")
S=System
Snr=onlyPlanet(S)
systemShow(S)


# ### Observational parameters

# In[4]:


#Cadence
tcad=15.0*MINUTE 
#Number of observations during transit
nobs=int(ceil((S.tmax-S.tmin)/tcad)) 
#Number of counts per cadence
Nexp=tcad*S.Flux
print(f"Number of observations per transit: {nobs}")
print(f"Average number of counts per observation: {Nexp}")


# ### Worker function

# In[13]:


def computeFlux(i,ts,Ts,s,ds):
    Sw=deepcopy(S)
    cs=random.poisson([gaussianQuadrature(transitFlux,t,t+tcad,args=(Sw,)) for t in tqdm(ts,position=i)])/Nexp
    cds=(cs/Nexp)**0.5
    s+=cs.tolist()
    ds+=cds.tolist()
    Ts+=ts.tolist()


# ### Parallel execution

# In[6]:


tini=0
#tend=2*YEAR
tend=4*S.Porb
times=arange(tini,tend,tcad)
print(f"Number of times: {len(times)}")


# In[ ]:


if QIPY:
    CONF.NP=1
    times=times[:100]


# In[7]:


ts=[]
s=[]
ds=[]
et=elapsed()
runParallel(computeFlux,times,ts,s,ds,CONF.NP,test=False)
et=elapsed()


# ### Save and plot light-curve

# In[8]:


light_curve=pd.DataFrame()
ts=array(ts)
isort=ts.argsort()
light_curve["#Time (BKJD)"]=S.tbkjd+ts[isort]/DAY
light_curve["Normalized PDCSAP_Flux"]=array(s)[isort]
ds=array(ds)[isort]


# In[9]:


light_curve.to_csv(f"{CONF.SYSTEM}.txt",index=False,header=True)


# In[10]:


fig=plt.figure()
ax=fig.gca()
ax.plot(light_curve["#Time (BKJD)"],light_curve["Normalized PDCSAP_Flux"],'k.')
fig.savefig(f"{CONF.SYSTEM}.png")


# In[ ]:




