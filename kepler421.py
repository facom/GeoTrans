from geotrans import *

System=dict2obj(dict(\

#########################################
#SYSTEM PRIMARY PARAMETERS
#########################################

#//////////////////////////////
#DETECTOR
#//////////////////////////////

#Aperture, m
Ddet=0.95,

#QuantumEfficiency
QDetector=QKEPLER,

#//////////////////////////////
#STAR
#//////////////////////////////

#Name
StarName="Kepler-421",
CatName="KIC8800954",

#Kepler data source
#Search: https://archive.stsci.edu/kepler/data_search/search.php
DataSource="https://archive.stsci.edu/kepler/preview.php?type=lc&dsn=KPLR008800954-2009131105131",

#Physical properties
Mstar=0.794*MSUN,
Rstar=0.819*RSUN,
Tstar=5141, #K, MAST Archive
Dstar=0.352100278159*KILO*PARSEC, #Simbad
Lstar=-1, #Compute using R and T

#Position
RA="18:53:1.6", #hours
Dec="+45:05:16", #degrees

#Limb darkening parameters
c1=0.70,#Limb Darkening
c2=-0.24,#Limb Darkening

#//////////////////////////////
#ORBIT
#//////////////////////////////

#Orbital properties
ap=1.0*AU,
ep=0.0,

#Time of passage by periapsis
tp="2009-05-22T00:19:58", #UTC date
tbkjd=140.6540, #BKJD epoch

#Orbit inclination with respect to the sky plane
iorb=89.95*DEG,

#Argument of the periapsis
wp=0.0*DEG,

#//////////////////////////////
#PLANET
#//////////////////////////////

#Properties of the planet
Mplanet=1.0*MSAT,
Rplanet=1.0*RSAT,

#Geometrical properties
fp=0.0, #Oblateness

#//////////////////////////////
#RINGS
#//////////////////////////////

#Dimensions
fe=RSAT_ARING/RSAT, #Exterior ring (Rp)
fi=RSAT_BRING/RSAT, #Interior ring (Rp)

#Orientation
ir=30.0*DEG, #Ring inclination
phir=60.0*DEG, #Ring roll angle

#Opacity
tau=1.0, #Opacity

))

#########################################
#SYSTEM DERIVATIVE PARAMETERS
#########################################
derivedSystemProperties(System)
updatePlanetRings(System)
updatePosition(System,System.tcen)
