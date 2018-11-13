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
StarName="Sun",
CatName="Sun",

#Kepler data source
#Search: https://archive.stsci.edu/kepler/data_search/search.php
DataSource="https://archive.stsci.edu/kepler/data_search/search.php",

#Physical properties
Mstar=1.0*MSUN,
Rstar=1.0*RSUN,
Lstar=1.0*LSUN,
Tstar=1.0*TSUN,
Dstar=1*KILO*PARSEC,

#Position
RA="00:00:00.00", #hours
Dec="+00:00:00.0", #degrees

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
tp="2018-01-01T00:00:00", #UTC date
tbkjd=3286.499785155058, #BKJD date

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
