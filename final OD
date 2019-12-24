

from math import *
from __future__ import division
import numpy as np

#input values
R1 = [-2.425884330693873E-01, 9.058984039584287E-01, 3.926679094006759E-01]
R2 = [-2.915566543679477E-01, 8.936487990469925E-01, 3.873574976180230E-01]
R3 = [-3.558343435190608E-01, 8.737323629744292E-01, 3.787253397394694E-01]
RA1 = [14, 45, 19.27]
RA2 = [14, 48, 52.49]
RA3 = [14, 54, 34.98]
DEC1 = [6, 38, 38.1]
DEC2 = [5, 16, 14]
DEC3 = [3, 16, 19.4]
t1 = 2458305.663194444
t2 = 2458308.663194
t3 = 2458312.672222222

#other stuff
k = 0.01720209895
oblq = (pi/180)*(23.4368812082)
c = 299792458

#-----FUNCTIONS-----
#dot product
def dot(Vector1, Vector2):
    
    
    dot = (Vector1[0]*Vector2[0]) + (Vector1[1]*Vector2[1]) + (Vector1[2]*Vector2[2])

    return dot

#magnitude
def mag(Vector):
    
    mag = sqrt(Vector[0]**2 + Vector[1]**2 + Vector[2]**2)
    
    return mag

#quad solver
def sincos(sin, cos):
    theta = atan(sin/cos) 
    if(sin and cos >0):
        quad = 1
    if (sin>0 and cos<0):
        theta = theta + (0.5*pi)
        quad = 2
    if (sin<0 and cos<0):
        theta = theta + pi
        quad = 3
    if (sin<0 and cos>0):
        theta = theta + (1.5*pi)
        quad = 4
    return [theta, quad]

#RA HMS to degrees
def HMStoDeg(h, m, s):
    RA = 15.0*h + (15.0*m)/60.0 + (15.0*s)/3600.0
    return RA

#DEC DMS to decimals
def DMStoDec(deg, arcmin, arcsec): 
    dec = abs(deg) + abs(arcmin/60.0) + abs(arcsec/3600.0) 
    if(deg) < 0:
        dec = -dec
    #if deg = 0
    if(arcmin) < 0:
        dec = -dec
    if(arcsec) < 0:
        dec = -dec
    return dec

#RA Decimals to HMS
def DegtoHMS(decimal): 
    
    a = decimal/15.0
    hours = int(a)
    b = a - hours
    c = b*60
    minutes = int(c)
    d = c - minutes
    seconds = d*60
    
    return [hours, minutes, seconds]

#DEC Decimals to DMS
def DectoDMS(decimal):
    
    degrees = int(decimal)
    v = decimal - degrees
    x = v*60.0
    decmin = int(x)
    y = x - decmin
    decsec = y*60.0
    
    return [degrees, decmin, decsec]

#-----THE OD-----

ra1 = HMStoDeg(RA1[0], RA1[1], RA1[2]) *pi/180. 
ra2 = HMStoDeg(RA2[0], RA2[1], RA2[2]) *pi/180. 
ra3 = HMStoDeg(RA3[0], RA3[1], RA3[2]) *pi/180.

dec1 = DMStoDec(DEC1[0], DEC1[1], DEC1[2])* pi/180. 
dec2 = DMStoDec(DEC2[0], DEC2[1], DEC2[2])* pi/180. 
dec3 = (DMStoDec(DEC3[0], DEC3[1], DEC3[2])) * pi/180.

#print ra1, ra2, ra3, dec1, dec2, dec3

#RA and Dec to cartesian
def rhohat(RA, dec): 
    x = cos(RA)*cos(dec) 
    y = sin(RA)*cos(dec) 
    z = sin(dec) 
    rhohat = [x, y, z] 
    return rhohat

rhohat1 = rhohat(ra1, dec1) 
rhohat2 = rhohat(ra2, dec2) 
rhohat3 = rhohat(ra3, dec3) 
print "rhohat:", rhohat1, '\n', rhohat2, '\n', rhohat3, '\n'

#Gaussian time intervals
T0 = k*(t3-t1)
T1 = k*(t1-t2)
T3 = k*(t3-t2)
print "T0:", T0, "T1:", T1, "T3:", T3

#D0
D0 = np.dot(rhohat1, np.cross(rhohat2, rhohat3))
if D0 == 0:
    print "error"
else:
    print "D0:", D0

#initial estimates of a1 and a3
a1est = T3/T0
a3est = -T1/T0
print "a1 and a3 est.:", a1est, a3est

#initial estimates of Dijs
D11 = np.dot(np.cross(R1, rhohat2), rhohat3)
D12 = np.dot(np.cross(R2, rhohat2), rhohat3)
D13 = np.dot(np.cross(R3, rhohat2), rhohat3)

D21 = np.dot(np.cross(rhohat1, R1), rhohat3)
D22 = np.dot(np.cross(rhohat1, R2), rhohat3)
D23 = np.dot(np.cross(rhohat1, R3), rhohat3)

D31 = np.dot(np.cross(rhohat2, R1), rhohat1)
D32 = np.dot(np.cross(rhohat2, R2), rhohat1)
D33 = np.dot(np.cross(rhohat2, R3), rhohat1)
D30 = np.dot(np.cross(rhohat2, rhohat3), rhohat1)

#initial estimates of rho
rho1est = ((a1est*D11)-D12+(a3est*D13))/(a1est*D0)
rho2est = ((a1est*D21)-D22+(a3est*D23))/(-1*D0)
rho3est = ((a1est*D31)-D32+(a3est*D33))/(a3est*D30)
print "est. rho:", rho1est, rho2est, rho3est

#r1 and r3
r1 = rho1est*np.array(rhohat1)-np.array(R1)
r3 = rho3est*np.array(rhohat3)-np.array(R3)
print "r1 and r3:", r1, r3

#r2
r2old = rho2est*np.array(rhohat2)-np.array(R2)
print "r2old:", r2old

#r2dot
r2dotold = (r3-r1)/T0
print "r2dotold:", r2dotold

#iteration loop
niter = 0
nmax = 1000
tiny = 1E-9

diff1 = 1000
diff2 = 1000


f1= 1- (T1)**2/(2*np.linalg.norm(r2old)**3)+dot(r2old, r2dotold)*T1**3/(2*np.linalg.norm(r2old)**5)+T1**4/(24*np.linalg.norm(r2old)**3)*(3*(dot(r2dotold, r2dotold)/(np.linalg.norm(r2old))**2-1/(np.linalg.norm(r2old)**3))-15*(dot(r2old, r2dotold)/(np.linalg.norm(r2old))**2)**2+1/(np.linalg.norm(r2old)**3))
f3= 1- (T3)**2/(2*np.linalg.norm(r2old)**3)+dot(r2old, r2dotold)*T3**3/(2*np.linalg.norm(r2old)**5)+T3**4/(24*np.linalg.norm(r2old)**3)*(3*(dot(r2dotold, r2dotold)/(np.linalg.norm(r2old))**2-1/(np.linalg.norm(r2old)**3))-15*(dot(r2old, r2dotold)/(np.linalg.norm(r2old))**2)**2+1/(np.linalg.norm(r2old)**3))
g1= T1-T1**3/(6*np.linalg.norm(r2old)**3)+dot(r2old, r2dotold)*T1**4/(4*np.linalg.norm(r2old)**5)
g3= T3-T3**3/(6*np.linalg.norm(r2old)**3)+dot(r2old, r2dotold)*T3**4/(4*np.linalg.norm(r2old)**5)


print "f1:", f1, "f3:", f3, "g1:", g1, "g3:", g3

a1 = g3/((f1*g3)-(f3*g1))
a3 = -1*g1/((f1*g3)-(f3*g1))


print "a1:", a1, "a3:", a3

r2 = (g3*r1-g1*r3)/(f1*g3-f3*g1)
r2dot = (f3*r1-f1*r3)/(f3*g1-f1*g3)

while diff1 > tiny and diff2 > tiny:
    niter = niter + 1

    r2old = r2
    r2dotold = r2dot

    T0 = k*((t3-(rho3est/c)-(t1-(rho1est/c))))
    T1 = k*((t1-(rho1est/c)-(t2-(rho2est/c))))
    T3 = k*((t3-(rho3est/c)-(t2-(rho2est/c))))
    
    rho1 = ((a1*D11)-D12+(a3*D13))/(a1*D0)
    rho2 = ((a1*D21)-D22+(a3*D23))/(-1*D0)
    rho3 = ((a1*D31)-D32+(a3*D33))/(a3*D30)
        
    r1 = rho1*np.array(rhohat1)-np.array(R1)
    r2 = rho2*np.array(rhohat2)-np.array(R2)
    r3 = rho3*np.array(rhohat3)-np.array(R3)

    r2dot = np.array((f3*r1-f1*r3)/(f3*g1-f1*g3))
    r2 = np.array((g3*r1-g1*r3)/(f1*g3-f3*g1))
    #print r2, r2dot

    f1= 1- (T1)**2/(2*np.linalg.norm(r2)**3)+dot(r2, r2dot)*T1**3/(2*np.linalg.norm(r2)**5)+T1**4/(24*np.linalg.norm(r2)**3)*(3*(dot(r2dot, r2dot)/(np.linalg.norm(r2))**2-1/(np.linalg.norm(r2)**3))-15*(dot(r2, r2dot)/(np.linalg.norm(r2))**2)**2+1/(np.linalg.norm(r2)**3))
    f3= 1- (T3)**2/(2*np.linalg.norm(r2)**3)+dot(r2, r2dot)*T3**3/(2*np.linalg.norm(r2)**5)+T3**4/(24*np.linalg.norm(r2)**3)*(3*(dot(r2dot, r2dot)/(np.linalg.norm(r2))**2-1/(np.linalg.norm(r2)**3))-15*(dot(r2, r2dot)/(np.linalg.norm(r2))**2)**2+1/(np.linalg.norm(r2)**3))
    g1= T1-T1**3/(6*np.linalg.norm(r2)**3)+dot(r2, r2dot)*T1**4/(4*np.linalg.norm(r2)**5)
    g3= T3-T3**3/(6*np.linalg.norm(r2)**3)+dot(r2, r2dot)*T3**4/(4*np.linalg.norm(r2)**5)

    a1 = g3/((f1*g3)-(f3*g1))
    a3 = -1*g1/((f1*g3)-(f3*g1))

    diff1 = mag(r2 - r2old)
    diff2 = mag(r2dot - r2dotold)
    niter = niter +1
    
print niter, "r2:", r2, "r2dot:", r2dot

#convert to ecliptic
m = [[1., 0., 0.], [0., cos(oblq), sin(oblq)], [0., -sin(oblq), cos(oblq)]]
ecr = np.matmul(m, r2)
ecrdot = np.matmul(m, r2dot)
#print "ecliptic r:", ecr

#magnitude of r
rmag = mag(ecr)
#print "r mag:", rmag

#h and magnitude
h = np.cross(ecr, ecrdot)
hmag = mag(h)
#print "h:", h, "h mag:", hmag

#v squared
v2 = np.dot(ecrdot, ecrdot)
#print "v2:", v2

#a
a = 1./((2./rmag)-v2)
print "a:", a

#e
e = sqrt(1.-((hmag**2.)/a))
print "e:", e

#i
i = acos(abs(h[2])/hmag)
i = i*(360./(2*pi))
print "i:", i

#omega
coso = -(h[1]/(hmag*sin(i)))
sino = h[0]/(hmag*sin(i))
o = sincos(sino, coso)
o1 = o[0]*(180/pi) + 90
print "o:", o1


#u
cosu = ((ecr[0]*coso) + (ecr[1]*sino))/rmag
sinu = ecr[2]/(rmag*sin(i))
u = sincos(sinu, cosu)
u1 = u[0]*(360/(2*pi))
#print "u:", u1

#v
inve = 1./e
ae = a*(1.-(e**2))
dot = np.dot(ecr, ecrdot)
cosv = inve*((ae/rmag)-1)
sinv = inve*((ae/hmag)*(dot/rmag))
v = sincos(sinv, cosv)
v1 = v[0]*(360/(2*pi)) + 90
#print "v:", v1, v[1]

    
#w
w = u1 - v1 + 360 + 90
print "w:", w

#E
E1 = acos(inve*(1-(rmag/a)))
E = E1*(180/pi)
#print "E:", E

#M
M = 270 - (E - e*sin(E))
print "M:", M

#T
n = sqrt(1./(a**3))
T = t2 - (M/(n*k))
#print "T:", T

#expected values
jpla = 1.863446788208271
jple = 0.396828542104546
jpli = 8.448189874886131
jplo = 133.1919441887713
jplw = 167.6468863346789
jplM = 256.4687547825516

eerror = ((e-jple)/jple)*100.
oerror = ((o1-jplo)/jplo)*100.
aerror = ((a-jpla)/jpla)*100.
werror = ((w-jplw)/jplw)*100.
Merror = ((M-jplM)/jplM)*100.
ierror = ((i-jpli)/jpli)*100.

print "expected a:", jpla, "% error:", aerror
print "expected e:", jple, "% error:", eerror
print "expected i:", jpli, "% error:", ierror
print "expected o:", jplo, "% error:", oerror
print "expected w:", jplw, "% error:", werror
print "expected M:", jplM, "% error:", Merror




