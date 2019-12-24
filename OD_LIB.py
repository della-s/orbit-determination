#!/usr/bin/env python
# coding: utf-8

# In[2]:


def dot(Vector1, Vector2): #dot product
    
    
    dot = (Vector1[0]*Vector2[0]) + (Vector1[1]*Vector2[1]) + (Vector1[2]*Vector2[2])

    print dot
    
Vector1 = [ 0.33616722, -1.61421096,  0.12680971]
Vector2 = [ 0.47734576, -0.19738789,  0.38962745]
dot(Vector1, Vector2) #expected result 32


# In[37]:


def cross(Vector1, Vector2): #cross product
    
    
    return (((Vector1[1]*Vector2[2])-(Vector1[2]*Vector2[1])), ((Vector1[2]*Vector2[0])-(Vector1[0]*Vector2[2])), ((Vector1[0]*Vector2[1])-(Vector1[1]*Vector2[0])))
    

cross([1, 2, 3], [4, 5, 6]) #expected result (-3, 6, -3)


# In[41]:


def triple(Vector1, Vector2, Vector3): #triple product A dot (B x C)
    
    a = cross(Vector2,Vector3)
    d = dot(Vector1,a)
    
    return d
    
triple ([1, 2, 3], [4, 5, 6], [7, 8, 9]) #expected result 0


# In[19]:


def HMStoDeg(h, m, s): #RA HMS to Degrees
    
    RAdeg = 15.0*h + (15.0*m)/60.0 + (15.0*s)/3600.0
    
    return RAdeg

HMStoDeg (12, 0, 0) #expected result 180.0


# In[23]:


def DMStoDec(deg, arcmin, arcsec): #DEC DMS to Decimals
    
    dec = abs(deg) + abs(arcmin/60.0) + abs(arcsec/3600.0)
    
    if(deg) < 0:  
        dec = -dec
        
    #if deg = 0
    
    if(arcmin) < 0:
        dec = -dec
        
    if(arcsec) < 0:
        dec = -dec
    
    return dec

#input format: put sign only in front of first number, eg. ((+/-)x, y, z) 

DMStoDec (-10, -11, 12) #expected result 10.1867


# In[1]:


def DegtoHMS(decimal): #RA Decimals to HMS
    
    a = decimal/15.0
    hours = int(a)
    b = a - hours
    c = b*60
    minutes = int(c)
    d = c - minutes
    seconds = d*60
    
    return [hours, minutes, seconds]

DegtoHMS (222.125077) #expected result (12, 0, 0)


# In[18]:


def DectoDMS(decimal): #DEC Decimals to DMS
    
    degrees = int(decimal)
    v = decimal - degrees
    x = v*60.0
    decmin = int(x)
    y = x - decmin
    z = y*60.0
    decsec = int(z)
    
    return [degrees, decmin, decsec]

DectoDMS(10.186666666666667) #expected result (10, 11, 12)


# In[3]:


from math import *
def magnitude(Vector1): #Magnitude
    
    mag = sqrt(Vector1[0]**2 + Vector1[1]**2 + Vector1[2]**2)
    
    return mag

Vector1 = [ 0.33616722, -1.61421096,  0.12680971]
magnitude(Vector1) #expected result 3.74


# In[ ]:




