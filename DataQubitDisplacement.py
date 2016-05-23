from math import sin, cos, sqrt, pi
from numpy import array, random

def random_uniform_cylinder(radius,half_height):
    '''
    Generates a uniform random offset within a cylinder along the z-axis
    of given radius and half-height
    '''
    r = random.uniform(0,radius)
    theta = random.uniform(0,2*pi)
    z = random.uniform(-half_height,half_height)
    
    x = r * cos(theta)
    y = r * sin(theta)
    
    return array([x,y,z])
    
def random_gaussian_cylinder(sd_radius,sd_height,convert=False):
    '''
    Generates a gaussian random offset, a random radius from the z-axis and
    a random distance along the z-axis
    
    Can optionally convert from uniform values such that the half-probability
    distance for the uniform cylinder corresponds to the half-probability distance
    for the gaussian distribution
    '''
    if convert == True: # A sort of conversion from radius and half-height to standard deviation
        sd_radius = sd_radius/2.388
        sd_height = sd_height/2.388
        
    r = random.gauss(0,sd_radius)
    theta = random.uniform(0,2*pi)
    z = random.gauss(0,sd_height)
    
    x = r * cos(theta)
    y = r * sin(theta)
    
    return array([x,y,z])