#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
surface_utils.py
'''

import numpy as np

def rotation_matrix(theta, phi):
	return np.array([[np.cos(theta) * np.cos(phi), -np.sin(phi), np.sin(theta) * np.cos(phi)],
					 [np.cos(theta) * np.sin(phi), np.cos(phi), np.sin(theta) * np.sin(phi)],
					 [-np.sin(theta), 0., np.cos(theta)]])
					 
def latlon2thetaphi(lat, lon, flattening):
	temp = (1. - flattening) * (1. - flattening)
	return np.pi / 2. - np.arctan(temp * np.tan(np.radians(lat))), np.radians(lon)

def thetaphi2latlon(theta, phi, flattening):
	temp = (1. - flattening) * (1. - flattening)
	return np.degrees(np.arctan(temp * np.tan(np.pi / 2. - theta))), np.degrees(phi)
		   
def thetaphi2xyz(theta, phi):
	return np.array([np.sin(theta) * np.cos(phi),
					 np.sin(theta) * np.sin(phi),
		 			 np.cos(theta)])

def xyz2thetaphi(xyz):
	theta = np.arccos(xyz[2])
	phi = np.arctan2(xyz[1], xyz[0])
	return theta, phi 

def interpLagrange(target, bases):
	nbases = len(bases)
	results = np.zeros(nbases)
	for dgr in np.arange(0, nbases):
		prod1 = 1.
		prod2 = 1.
		for i in np.arange(0, nbases):
			if i != dgr:
				prod1 *= target - bases[i]
				prod2 *= bases[dgr] - bases[i]
		results[dgr] = prod1 / prod2
	return results	
	
class SurfaceStation:
	def __init__(self, network, name):
		self.network = network
		self.name = name

	def setloc_geographic(self, lat, lon, flattening, srclat, srclon, srcflattening):
		self.lat = lat
		self.lon = lon
		srctheta, srcphi = latlon2thetaphi(srclat, srclon, srcflattening)
		rmat = rotation_matrix(srctheta, srcphi)
		theta, phi = latlon2thetaphi(lat, lon, flattening)
		xglb = thetaphi2xyz(theta, phi)
		xsrc = rmat.T.dot(xglb)
		self.dist, self.azimuth = xyz2thetaphi(xsrc)
	
	def setloc_source_centered(self, dist, azimuth, flattening, srclat, srclon, srcflattening):
		self.dist = dist
		self.azimuth = azimuth
		srctheta, srcphi = latlon2thetaphi(srclat, srclon, srcflattening)
		rmat = rotation_matrix(srctheta, srcphi)
		xsrc = thetaphi2xyz(dist, azimuth)
		xglb = rmat.dot(xsrc)
		theta, phi = xyz2thetaphi(xglb)
		self.lat, self.lon = thetaphi2latlon(theta, phi, flattening)
		
		