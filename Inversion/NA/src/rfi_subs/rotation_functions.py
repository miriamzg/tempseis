
import numpy as np
import sys

def lonlatdep_2_xy(lon, lat, depth, rc, dip, strike):
	R = 6371.
	L = (2*np.pi*R)/360.
	dip = np.deg2rad(dip)
	strike = np.deg2rad(strike)	
	teta = (np.pi/2.)-strike

	rc_lon = rc[0]
	rc_lat = rc[1]
	rc_depth = rc[2]
	x = lon - rc_lon
	y = lat - rc_lat
	x = x * (L*np.cos(np.deg2rad(rc_lat)))
	y = y * L
	z = depth

	R1 = [[np.cos(teta), -np.sin(teta),0],[np.sin(teta),np.cos(teta), 0],[0,0,1]]
	r = [x,y,z]
	r1 = np.dot(r,R1)
	R2 = [[1,0,0],[0,np.cos(-dip), -np.sin(-dip)],[0,np.sin(-dip),np.cos(-dip)]]
	r2 = np.dot(r1, R2)
	r2[2] = r2[2] - rc_depth

	xrot, yrot, zrot = r2[0], r2[1], r2[2]

	return xrot, yrot, zrot


def xy_2_lonlatdep(x, y, z, rc, dip, strike):
	R = 6371.
	L = (2*np.pi*R)/360.
	dip = np.deg2rad(dip)
	strike = np.deg2rad(strike)	
	
	rc_lon = rc[0] 
	rc_lat = rc[1]
	rc_depth = rc[2]

	r = [x,y,z]
	R1 = [[1,0,0],[0,np.cos(dip), -np.sin(dip)],[0,np.sin(dip),np.cos(dip)]]
	r1 = np.dot(r, R1)

	teta = (np.pi/2.)-strike
	R2 = [[np.cos(-teta), -np.sin(-teta),0],[np.sin(-teta),np.cos(-teta), 0],[0,0,1]]
	r2 = np.dot(r1, R2)

	x = r2[0]
	y = r2[1]
	z = r2[2]


	#print dep

	lon = x / (L*np.cos(np.deg2rad(rc_lat))) 
	lat = y / L
	lon = lon + rc_lon
	lat = lat + rc_lat
	dep = z + rc_depth

	return lon, lat, dep