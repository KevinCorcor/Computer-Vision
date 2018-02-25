'''Complete CS410 Assignment 07/11/2017

DRIVER CODE:
	To run all three desired functions:> run CalibrateCamera.py

author: Kevin Corcoran - 14301776
'''
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d

def calibrateCamera3D(data):											#outputs Perspective Projection Matrix

	XYZ1 = np.c_[data[:,[0,1,2]],np.ones(491)] 							#concatenate XYZ with 1

	L = np.zeros((982,4)) 												#Left 4 columns of A of Ax=y
	L[0::2] = XYZ1  													#every second row gets initialised with XYZ1 from first row

	M = np.zeros((982,4)) 												#Middle 4 columns of A of Ax=y
	M[1::2] = XYZ1 														#every second row gets initialised with XYZ1 from second row

	bx = (-1) * np.array(data[:,3])										#-x coords
	cx = XYZ1 * bx[:,None]												#multiply XYZ1 by -x

	by = (-1) * np.array(data[:,4])										#-y coords
	cy = XYZ1 * by[:,None]												#multiply XYZ1 by -y

	R = np.zeros((982,4))												#Right hand side of Ax=y
	R[0::2] = cx														#every second row gets initialised with -x from first row
	R[1::2] = cy														#every second row gets initialised with -y from second row

	LMR = np.asmatrix(np.c_[ np.c_[ L, M ], R ])						#concatenated left, middle and right 4 columns

	eigenvals, eigenvectors  = np.linalg.eig(np.asmatrix(LMR.T)*LMR)	#calculate eigenvalues and eigenvectors

	return eigenvectors[:,eigenvals.argmin()].reshape(3,4)				#return the eigenvector with the lowest corresponding eigenvalue


def visualiseCameraCalibration3D(data, P):								#generate the actual 2D image on top of the reprojection from the perspective projective matrix and actual 3D points

	XYZ1_T = (np.c_[np.asmatrix(data[:,[0,1,2]]),np.ones(491)]).T 		#traspose of concatenate XYZ with 1
	xyw = (P *	XYZ1_T).T 												#transpose of the perspective projection matrix multiplied by XYZ1_Transposed to produce xyw

	fig = plt.figure()
	ax = fig.gca()

	ax.plot(data[:,3], data[:,4],'gx', label = '2D image')
	ax.plot(xyw[:,0]/xyw[:,2], xyw[:,1]/xyw[:,2],'b+',label = 'Reprojection')#plot x' = x/w, y' = y/w

	fig.legend(loc='upper center')
	plt.ion()
	plt.show()															#plot actual image points


def evaluateCameraCalibration3D(data, P):								#gemerate reprojection statistics regarding distance of corresponding points

	XYZ1_T = (np.c_[np.asmatrix(data[:,[0,1,2]]),np.ones(491)]).T 		#transpose of XYZ concatenated with 1
	xyw = (P *	XYZ1_T).T 												#transpose of the perspective projection matrix multiplied by XYZ1_Transposed to produce xyw

	x = np.array(xyw[:,0]/xyw[:,2]) 									#reprojected imgage x coordinates
	y = np.array(xyw[:,1]/xyw[:,2]) 									#reprojected imgage y coordinates
	d = np.sqrt((x[:,0] - data[:,3])**2 + (y[:,0] - data[:,4])**2)  	#calculate the euclidean distance between corresponding points

	print "\nDistance Statistics: \n---------------------------------",\
	"\nMean\t=\t", np.mean(d),\
	"\nVariance=\t", np.var(d),\
	"\nMinimum\t=\t", d.min(),\
	"\nMaximum\t=\t", d.max()


data = np.loadtxt('data.txt')

P = calibrateCamera3D(data)

visualiseCameraCalibration3D(data, P)

evaluateCameraCalibration3D(data, P)
