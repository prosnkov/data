#!/usr/bin/python
# -*- coding: utf-8 -*-
#My module with some mini programs for my science work
#Upd: 29.11.2018
from __future__ import unicode_literals
import os 
import subprocess
import numpy as np
import requests
import re
import math
import scipy.interpolate as interpolate
import scipy.integrate as integrate
from collections import namedtuple
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rc('text', usetex=True)
mpl.rcdefaults()
mpl.font = mpl.font_manager.FontProperties(fname="/usr/share/matplotlib/mpl-data/fonts/ttf/cmb10.ttf")

nucl=['neut','H','He','Li','Be','B','C','N','O','F','Ne','Na', \
	'Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti', \
	'V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge', \
	'As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc', \
	'Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe', \
	'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd', \
	'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re', \
	'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At', \
	'Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm', \
	'Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg', \
	'Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og']

#initialise talys, added Apr 2017
# example:
#	talys_init(Z,A,Emin=5,Emax=40,dE=0.5)
def talys_init(element, mass, Emin=1, Emax=55, dE=1):
#	d=os.getcwd() 
#	os.chdir(d)
#	os.mkdir(str(mass))
#	os.chdir(str(mass))
	f0 = open('input','w')
	f0.write('projectile g'+'\n')
	f0.write('element '+str(element)+'\n')
	f0.write('mass '+str(mass)+'\n')
	f0.write('energy '+str(Emin)+' '+str(Emax)+' '+str(dE)+'\n')
	f0.write('ejectiles n p'+'\n')
	f0.close()
	subprocess.call("talys<input>output", shell=True)
#	os.chdir('..')
	return 0

#take cross sections from talys, added Mar 2018
#element - charge number
#mass - mass number
#prot - quantity of emitted protons
#neut - quantity of emitter neutrons
# example:
#	E=sinp.xs_talys(Z,A,0,1).x
#	xs=sinp.xs_talys(Z,A,0,1).y
def xs_talys(element, mass, prot, neut):
	if os.path.isfile('rp0'+str(element-prot)+str(mass-neut-prot)+'.tot'):
		lines = open('rp0'+str(element-prot)+str(mass-neut-prot)+'.tot').readlines()			
	if os.path.isfile('rp0'+str(element-prot)+'0'+str(mass-neut-prot)+'.tot'):
		lines = open('rp0'+str(element-prot)+'0'+str(mass-neut-prot)+'.tot').readlines()
	for s in range(5):
		del lines[0]
	table = np.array([row.split() for row in lines])
	x = np.float64(np.delete(table,1,1))
	y = np.float64(np.delete(table,0,1))
	result = namedtuple('arrays', ['x','y'])
	return result(x,y)

#interpolate data, added Apr 2018
#interp_val must be used for changing dE with the same boundaries
#example:
#	E_int=sinp.interp_val(E,xs,dE=0.001).x
#	xs_int=sinp.interp_val(E,xs,dE=0.001).y
#interp_arr must be used for changing array of energies
#example:
#	xs_int=sinp.interp_arr(E,xs,E_int)
def interp_val(arr_x, arr_y, dE=0.01):
	x = np.arange(arr_x[0], arr_x[-1]+dE, dE)
	f = interpolate.InterpolatedUnivariateSpline(arr_x, arr_y)
	y = f(x)
	result = namedtuple('arrays', ['x','y'])
	return result(x,y)
def interp_arr(x, y, x_new):
	f = interpolate.InterpolatedUnivariateSpline(x, y)
	y = f(x_new)
	return y

#parsing data from EXFOR database at cdfe.sinp.msu.ru, added Nov 2018
# example:
#	E=sinp.exfor(link).x
#	xs=sinp.exfor(link).y
#	err=sinp.exfor(link).err
def exfor(link):
	r = requests.get(link)
	lines = r.text.split('\n')
	i=0
	for line in lines:
		if line.startswith('MEV        MB         MB') == True:
			break
		i=i+1
	arr = r.text.split('\n')[i+1:-3]
	arr = [re.sub(r'\s+', r' ', elem.lstrip().rstrip()).split() for elem in arr]
	table = []
	for elem in arr:
		if len(elem) == 5:
			table.append(elem[:2] + [0])
		
		else:
			table.append(elem[:3])
	x = np.float64(np.delete(table,[1,2],1)).ravel()
	y = np.float64(np.delete(table,[0,2],1)).ravel()
	err = np.float64(np.delete(table,[0,1],1)).ravel()
	result = namedtuple('arrays', ['x','y','err'])
	return result(x,y,err)
#similar to exfor, added Nov 2018
def endf(link):
	r = requests.get(link)
	lines = r.text.split('\n')
	i=0
	for line in lines:
		if line.startswith('#E,eV        Sig,b        ') == True:
			break
		i=i+1
	arr = r.text.split('\n')[i+1:-2]
	arr = [re.sub(r'\s+', r' ', elem.lstrip().rstrip()).split() for elem in arr]
	table = []
	for elem in arr:
		table.append(elem[:2])
	x = np.float64(np.delete(table,1,1)).ravel()
	y = np.float64(np.delete(table,0,1)).ravel()
	x = np.multiply(x,math.pow(10,-6))
	y = np.multiply(y,math.pow(10,3))
	result = namedtuple('arrays', ['x','y'])
	return result(x,y)


#make one plot, added Nov 2018
#Can be used for plotting cross sections
#Beta version
def plot(x,y,err,lab,name,param):
	fig = plt.figure()
	fig.set_rasterized(False)
	if param == 'exp':
		plt.scatter(x,y,c='k',s=20,marker='o',label=lab)
		plt.errorbar(x, y, yerr=err, ls='none',c='k',elinewidth=1) 
	if param == 'th':
		plt.plot(x,y,c='k',linewidth=2,label=lab)
	plt.tight_layout()
	plt.legend()
	plt.xlabel('$E_{\gamma}$, МэВ')
	plt.ylabel('$\sigma$, мб')
	plt.xlim(5,)
	plt.ylim(0,)
	default_size = fig.get_size_inches() 
	fig.set_size_inches(default_size[0]*1.5, default_size[1]*1.5,forward=True)
	plt.savefig(name+'.eps',fmt='eps',rasterized=False)
	return fig
