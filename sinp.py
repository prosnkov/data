#!/usr/bin/python
# -*- coding: utf-8 -*-
#Different mini programs for my science work
from __future__ import unicode_literals
import os 
import subprocess
import numpy as np
import requests
import re
import math
import scipy.interpolate as interpolate
import scipy.integrate as integrate
import shutil
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rc('text', usetex=True)
mpl.rcdefaults()
mpl.font = mpl.font_manager.FontProperties(fname="/usr/share/matplotlib/mpl-data/fonts/ttf/cmb10.ttf")

#initialise talys, added Apr 2017
def talys_init(element, mass, Emin, Emax, dE):
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
	os.chdir('..')
	return 0

#take cross sections from talys, added Mar 2018
def xs_talys(element, mass, prot, neut, param):
	if os.path.isfile('rp0'+str(element-prot)+str(mass-neut-prot)+'.tot'):
		lines = open('rp0'+str(element-prot)+str(mass-neut-prot)+'.tot').readlines()			
	if os.path.isfile('rp0'+str(element-prot)+'0'+str(mass-neut-prot)+'.tot'):
		lines = open('rp0'+str(element-prot)+'0'+str(mass-neut-prot)+'.tot').readlines()
	for s in range(5):
		del lines[0]
	table = np.array([row.split() for row in lines])
	x = np.float64(np.delete(table,1,1))
	y = np.float64(np.delete(table,0,1))
	if param == 'x':
		return x
	else:
		return y

#interpolate data, added Apr 2018
def interp_val(arr_x, arr_y, dE, param):
	x = np.arange(arr_x[0], arr_x[-1]+dE, dE)
	f = interpolate.InterpolatedUnivariateSpline(arr_x, arr_y)
	y = f(x)
	if param == 'x':
		return x
	if param == 'y':
		return y
def interp_arr(arr_x, arr_y, new_arr):
	f = interpolate.InterpolatedUnivariateSpline(arr_x, arr_y)
	y = f(new_arr)
	return y

#parsing data from EXFOR database at cdfe.sinp.msu.ru, added Nov 2018
def parse(link,param):
	r = requests.get(str(link))
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
	if param == 'x':
		return x
	if param == 'y':
		return y
	if param == 'err':
		return err


#make one plot, added Nov 2018
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
	plt.xlim(0,)
	plt.ylim(0,)
	default_size = fig.get_size_inches() 
	fig.set_size_inches(default_size[0]*1.5, default_size[1]*1.5,forward=True)
	plt.savefig(name+'.eps',fmt='eps',rasterized=False)
	return fig
