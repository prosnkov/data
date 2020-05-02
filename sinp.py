#!/usr/bin/python3
# -*- coding: utf-8 -*-
#My module with some mini programs for my science work
#Version: 0.3.9
#Upd: 02.05.2020
from __future__ import unicode_literals
import numpy as np
import requests
import re
import math
from collections import namedtuple
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rc('text', usetex=True)
mpl.rcdefaults()
mpl.font = mpl.font_manager.FontProperties(fname="/usr/share/matplotlib/mpl-data/fonts/ttf/cmb10.ttf")
mpl.rc('figure', max_open_warning = 0)

#array of periodic table names
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

#take cross sections from talys, added Mar 2018
#element - charge number
#mass - mass number
#prot - quantity of emitted protons
#neut - quantity of emitter neutrons
#example:
#	E=sinp.xstal(Z,A,0,1).x
#	xs=sinp.xstal(Z,A,0,1).y
def xstal(element, mass, prot, neut):
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

#parsing data from CDFE database at cdfe.sinp.msu.ru, added Nov 2018
#example:
#	E=sinp.cdfe(link).x
#	xs=sinp.cdfe(link).y
#	err=sinp.cdfe(link).err
def cdfe(link):
	if link.endswith('&SOURCE=ON') == False:
		link=link+'&SOURCE=ON'
	if link.startswith('http://') == False:
		link='http://'+link
	r = requests.get(link)
	lines = r.text.split('\n')
	i=0
	key=0
	col1=['KEV','MEV','GEV','ADEG','NO-DIM']
	col2=[' B',' MB',' MICRO-B',' B/SR',' MB/SR',' MU-B/SR',' NB/SR',' ARB-UNITS',' NO-DIM',' PC/FIS']
	col3=[' B',' MB',' MICRO-B',' B/SR',' MB/SR',' MU-B/SR',' NB/SR',' ARB-UNITS',' NO-DIM',' PC/FIS','']
	col4=[' NO-DIM',' MEV','']
	for a in range(len(col1)):
		i=0
		for b in range(len(col2)):
			i=0
			for c in range(len(col3)):
				i=0
				for d in range(len(col4)):
					i=0
					for line in lines:
						if " ".join(line.split()).startswith(f'{col1[a]}{col2[b]}{col3[c]}{col4[d]}') == True or \
						" ".join(line.split()).startswith(f'{col1[a]}{col4[d]}{col2[b]}{col3[c]}') == True:
							if d==0:key=1.1
							if d==1:key=2
							if d==2:key=1
							break
						i=i+1
					if key!=0: break
				if key!=0: break
			if key!=0: break
		if key!=0: break
	if key==1 or key==1.1:
		arr = r.text.split('\n')[i+1:-3]
		arr = [re.sub(r'\s+', r' ', elem.lstrip().rstrip()).split() for elem in arr]
		table = []
		if key==1:
			for elem in arr:
				if len(elem)%2 != 0:
					table.append(elem[:2] + [0])
				if len(elem)%2 == 0:
					table.append(elem[:3])
		if key==1.1:
			for elem in arr:
				if len(elem) == 4:
					table.append(elem[:2] + [0])
				if len(elem) >= 5:
					table.append(elem[:3])
		x = np.float64(np.delete(table,[1,2],1)).ravel()
		y = np.float64(np.delete(table,[0,2],1)).ravel()
		yerr=np.float64(np.delete(table,[0,1],1)).ravel()
		result = namedtuple('arrays', ['x','y','yerr','key'])
		key=1
		return result(x,y,yerr,key)
	if key==2:
		arr = r.text.split('\n')[i+1:-3]
		arr = [re.sub(r'\s+', r' ', elem.lstrip().rstrip()).split() for elem in arr]
		table = []
		for elem in arr:
				if len(elem) == 5:
					table.append(elem[:1] + [0] + elem[1:2] + [0])
				if len(elem) == 6:
					table.append(elem[:3] + [0])				
				if len(elem) >= 7:
					table.append(elem[:4])
		x = np.float64(np.delete(table,[1,2,3],1)).ravel()
		xerr = np.float64(np.delete(table,[0,2,3],1)).ravel()
		y = np.float64(np.delete(table,[0,1,3],1)).ravel()
		yerr = np.float64(np.delete(table,[0,1,2],1)).ravel()
		result = namedtuple('arrays', ['x','xerr','y','yerr','key'])
		return result(x,xerr,y,yerr,key)

#parsing data from ENDF database at nndc.bnl.gov, added Nov 2018
#	E=sinp.endf(link).x
#	xs=sinp.endf(link).y
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
	x = np.multiply(np.float64(np.delete(table,1,1)).ravel(),math.pow(10,-6))
	y = np.multiply(np.float64(np.delete(table,0,1)).ravel(),math.pow(10,3))
	result = namedtuple('arrays', ['x','y'])
	return result(x,y)

#parsing data from non-smoker database at nucastro.org, added Nov 2018
#example:
#	E=sinp.nonsmok(Z,A,'p').x
#	xs=sinp.nonsmok(Z,A,'p').y
#param = 'n', 'p', 'a'
def nonsmok(element,mass,param='n'):
	x=[]
	y=[]
	if param == 'n':
		link = 'http://download.nucastro.org/astro/photon/table1.asc'
	if param == 'p':
		link = 'http://download.nucastro.org/astro/photon/table2.asc'
	if param == 'a':
		link = 'http://download.nucastro.org/astro/photon/table3.asc'
	r = requests.get(link)
	lines = r.text.split('\n')
	i=0
	for line in lines:
		if line.startswith(nucl[element].lower()+str(mass)+'  ') == True or line.startswith(' '+nucl[element].lower()+str(mass)+'  ') == True:
			break
		i=i+1
	arr = r.text.split('\n')[i+1:i+20]
	arr = [re.sub(r'\s+', r' ', elem.lstrip().rstrip()).split() for elem in arr]
	for line in arr:
		for i in range(len(line)):
			if i % 2 == 0:
				x.append(line[i])
			else:
				y.append(line[i])
	x = np.multiply(np.float64(x),math.pow(10,-3))
	y = np.multiply(np.float64(y),math.pow(10,3))
	result = namedtuple('arrays', ['x','y'])
	return result(x,y)

#take cross sections from CMPR(Orlin), added Feb 2020
#inputs:
#prot - quantity of emitted protons
#neut - quantity of emitter neutrons
#outputs:
#xsmax - maximum of cross section
#Emax - energy at maximum of cross section
#Eg - energies of gamma
#xsg - absolute cross section
#T0 - T_< component
#T1 - T_> component
#quad - quadrupole component
#ober - obertone component
#deut - quasideutrone component
#example:
#	E=sinp.xsor(0,1).Eg
#	xs=sinp.xsor(0,1).xsg
def xsor(prot, neut):
	lines = open('key').readlines()
	table = np.array([row.split() for row in lines])
	k = np.intc(np.delete(table,[1,2,3,4],1)).ravel()
	sec = np.float64(np.delete(table,[0,2,3,4],1)).ravel()
	energ = np.float64(np.delete(table,[0,1,3,4],1)).ravel()
	p = np.intc(np.delete(table,[0,1,2,4],1)).ravel()
	n = np.intc(np.delete(table,[0,1,2,3],1)).ravel()
	if prot=='abs' and neut=='abs':
		j=0
		lines = open('sct'+str(k[j])+'.dat').readlines()			
		table = np.array([row.split() for row in lines])
		Eg = np.float64(np.delete(table,[1,2,3,4,5,6],1)).ravel()
		xsg = np.float64(np.delete(table,[0,2,3,4,5,6],1)).ravel()
		T0 = np.float64(np.delete(table,[0,1,3,4,5,6],1)).ravel()
		T1 = np.float64(np.delete(table,[0,1,2,4,5,6],1)).ravel()
		quad = np.float64(np.delete(table,[0,1,2,3,5,6],1)).ravel()
		ober = np.float64(np.delete(table,[0,1,2,3,4,6],1)).ravel()
		deut = np.float64(np.delete(table,[0,1,2,3,4,5],1)).ravel()
		j=j+1
		while j<len(k):
			lines = open('sct'+str(k[j])+'.dat').readlines()			
			table = np.array([row.split() for row in lines])
			xsg = xsg + np.float64(np.delete(table,[0,2,3,4,5,6],1)).ravel()
			T0 = T0 + np.float64(np.delete(table,[0,1,3,4,5,6],1)).ravel()
			T1 = T1 + np.float64(np.delete(table,[0,1,2,4,5,6],1)).ravel()
			quad = quad + np.float64(np.delete(table,[0,1,2,3,5,6],1)).ravel()
			ober = ober + np.float64(np.delete(table,[0,1,2,3,4,6],1)).ravel()
			deut = deut + np.float64(np.delete(table,[0,1,2,3,4,5],1)).ravel()
			j=j+1
		xsmax=max(xsg)
		for j in range(len(xsg)):
			if xsg[j]==max(xsg):
				Emax=Eg[j]
	else:
		for i in range(len(k)):
			if p[i]==prot and n[i]==neut:
				break
		xsmax=sec[i]
		Emax=energ[i]
		lines = open('sct'+str(k[i])+'.dat').readlines()			
		table = np.array([row.split() for row in lines])
		Eg = np.float64(np.delete(table,[1,2,3,4,5,6],1)).ravel()
		xsg = np.float64(np.delete(table,[0,2,3,4,5,6],1)).ravel()
		T0 = np.float64(np.delete(table,[0,1,3,4,5,6],1)).ravel()
		T1 = np.float64(np.delete(table,[0,1,2,4,5,6],1)).ravel()
		quad = np.float64(np.delete(table,[0,1,2,3,5,6],1)).ravel()
		ober = np.float64(np.delete(table,[0,1,2,3,4,6],1)).ravel()
		deut = np.float64(np.delete(table,[0,1,2,3,4,5],1)).ravel()
	result = namedtuple('arrays', ['xsmax','Emax','Eg','xsg','T0','T1','quad','ober','deut'])
	return result(xsmax,Emax,Eg,xsg,T0,T1,quad,ober,deut)

#take spectra from CMPR(Orlin), added Feb 2020
#inputs:
#param - emitted particle. param = 'n' or param = 'p'
#E - gamma energy
#outputs:
#Enuc - emitted particle energy
#spec - full nucleon spectrum
#gdr - GDR contribution to spectrum
#ober - obertone contribution to spectrum
#quad - quadrupole contribution to spectrum
#deut - quasideutrone contribution to spectrum
#example:
#	E=sinp.spor('p',25).Enuc
#	xs=sinp.spor('p',25).spec
def spor(param,E):
	lines = open(param+'sp.dat').readlines()
	table = np.array([row.split() for row in lines])
	Eg = np.float64(np.delete(table,[1,2,3,4,5,6],1)).ravel()
	Enucf = np.float64(np.delete(table,[0,2,3,4,5,6],1)).ravel()
	specf = np.float64(np.delete(table,[0,1,3,4,5,6],1)).ravel()
	gdrf = np.float64(np.delete(table,[0,1,2,4,5,6],1)).ravel()
	oberf = np.float64(np.delete(table,[0,1,2,3,5,6],1)).ravel()
	quadf = np.float64(np.delete(table,[0,1,2,3,4,6],1)).ravel()
	deutf = np.float64(np.delete(table,[0,1,2,3,4,5],1)).ravel()
	j=0
	for i in range(len(Eg)):
		if Eg[i]==E:
			j=j+1
	Enuc=np.zeros(j)
	spec=np.zeros(j)
	gdr=np.zeros(j)
	ober=np.zeros(j)
	quad=np.zeros(j)
	deut=np.zeros(j)
	j=0
	for i in range(len(Eg)):
		if Eg[i]==E:
			Enuc[j]=Enucf[i]
			spec[j]=specf[i]
			gdr[j]=specf[i]
			ober[j]=oberf[i]
			quad[j]=quadf[i]
			deut[j]=deutf[i]
			j=j+1
	result = namedtuple('arrays', ['Enuc','spec','gdr','ober','quad','deut'])
	return result(Enuc,spec,gdr,ober,quad,deut)

#make one plot, added Nov 2018
#Can be used for plotting cross sections
#Beta version
#if there is no error data, use err=0
def plot(param,lab,name,x,y,err,xlab='$E$, МэВ',ylab='$\sigma$, мб'):
	fig = plt.figure()
	fig.set_rasterized(False)
	if param == 'exp':
		plt.scatter(x,y,c='k',s=20,marker='o',label=lab)
		plt.errorbar(x, y, yerr=err, ls='none',c='k',elinewidth=1) 
	if param == 'th':
		plt.plot(x,y,c='k',linewidth=2,label=lab)
	plt.tight_layout()
	plt.legend()
	plt.xlabel(xlab,fontsize=14)
	plt.ylabel(ylab,fontsize=14)
	plt.xlim(np.float64(x[0]),)
	plt.ylim(0,)
	default_size = fig.get_size_inches() 
	fig.set_size_inches(default_size[0]*1.5, default_size[1]*1.5,forward=True)
	plt.savefig(name+'.eps',fmt='eps',rasterized=False)
	return fig
