from tabulate import tabulate
from operator import itemgetter
import numpy as np
import matplotlib.pyplot as plt

def Dic2LatexTable(dic,tablefmt='latex'):
	'''Pre: dictionary for C and N contribution by each substrate (obtained from PercentCont)
	Post: returns latex table'''
	rv = []	
	for s in dic.keys():
		lst = [s,dic[s]['C'],dic[s]['N']]
		rv.append(lst)
	table = tabulate(sorted(rv, key=itemgetter(2), reverse=True), headers=["Substrate", "Carbon", "Nitrogen"],tablefmt=tablefmt)
	print (table)

def Plot(dic,xaxis='Percent Carbon Contribution',yaxis='Percent Nitrogen Contribution',filename='percent_cont.pdf'):
	'''Pre: dictionary for C and N contribution by each substrate (obtained from PercentCont)
	Post: returns Plot'''
	res = {}
	for sub in dic.keys():
		res[sub] = [dic[sub]['C'],dic[sub]['N']]

	# repackage data into array-like for matplotlib, pythonically
	xs,ys = zip(*res.values())
	labels = res.keys()   
	# display
	plt.xlabel(xaxis, fontsize=15)
	plt.ylabel(yaxis, fontsize=15)
	plt.xscale('symlog', base=2)
	plt.yscale('symlog', base=2)
	plt.scatter(xs, ys, marker = 'x')
	for label, x, y in zip(labels, xs, ys):
		plt.annotate(label, xy = (x, y))
	plt.savefig('../../Figs/'+filename)
	plt.show()

def PieChart(dic,byprod=True,remove=['NIACINE','CYS'],filename='percent_cont.svg'):
	fig, (ax1,ax2) = plt.subplots(1,2,figsize=(20,20))
	keys = []
	valsC = []
	valsN = []
	explode = []
	extC = 0.0
	extN = 0.0
	propdic = dic.copy()
	for s in dic.keys():
		if dic[s]['C'] < 0.0 or dic[s]['N'] < 0.0:
			extC += abs(dic[s]['C'])
			extN += abs(dic[s]['N'])
			propdic.pop(s)
	propdic['By-products'] = {'C':extC,'N':extN}
	for r in remove:
		propdic.pop(r)
	if byprod == False:
		propdic.pop('By-products')
	for s in propdic.keys():
		keys.append(s)
		valsC.append(propdic[s]['C'])
		valsN.append(propdic[s]['N'])
		if s == 'By-products':
			explode.append(0.1)
		else:
			explode.append(0.0)
	yc = np.array(valsC)
	yn = np.array(valsN)
	labels = keys
	ax1.pie(yc, labels =labels,explode = explode, autopct='%.2f')
	ax1.set_title('Percent C contribution')
	ax2.pie(yn, labels =labels,explode = explode, autopct='%.2f')
	ax2.set_title('Percent N contribution')
	plt.savefig('../../Figs/'+filename)
	plt.show() 
