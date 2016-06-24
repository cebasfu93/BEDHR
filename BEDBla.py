
# coding: utf-8

# # BEDbla: Boltzmann Enhanced Discrimination bla

# In[1]:

#Importa los paquetes necesarios
from __future__ import division
from pylab import *
#%pylab inline 
import csv, os, random, sys, shutil, time
from optparse import OptionParser
from scipy import stats
from scipy import optimize
from optparse import OptionParser
import numpy
numpy.set_printoptions(threshold=numpy.nan)

parser=OptionParser()
tamano=18


# In[2]:

parser.add_option("-i","--input", action="store", type="string", dest="InputFile", help="Path of the .csv input file. Special requirements: There must be a column with the numerical activity values and another one whose values will be used for the sorting (scoring function). The first row will be ignored.")
parser.add_option("-d", "--threshold", action="store", type="int", dest="Threshold", help="Entries with an activity lower or equal than this value will be considered as active molecules.")
parser.add_option("-a", "--cactivity", action="store", type="int", dest="NumericActivity", help="Number of the column (starting from 0) that contains the numeric activity of each entry.")
parser.add_option("-f", "--nbedef", action="store", type="int", dest="NBEDEF", help="Percentage of screened compounds responsible for 80% of BEDEF value. By default BEDHR(1,5,20%) will be displayed.")
parser.add_option("-s", "--scoring", action="store", type="int", dest="ScoreColumn", default=-1, help="Number of the column (starting from 0) that contains the values of the scoring function used for the sorting. Default: The last one")
parser.add_option("-e", "--curves", action="store_true", dest="Curves", default=False, help="Show enrichment curve and plot of enrichment factors as a function of N. Default: False")

(options, args) = parser.parse_args()

#args_test=["-i", "./BEDHR-Results/Dock10.csv", "-d", "1", "-a", "1", "-f", "20"]
#(options, args) = parser.parse_args(args_test)


# In[3]:

inputfile=options.InputFile

thresAct=options.Threshold #Por debajo (o igual) a este valor se considera el compuesto activo
coluMIC=options.NumericActivity
Nbedef=options.NBEDEF
curves=options.Curves

#########################################################################################################################

Datos=genfromtxt(inputfile, delimiter=',', skip_header=1)  


# In[4]:

def arreglar(array):
    filas=len(array[:,0])
    columnas=len(array[0,:])
    array_fix=zeros((filas , columnas+1))
    array_fix[:,1:columnas+1]=array
    
    #Asigna valor boolean a columna 0    
    for o in range(0,filas):
        if (array_fix[o, coluMIC+1]<=thresAct):
            array_fix[o, 0]=1
        if (array_fix[o, coluMIC+1]>thresAct):
            array_fix[o,0]=2
    return array_fix

#Arregla input
DatosArreg=arreglar(Datos)

if options.ScoreColumn == -1:
    scorecolu=len(DatosArreg[0,:])-1
else:
    scorecolu=options.ScoreColumn


# ## Enriquecimiento

# In[5]:

#Ordena el array de menor a mayor segun los valores de numcol
def Ordenador(array, numcol):
    
    Ordered=array[array[:,numcol].argsort()]
    
    return Ordered


# In[6]:

#Crea array con la cantidad de compuestos screened en cada slot
def CompScreened(arrayinp):
    
    compuestos_totales=len(arrayinp[:,0])
    compuestos_screened=1
    porcentaje_screened=zeros(compuestos_totales)
    
    for i in range(compuestos_totales):
        nuevo_procentaje=(compuestos_screened/compuestos_totales)*100
        compuestos_screened +=1
        porcentaje_screened[i]=nuevo_procentaje
        
    return porcentaje_screened

#Crea array con la cantidad de compuestos encontrados hasta cada slot
def ActScreened(array, numcol):
    
    arrayord=Ordenador(array, numcol)
    acttot=0
    actfound=0
    actscreen=zeros(len(array[:,0]))
    
    for h in range(0, len(arrayord[:,0])):
        if (arrayord[h,0] == 1): #Aca se busca el 1 en la columna 0
            acttot += 1
            
    for i in range(0, len(arrayord[:,0])):
        if (arrayord[i, 0] == 1):
            actfound += 1
        actpercent = actfound/acttot * 100
        actscreen[i]=actpercent
        
    return actscreen


# In[7]:

#Grafica curva de enriquecimiento ordenando segun valores de CODE
if curves == True:
    figure(figsize=(8,8)); hold=True
    p=arange(0,101,1)
    lin=p
    plot(p,lin, 'k-', lw=1, label=r'Random Ranking')
    plot(CompScreened(DatosArreg), ActScreened(DatosArreg, scorecolu), lw=1, label=r'Scoring Function')
    plot(CompScreened(DatosArreg), ActScreened(DatosArreg, 0), 'm-', lw=1, label=r'Ideal Behaviour')
    title(r'Enrichment Curve', fontsize=tamano)
    xlabel(r'Screened Compounds (%)', fontsize=tamano)
    ylabel(r'Active Compounds Found (%)', fontsize=tamano)
    legend(fontsize=tamano, fancybox=True, loc=2).get_frame().set_alpha(0.5)
    grid()
    savefig("EnrichmentCurve.png")


# ## BEDHR

# In[8]:

def EnrichFacN(screen_array, act_array):
    pend_array=[]
    for i in range (0, len(screen_array)):
        if act_array[i]!=100:
            pendiente=act_array[i]/screen_array[i]
            pend_array.append(pendiente)
        if act_array[i]==100:
            pend_array.append(pendiente)
    return pend_array

pendiente_CODE=EnrichFacN(CompScreened(DatosArreg), ActScreened(DatosArreg, scorecolu))
pendiente_ideal=EnrichFacN(CompScreened(DatosArreg), ActScreened(DatosArreg, 0))


# In[9]:

#Curva de enriquecimiento para todos los N
if curves==True:
    figure(figsize=(8,8)); hold=True
    plot(CompScreened(DatosArreg), pendiente_CODE, label="Scoring Function", lw=2)
    plot(CompScreened(DatosArreg), pendiente_ideal, label="Ideal Behaviour", lw=2)
    title(r'Enrichment Factors at All N Values', fontsize=tamano)
    xlabel("Screened Compounds (N%)", fontsize=tamano)
    ylabel("Enrichment Factor at the N%", fontsize=tamano)
    legend(fontsize=tamano, fancybox=True, loc=1).get_frame().set_alpha(0.5)
    grid()
    savefig('EnrichmentFactors.png')


# In[10]:

def pendiente(screen_array, act_array):
    pend_array=[]
    for i in range (0, len(screen_array)):
        if act_array[i]!=100:
            pendiente=(act_array[i+1]-act_array[i])/(screen_array[i+1]-screen_array[i])
            pend_array.append(pendiente)
        if act_array[i]==100:
            pend_array.append(pendiente)
    return pend_array

pend_inst_code=pendiente(CompScreened(DatosArreg), ActScreened(DatosArreg, scorecolu))
pend_inst_ideal=pendiente(CompScreened(DatosArreg), ActScreened(DatosArreg, 0))


# In[11]:

#Curva de enriquecimientos no acumulativos
if curves==True:
    figure(figsize=(8,8)); hold=True
    plot(CompScreened(DatosArreg), pend_inst_code, label="Scoring Function")
    plot(CompScreened(DatosArreg), pend_inst_ideal, label="Ideal Behaviour", lw=2)
    title(r'Non-Acumulative Enrichment Factors', fontsize=tamano)
    xlabel("Screened Compounds (N%)", fontsize=tamano)
    ylabel("Instant Slope of Enrichment Curve at N%", fontsize=tamano)
    legend(fontsize=tamano, fancybox=True, loc=1).get_frame().set_alpha(0.5)
    grid()
    savefig("Non-acumulativeEnrichments.png")


# In[12]:

def alfacalc(x, Nimp):
    return exp(-x*Nimp/100)-0.8*exp(-x)-0.2


# In[13]:

def bedef(screened_comp, Nimp):
    
    bedef=0
    
    alfa = optimize.brentq(alfacalc, 0.0001, 1000000, args=(Nimp))
    
    weight=exp(-alfa*screened_comp/100)
    bedef_temp=dot(weight, pend_inst_code)
    norm=dot(weight, pend_inst_ideal)
    bedef=bedef_temp/norm
    return bedef

BEDEF20=bedef(CompScreened(DatosArreg), 20)
BEDEF5=bedef(CompScreened(DatosArreg), 5)
BEDEF1=bedef(CompScreened(DatosArreg), 1)
BEDEFN=bedef(CompScreened(DatosArreg), Nbedef)


# ## Log File

# In[14]:

arch = open('BEDbla.log', "a")
    
arch.write('Input file: ' + inputfile)

arch.write('\n BEDHR(20%): ' + str(BEDEF20))
arch.write('\n BEDHR(5%): ' + str(BEDEF5))
arch.write('\n BEDHR(1%): ' + str(BEDEF1))
arch.write('\n BEDHR(' + str(Nbedef) + '%): ' + str(BEDEFN))

arch.close()
    
msj='The BEDHR values were written to BEDHR.log'
print msj


# In[ ]:



