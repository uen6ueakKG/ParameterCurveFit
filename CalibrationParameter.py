# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 12:09:42 2021

@author: lowkg
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

Te = ( 38.58, 95.09,211.81,420.01,855.64)
Fc= ( 3.42,14.70,20.28,25.21,34.5)
'''
Te = ( 43.18, 100.69,217.11,423.74,855.64)
Fc= ( 13.98,19.23,38.54,45.58,52.05)
'''
Eq_Option = [0,1]
legend =[]

def Func_ScandinavianCode(Te, Fc_u, tau, a):
    # ScandinavianCode Equation to calculate Compressive Strength
    return Fc_u * (np.exp(-np.power((tau / Te), a)))

def Func_AmericanCode(Te, Fc_u, k, Te_0):
    # AmericanCode Equation to calculate Compressive Strength
    return Fc_u*((k*(Te-Te_0))/(1+k*(Te-Te_0)))

# array with 3 elements in order of Fc_u,tau,a
def Fc_ParameterGraph(pars,option): 
    df = pd.DataFrame()
    Te_graph = []
    Fc_graph = []
    
    Fc_u = round(pars[0],4)
    par1 = round(pars[1],4)
    par2 = round(pars[2],4)
    
    iLim = max(Te)
    i = 0.5
    while i < iLim:
        Te_plot = i
        if option == 0:
            Fc_cal = Fc_u * (np.exp(-np.power((par1 / Te_plot), par2)))
            PlotLegend = 'ScandinavianEq'
        elif option == 1:
            Fc_cal = Fc_u*((par1*(Te_plot-par2))/(1+par1*(Te_plot-par2)))
            PlotLegend = 'AmericanEq'
        Te_graph.append(i)
        Fc_graph.append(Fc_cal)
        i+=0.5

    df['Te'] = Te_graph
    df['Fc'] = Fc_graph
    
    plt.scatter(Te, Fc)
    legend.append(PlotLegend)
    plt.plot(df['Te'], df['Fc'],linestyle='dashdot')
    plt.legend(legend)

def DeviationCalculation(pars,Te,Fc,option):
    Fc_u = round(pars[0],4)
    par1 =round(pars[1],4)
    par2 = round(pars[2],4)
    
    for i in range(len(Te)):
        Te_input = Te[i]
        if option == 0:
            Fc_cal = Fc_u * (np.exp(-np.power((par1 /Te_input), par2)))
        elif option == 1:
            Fc_cal = Fc_u*((par1*(Te_input-par2))/(1+par1*(Te_input-par2)))
            
        Dev = round((Fc_cal - Fc[i])/Fc[i]*100,2)
        print("Fc: " + str(Fc[i]) + "\t Deviation: " + str(Dev))

def getStrengthParameter(Te,Fc):
    F_strength =[Func_ScandinavianCode,Func_AmericanCode]
    Fc_0 = 25
    for option in Eq_Option:
        # Guestimate initial parameter
        if option == 0:
            par1 = 10 ; par2 = 0.6;
            param = ['Fc_u','tau','a']
        
        elif option == 1:
            par1 = 0 ; par2 = 0
            param = ['Fc_u','k','Te_0']
            
        pars, cov = curve_fit(f=F_strength[option], xdata=Te, ydata=Fc,
                          p0=[Fc_0,par1,par2])
        
        for i in range(len(param)):
            print(param[i] + ':' + str(round(pars[i],2)))
            
        Fc_ParameterGraph(pars,option)
        DeviationCalculation(pars, Te, Fc,option)
        
    return pars

plt.figure(figsize=(6, 4))
pars = getStrengthParameter(Te,Fc)



