
import sympy as sp
import numpy as np
import pandas as pd
import pint
u = pint.UnitRegistry()
Q= u.Quantity
import scipy.constants as cnst
import scipy

from matplotlib import pyplot as plt 
from sympy.interactive import printing
sp.init_printing(use_latex=True) # doctest: +SKIP

from IPython.display import HTML
from iapws import IAPWS97

def Presion_Vapor(Sustancia,T):
    """
    Sustancia en forma vector de libro Perry's
    Temperatura en KELVIN
    
    """
    C1=Sustancia[4][0]
    C2=float(Sustancia[5][0].replace(u",", "").replace(u"\u2212", "-"))
    C3=float(Sustancia[6][0].replace(u"\u2212", "-"))
    if(type(Sustancia[7][0])==str):
        C4= float(Sustancia[7][0].replace(u"\u2212", "-"))
    else:
        C4=Sustancia[7][0]
    C5=Sustancia[8][0]

    Pvap=Q(np.exp(C1 + C2/T + C3*np.log(T) + C4*T**C5),"Pa")
    return(Pvap)


def Funcion_Caracteristica_Burbuja (T,Sustancias,Composiciones,PresionSistema):
    
    zA,zB,zC,zD,zE,zF= Composiciones
    S1,S2,S3,S4,S5=Sustancias
    P1=zA*Presion_Vapor(S1,T)
    P2=zB*Presion_Vapor(S2,T)
    P3=zC*Presion_Vapor(S3,T)
    P4=zD*Presion_Vapor(S4,T)
    P5=zE*Presion_Vapor(S5,T)
    P6=zF*Presion_Vapor(S6,T)
    return((P1+P2+P3+P4+P5+P6-PresionSistema).to("Pa").magnitude)


def Funcion_Caracteristica_Rocio (T,Sustancias,Composiciones,PresionSistema):
    zA,zB,zC,zD,zE,zF= Composiciones
    S1,S2,S3,S4,S5,S6=Sustancias
    P1=zA/Presion_Vapor(S1,T)
    P2=zB/Presion_Vapor(S2,T)
    P3=zC/Presion_Vapor(S3,T)
    P4=zD/Presion_Vapor(S4,T)
    P5=zE/Presion_Vapor(S5,T)
    P6=zF/Presion_Vapor(S6,T)
    return(((PresionSistema*(P1+P2+P3+P4+P5+P6)-1)).magnitude)

    