import numpy as np
from scipy.interpolate import CubicSpline

def getsub(f, x, xref):
    f_cs = CubicSpline(x, f)
    return f_cs(xref)

def Qstream(ne, Te):
    me = 9.1094e-28 # [g]
    eV2K = 1.1604e4 # K = eV2K * eV
    erg2J = 1e-7
    kB = 1.3807e-16 * eV2K # [erg/eV]
    # Local thermal energy density
    eTh = ne * kB * Te
    # Thermal velocity
    vTh = (kB * Te / me)**0.5
    # Free-streaming heat flux [cm/s*erg/cm3]
    Qfs = vTh * eTh
    return erg2J * Qfs

# Limited heat flux (logistic weighting of Qloc and Qfs) 
def Qeff(X, flim):
    #fit function for Qloc profile
    ne, Z, Te, gradTe = X
    kQSH = 6.1e+02 # scaling constant corresponding to the SHICK local heat flux
    Qloc = -(kQSH/Z)*((Z+0.24)/(Z+4.2))*Te**2.5*gradTe
    Qfs = Qstream(ne, Te)
    Qeff = flim * Qfs * (1.0 - np.exp(-Qloc/(flim*Qfs)))
    return Qeff

def QSHlimited(x, ne, Zbar, Te, flim):
    #calculating Te gradient
    N = int(len(x))
    xref = np.linspace(min(x), max(x), N)
    Te = getsub(Te, x, xref)
    ne = getsub(ne, x, xref)
    Zbar = getsub(Zbar, x, xref)
    gradTe = np.gradient(Te, xref)
    return xref, Qeff((ne, Zbar, Te, gradTe), flim)

def Tselfsimilar_beta(x0, T0, x1, T1, x, beta):
    # TODO find a way of passing x0, T0, x1, T1 along with x
    def Tselfsimilar(X, beta):
        x = X
        b = ( (1.0 - (T1 / T0)**beta) / 
            (x1**2.0 - x0**2.0 * (T1 / T0)**beta) )
        c = T0 / (1.0 - b * x0**2.0)**(1.0 / beta)
        #print(f'T0 {T0} T0(b, c) {c * (1.0 - b * x0**2.0)**(1.0/beta)}')
        #print(f'T1 {T1} T1(b, c) {c * (1.0 - b * x1**2.0)**(1.0/beta)}')
        return c * (1.0 - b * x**2.0)**(1.0 / abs(beta))
    # par, cov = curve_fit(Tselfsimilar, xpoints, Tpoints, maxfev = 1000)
    return Tselfsimilar(x, beta)
