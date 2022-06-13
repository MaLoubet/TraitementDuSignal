from scipy.integrate import quad
from numpy import cos
from numpy import sin 
from math import pi
from numpy import linspace
from math import floor
from scipy.fftpack import fft
from cmath import exp
from matplotlib.widgets import Slider, RadioButtons, Button
import matplotlib.pyplot as plt


def filtrage(e, T, H, x0, Q, g, d, p) :
    C = calculCoeffs(e, T, 70)
    S = [H(0, x0, Q) * C[0][0] / 2] * p
    s = [0] * p
    i = 0
    for t in linspace(g, d, p) :
        for n in range(70) :
            S[i] += complex(C[n][0],- C[n][1]) * H(n * 2 * pi / T, x0, Q) * exp(complex(0,n * 2 * pi * t / T))
        s[i] = S[i].real
        i += 1
    return s


def calculCoeffs(e, T, N) :
    coeffs = []
    a = 0
    b = 0
    for n in range(N) :
        a = 2/T * quad(lambda t : e(t)*cos(n*2*pi*t/T), 0, T)[0]
        b = 2/T * quad(lambda t : e(t)*sin(n*2*pi*t/T), 0, T)[0]
        coeffs.append([a,b])
    return coeffs

def H1(x, x0, Q) :
    return 1/ complex(1, x/x0)
    
def H2(x, x0, Q) :
    return complex(0, x/x0) / complex(1, x/x0)
    
def H3(x, x0, Q) :
    return 1 / complex(1 - (x/x0) ** 2, Q*x/x0)
    
def H4(x, x0, Q) :
    return complex(0, Q*x/x0) / complex(1 - (x/x0) ** 2, Q*x/x0)
    
def H5(x, x0, Q) :
    return complex(- (x/x0)**2,0) / complex(1 - (x/x0) ** 2, Q*x/x0)
    
def H6(x, x0, Q) :
    return complex(1 - (x/x0) ** 2, 0) / complex(1 - (x/x0) ** 2, Q*x/x0)
    
def triang(t) :
    return (1 - (floor(t) % 2)) * 2 * (t - floor(t)) - (floor(t) % 2) * 2 * (t -1 - floor(t)) - 1 
    
def creneaux(t) :
    return (-1) ** (floor(t) % 2)
   
def sinp(t) :
    return sin(pi * t)

def partiel(e, T, g, d, n, p) :
    E = [0] * n
    i = 0
    coeffs = calculCoeffs(e, T, p)
    for t in linspace(g, d, n) :
        E[i] += coeffs[0][0] / 2
        for v in range(p) :
            E[i] += (coeffs[v][0] * cos(v * 2 * pi * t / T) + coeffs[v][1] * sin(v * 2 * pi * t / T))
        i+=1
    return E

axR = plt.axes([0.01, 0.94, 0.1, 0.05])
bR = Button(axR, 'Reset')

axN = plt.axes([0.10, 0.03, 0.31, 0.03])
plt.title('Nombre de coefficients du signal partiel', size = 9)
sN = Slider(axN, 'N', 1, 200, valinit = 70)

axx0 = plt.axes([0.57, 0.05, 0.35, 0.02])
sx0 = Slider(axx0, 'x0', 0.1, 10, valinit = 5)

axQ = plt.axes([0.57, 0.02, 0.35, 0.02])
sQ = Slider(axQ, 'Q', 0.01, 10, valinit = 5)
    
axsig = plt.axes([0.13, 0.7, 0.35, 0.25])
plt.title("Signal d'entrée")
rsig = RadioButtons(axsig, ('Triangulaire', 'Créneaux', 'Sinusoïdale'))

axfiltre = plt.axes([0.53, 0.7, 0.45, 0.25])
plt.title("Filtrage")
rfil = RadioButtons(axfiltre, ('Premier ordre : passe-bas', 'Premier ordre : passe-haut', 'Second ordre : passe-bas','Second ordre : passe-bande','Second ordre : passe-haut','Second ordre : coupe-bande'))

def updateD(val) :
    x0 = sx0.val
    Q = sQ.val
    if rfil.value_selected == "Premier ordre : passe-bas" :
        if rsig.value_selected == "Triangulaire" :
            S = filtrage(triang, 2, H1, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Créneaux" :
            S = filtrage(creneaux, 2, H1, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Sinusoïdale" :
            S = filtrage(sinp, 2, H1, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
    elif rfil.value_selected == "Premier ordre : passe-haut" :
        if rsig.value_selected == "Triangulaire" :
            S = filtrage(triang, 2, H2, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Créneaux" :
            S = filtrage(creneaux, 2, H2, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Sinusoïdale" :
            S = filtrage(sinp, 2, H2, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
    elif rfil.value_selected == "Second ordre : passe-bas" :
        if rsig.value_selected == "Triangulaire" :
            S = filtrage(triang, 2, H3, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Créneaux" :
            S = filtrage(creneaux, 2, H3, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Sinusoïdale" :
            S = filtrage(sinp, 2, H3, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
    elif rfil.value_selected == "Second ordre : passe-bande" :
        if rsig.value_selected == "Triangulaire" :
            S = filtrage(triang, 2, H4, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Créneaux" :
            S = filtrage(creneaux, 2, H4, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Sinusoïdale" :
            S = filtrage(sinp, 2, H4, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
    elif rfil.value_selected == "Second ordre : passe-haut" :
        if rsig.value_selected == "Triangulaire" :
            S = filtrage(triang, 2, H5, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Créneaux" :
            S = filtrage(creneaux, 2, H5, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Sinusoïdale" :
            S = filtrage(sinp, 2, H5, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
    elif rfil.value_selected == "Second ordre : coupe-bande" :
        if rsig.value_selected == "Triangulaire" :
            S = filtrage(triang, 2, H6, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Créneaux" :
            S = filtrage(creneaux, 2, H6, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Sinusoïdale" :
            S = filtrage(sinp, 2, H6, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
    plt.draw()
    return None
 
def updateG(val) :
    N = floor(sN.val)
    if rsig.value_selected == "Triangulaire" :
        G2.set_ydata(partiel(triang, 2, -10, 10, 1000, N))
    elif rsig.value_selected == "Créneaux" :
        G2.set_ydata(partiel(creneaux, 2, -10, 10, 1000, N))
    elif rsig.value_selected == "Sinusoïdale" :
        G2.set_ydata(partiel(sinp, 2, -10, 10, 1000, N))
    plt.draw() 
    
def reset(event) : 
    sQ.reset()
    sN.reset()
    sx0.reset()
    rsig.set_active(0)
    rfil.set_active(0)

def filtre(value_selected) :
    sQ.reset()
    sx0.reset()
    sN.reset()
    if value_selected == "Premier ordre : passe-bas" :
        if rsig.value_selected == "Triangulaire" :
            S = filtrage(triang, 2, H1, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Créneaux" :
            S = filtrage(creneaux, 2, H1, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Sinusoïdale" :
            S = filtrage(sinp, 2, H1, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
    elif value_selected == "Premier ordre : passe-haut" :
        if rsig.value_selected == "Triangulaire" :
            S = filtrage(triang, 2, H2, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Créneaux" :
            S = filtrage(creneaux, 2, H2, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Sinusoïdale" :
            S = filtrage(sinp, 2, H2, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
    elif value_selected == "Second ordre : passe-bas" :
        if rsig.value_selected == "Triangulaire" :
            S = filtrage(triang, 2, H3, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Créneaux" :
            S = filtrage(creneaux, 2, H3, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Sinusoïdale" :
            S = filtrage(sinp, 2, H3, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
    elif value_selected == "Second ordre : passe-bande" :
        if rsig.value_selected == "Triangulaire" :
            S = filtrage(triang, 2, H4, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Créneaux" :
            S = filtrage(creneaux, 2, H4, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Sinusoïdale" :
            S = filtrage(sinp, 2, H4, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
    elif value_selected == "Second ordre : passe-haut" :
        if rsig.value_selected == "Triangulaire" :
            S = filtrage(triang, 2, H5, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Créneaux" :
            S = filtrage(creneaux, 2, H5, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Sinusoïdale" :
            S = filtrage(sinp, 2, H5, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
    elif value_selected == "Second ordre : coupe-bande" :
        if rsig.value_selected == "Triangulaire" :
            S = filtrage(triang, 2, H6, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Créneaux" :
            S = filtrage(creneaux, 2, H6, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rsig.value_selected == "Sinusoïdale" :
            S = filtrage(sinp, 2, H6, x0, Q, -10, 10, 1000)
            Sf = fft(S)
        
            D1.set_ydata(S)
            D2.set_ydata(Sf)
    plt.draw()
    
    
def entree(value_selected) :
    sQ.reset()
    sx0.reset()
    sN.reset()
    if value_selected == "Triangulaire" :
        if rfil.value_selected == "Premier ordre : passe-bas" :
            E = [triang(n) for n in t]
            Ep = partiel(triang, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(triang, 2, H1, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rfil.value_selected == "Premier ordre : passe-haut" :
            E = [triang(n) for n in t]
            Ep = partiel(triang, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(triang, 2, H2, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rfil.value_selected == "Second ordre : passe-bas" :
            E = [triang(n) for n in t]
            Ep = partiel(triang, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(triang, 2, H3, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rfil.value_selected == "Second ordre : passe-bande" :
            E = [triang(n) for n in t]
            Ep = partiel(triang, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(triang, 2, H4, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rfil.value_selected == "Second ordre : passe-haut" :
            E = [triang(n) for n in t]
            Ep = partiel(triang, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(triang, 2, H5, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rfil.value_selected == "Second ordre : coupe-bande" :
            E = [triang(n) for n in t]
            Ep = partiel(triang, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(triang, 2, H6, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
    elif value_selected == "Créneaux" :
        if rfil.value_selected == "Premier ordre : passe-bas" :
            E = [creneaux(n) for n in t]
            Ep = partiel(creneaux, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(creneaux, 2, H1, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rfil.value_selected == "Premier ordre : passe-haut" :
            E = [creneaux(n) for n in t]
            Ep = partiel(creneaux, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(creneaux, 2, H2, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rfil.value_selected == "Second ordre : passe-bas" :
            E = [creneaux(n) for n in t]
            Ep = partiel(creneaux, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(creneaux, 2, H3, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rfil.value_selected == "Second ordre : passe-bande" :
            E = [creneaux(n) for n in t]
            Ep = partiel(creneaux, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(creneaux, 2, H4, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rfil.value_selected == "Second ordre : passe-haut" :
            E = [creneaux(n) for n in t]
            Ep = partiel(creneaux, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(creneaux, 2, H5, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rfil.value_selected == "Second ordre : coupe-bande" :
            E = [creneaux(n) for n in t]
            Ep = partiel(creneaux, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(creneaux, 2, H6, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
    elif value_selected == "Sinusoïdale" :
        if rfil.value_selected == "Premier ordre : passe-bas" :
            E = [sinp(n) for n in t]
            Ep = partiel(sinp, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(sinp, 2, H1, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rfil.value_selected == "Premier ordre : passe-haut" :
            E = [sinp(n) for n in t]
            Ep = partiel(sinp, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(sinp, 2, H2, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rfil.value_selected == "Second ordre : passe-bas" :
            E = [sinp(n) for n in t]
            Ep = partiel(sinp, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(sinp, 2, H3, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rfil.value_selected == "Second ordre : passe-bande" :
            E = [sinp(n) for n in t]
            Ep = partiel(sinp, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(sinp, 2, H4, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rfil.value_selected == "Second ordre : passe-haut" :
            E = [sinp(n) for n in t]
            Ep = partiel(sinp, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(sinp, 2, H5, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
        elif rfil.value_selected == "Second ordre : coupe-bande" :
            E = [sinp(n) for n in t]
            Ep = partiel(sinp, 2, -10, 10, 1000, N)
            Ef = fft(E)
                
            S = filtrage(sinp, 2, H6, x0, Q, -10, 10, 1000)
            Sf = fft(S)
                
            G1.set_ydata(E)
            G2.set_ydata(Ep)
            G3.set_ydata(Ef)
            D1.set_ydata(S)
            D2.set_ydata(Sf)
    plt.draw()
    
rfil.on_clicked(filtre)   
rsig.on_clicked(entree)
bR.on_clicked(reset)
sQ.on_changed(updateD)
sN.on_changed(updateG)
sx0.on_changed(updateD)

N = 70
Q = 5.0
x0 = 5.0
t = linspace(-10,10,1000)
E = [triang(n) for n in t]
Ep = partiel(triang, 2, -10, 10, 1000, N)
Ef = fft(E)

S = filtrage(triang, 2, H1, x0, Q, -10, 10, 1000)
Sf = fft(S)

plt.subplots_adjust(top = 0.65, bottom = 0.17, wspace = 0.5)

plt.subplot(3,2,1)
G1, = plt.plot(t, E)

plt.subplot(3,2,3)
plt.ylabel("signal partiel")
G2, = plt.plot(t, Ep)

plt.subplot(3,2,5)
plt.ylabel("spectre de fourier")
G3, = plt.plot(Ef)

plt.subplot(2,2,2)
plt.ylabel("signal de sortie")
D1, = plt.plot(t, S)

plt.subplot(2,2,4)
plt.ylabel("Spectre de Fourier")
D2, = plt.plot(t, Sf)
plt.show()