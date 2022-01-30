import glob, pandas as pd, matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cf

file=glob.glob('*.csv')

# data is in file[0]
df=pd.read_csv(file[0], delimiter=',').dropna()

# constants definition
q=1 
m0=9.1e-31
m=0.276*m0
h = 6.626e-34
h=h/1.6e-19 # eV.s
ħ=h/(2*np.pi) # eV.s


# As from julia it is seen that we can curtail at second term

def wkb(En, Eg, F, L): 
    ΔE1 = ( 0.75*ħ*q/(np.sqrt(2*m)) )**(2/3)
    ΔE2 = ( 1.875*ħ*(q**2)/(np.sqrt(2*m)) )**(2/5)
    
    return ( abs(Eg - En)/(ΔE1*(F**(2/3) ) ) )**(3/2)   +  ( abs(Eg - En)/(ΔE2*( ( L*(F**2) )**(2/5) ) ) )**(5/2) 


# curtail the input vectors at desired intervals

En, Yh=df.iloc[:,0], df.iloc[:,1]

#initial condition
p0=[5.0,1.5e+3,150e-9] # [Eg, F, L]

# fittiong

popt, pcov = cf(wkb, En, Yh, p0=p0)

plt.plot(En,wkb(En,*popt), label='fit.')
plt.scatter(En, Yh, label='exp.')
plt.show()


