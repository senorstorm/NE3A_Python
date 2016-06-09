import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

ka, kb = 0.0318, 0.05
aa, ab = 0.694, 0.5
Na = 6.02e23
vcell = 1e-3
Meg = 2.016 
Mo = 143.186
p = 0.995
x = np.linspace(0,50,50)

## Molecular weight of styrene
Mpst = np.exp(-0.3*x + 20)

## hydrodynamic volume
Vh = (ka*Mpst**(aa +1))/(2.5*Na)

## Molecular weight of polyurethane
Mpur = np.exp((1/(1+ab))*np.log((2.5*Vh*Na)/kb))

## Degree of polymerization
Xpur = (Mpur - Meg)/Mo

## 2.3b Probability distribution of Xn w.r.t Vel (Molar fraction)
sum_xpr = sum((p**(1+Xpur)*(1-p)**2))
Pr = (((p**(Xpur-1))*(1-p)**2)/sum_xpr)

##2.3a Normalizing DRI
sum_xpr_norm = sum((p**(1+Xpur)*(1-p)**2)*Mpur)
Pr_norm = (((p**(Xpur-1))*(1-p)**2)/sum_xpr_norm)*Mpur

## Plotting Normalized DRI and Molar fraction
fig = plt.figure()
ax = fig.add_subplot(111)
lns1 = ax.plot(x,Pr_norm,'-r', label='2.3a: Normalized DRI')
ax2 = ax.twinx()
lns2 = ax2.plot(x,Pr,'-', label = '2.3b: Molar Fraction')
ax.grid()
lns = lns1+lns2
labs = [i.get_label() for i in lns]
ax.legend(lns,labs,loc=0)
ax.set_xlabel("Elution volume (mL)")
ax.set_ylabel("2.3a")
ax2.set_ylabel("2.3b")
ax.set_ylim(0,14e-2)
ax2.set_ylim(0,9e-2)
plt.show()

#2.3d Vmax of Normalized DRI
print(max(Pr_norm))