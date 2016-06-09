import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

t_list = [600,2400,10800,28800,259200]
x = np.linspace(1,51,51)
x_ = np.exp(0.135*x)

## 3.3 - Find conversion factors at discrete times
def conversion_factor(t_):
     k = 5e-3
     Mo = 2
     return 1-((1+t_*2*k*Mo**2)**(-1/2))    

p = [conversion_factor(i) for i in t_list]

##3.3 - Find probability distribution of p
N_x = []
N_xi = []
Wt_ = []
Wt_n = []

for i in range(0,len(p)):
    N_x.insert(i,(p[i]**(x_-1))*((1-p[i])**2))
    N_xi.insert(i,(N_x[i]/sum(N_x[i])))
    Wt_.insert(i,((N_xi[i])*(246*x_+18.02)))
    Wt_n.insert(i,(Wt_[i]/sum(Wt_[i])))

##3.3 - Plotting for Molar Fraction v Deg of Poly 
fig = plt.figure()    
ax = fig.add_subplot(111)

lns1 = ax.plot(x_,N_xi[0],label='t=10min')
lns2 = ax.plot(x_,N_xi[1],label='t=40min')
lns3 = ax.plot(x_,N_xi[2],label='t=3hours')
lns4 = ax.plot(x_,N_xi[3],label='t=8hours')
lns5 = ax.plot(x_,N_xi[4],label='t=3days')

lns = lns1+lns2+lns3+lns4+lns5
labs = [i.get_label() for i in lns]
ax.set_title('Molar Fraction vs Degree of Polymerization')
ax.set_xlabel('Degree of polymerization (x)')
ax.set_ylabel('Molar Fraction (Nx)')
ax.set_xscale('log')
ax.legend(lns,labs,loc=0)
plt.show()

##3.4 
fig = plt.figure()    
ax = fig.add_subplot(111)

lns1 = ax.plot(x_,Wt_n[0],label='t=10min')
lns2 = ax.plot(x_,Wt_n[1],label='t=40min')
lns3 = ax.plot(x_,Wt_n[2],label='t=3hours')
lns4 = ax.plot(x_,Wt_n[3],label='t=8hours')
lns5 = ax.plot(x_,Wt_n[4],label='t=3days')

lns = lns1+lns2+lns3+lns4+lns5
labs = [i.get_label() for i in lns]
ax.set_title('Weight Distribution of Degree of Polymerization')
ax.set_xlabel('Degree of polymerization (x)')
ax.set_ylabel('Weight Fraction (wt%)')
ax.set_xscale('log')
ax.legend(lns,labs,loc=0)
plt.show()
