import numpy as np
import matplotlib.pyplot as plt
'''
N=10000
M=10
tlim = (0.,24.)
Ti_LS = 35.
Ti_tt = 37.
#k_LS = (500*1.9/(70*4186)*60)
#k_tt = (100*4.5/(290*205)*60)
Vm = 0.187249
day_H=24.
Tmin_H = -60.
r_max = 0.310037
r_max_luke = (1.9/(4*np.pi))**(1/2)

volume_max = 0.187249
volume_max_luke = ((3*r_max_luke)/(4*np.pi))**(1/2)
'''

def cooling_constant_computation(radius):## KLC: you need to pass in a radius variable
    r_max_tt = 0.310037
    r_max_luke = (1.9/(4*np.pi))**(1/2)
    volume_max_luke = ((3*r_max_luke)/(4*np.pi))**(1/2)
    volume_max_tt= 0.187249
    
    ## KLC: But you also nee an appropriately sized container
    cooling_constant = np.empty(len(radius))
    for rr in range(len(radius)):
       ## R = radius[rr] # KLC: doing this way more robust to diff in
                       # Python v2 vs v3 *and* enables rr to index
                       # cooling_constant (for R in radius would have
                       # R not applicable)
        height = r_max_luke + radius[rr]
        if radius[rr] <= r_max_luke: ##Below are computations for luke
            Surface_Area = 2*np.pi*radius[rr]*(radius[rr]+height)
            Volume = (4/3)*np.pi*(radius[rr]**3) ##Volume of a sphere for Luke
            Mass = 70*(Volume/volume_max_luke)
            cooling_constant[rr] = 60*(500*Surface_Area)/(Mass*4186)
        elif radius[rr] > r_max_luke and radius[rr] < r_max_tt:  ##Below are computations for tauntaun
            Surface_Area = 2*np.pi*radius[rr]*(radius[rr]+height)
            Volume = np.pi*(radius[rr]**2)*height ## Volume of a cylinder for tauntaun
            Mass = 290*(Volume/volume_max_tt)
            cooling_constant[rr] = 60*(100*Surface_Area)/(Mass*205)
    return cooling_constant # KLC: return the array (not within for-loop)


def dTdt(T_now,radius,Tmin_H = -60):
    nshells = len(T_now)
    dTdt_now = np.empty(nshells)
#    dTdt_now[0] = -cooling_constant_computation[0]*(T_now[0] - T_now[1])

    ## KLC: functions are accessed with parentheses, not square
    ## brackets but it seems you actually want an array of cooling
    ## constants so since the radius variable is defined and
    ## instantiated outside any function, it is available to all
    ## functions. So we'll use the function
    ## cooling_constant_computation() on the "global-ish" variable
    ## radius to define an array of cooling constants. Then use it in
    ## the calculations.
    cooling_const = cooling_constant_computation(radius) ## KLC: but really, pass in or make global variable e.e., k_r
    dTdt_now[0] = -cooling_const[0]*(T_now[0] - T_now[1])
    for ss in range(1,nshells-1): # up to but not including outer shell
#        print('ss = {0}'.format(ss))
        dTdt_now[ss] = -cooling_const[ss]*(T_now[ss] - T_now[ss+1])
    ## Here would be another for-loop or if-statement to handle next object
    dTdt_now[-1] = -cooling_const[-1]*(T_now[-1] - Tmin_H)# - T_now[4]) error in T_now[4]: 
    return dTdt_now            
    
#T_Hoth = 0.5*Tmin_H*(1+np.sin(2*np.pi*time/day_H))

'''
T_all = np.empty((N,M),float)

## Modeling initial temperature gradient (linear from core to surface)
T_inner = 37. # C; human core temp
T_outer = 32. # C; human freezing
slope = (T_outer - T_inner)/radius[-1] # C/m where origin at radius = 0 
T_init = slope*radius + T_inner # C

#T_all[0,:] = np.array([37,35,33,31],float) # arbitrary instantiation
T_all[0,:] = T_init 
#[0],radius[1],radius[2],radius[3]],float)
 
radius = np.linspace(0,0.310037,M+1) # +1 to make M shells
radius = radius[1:] # want radius of outer parts of shells
#above line eliminates first place 0, to avoid a division by 0
'''    
def answer(Ti_LS = 35., Ti_tt = 37., day_H = 24, tlim = (0.,24.), Tmin_H = -60,
           N=10000, M=10):
    radius = np.linspace(0,0.310037,M+1) # +1 to make M shells
    radius = radius[1:]
    
    time = np.linspace(tlim[0],tlim[1],N)
    dt = time[1] - time[0]
    
    T_all = np.empty((N,M),float)
    T_inner = 37. # C; human core temp
    T_outer = 32. # C; human freezing
    slope = (T_outer - T_inner)/radius[-1] # C/m where origin at radius = 0 
    T_init = slope*radius + T_inner # C

#T_all[0,:] = np.array([37,35,33,31],float) # arbitrary instantiation
    T_all[0,:] = T_init 

    
    ## KLC: make all the floating "preamble" definitions part of this
    ## function
    for ii in range(1,N):
        k1 = dt * dTdt(T_all[ii-1,:],radius)
        
        k2 = dt * dTdt(T_all[ii-1,:]+0.5*k1,radius)
        
        k3 = dt * dTdt(T_all[ii-1,:]+0.5*k2,radius)
        
        k4 = dt * dTdt(T_all[ii-1,:]+k3,radius)
        
        T_all[ii,:] = T_all[ii-1,:] + (k1+ 2*(k2+k3) + k4)/6.
        
    for aa in range(M):
        plt.plot(time,T_all[:,aa],label='M = {0}'.format(aa))
#plt.semilogy()
    plt.legend()
    plt.show()
    print(T_all)
    return T_all


T_all = answer()


for ii in range(M):
    plt.plot(time,T_all[:,ii],label='M = {0}'.format(ii))
#plt.semilogy()
plt.legend()
plt.show()




## KLC: Josh's To Do
## - fix dTdt() calling cooling_constant() function to accept k_r array
## - fix answer() to have all the "preamble" stuff and allow user to
##   change parameters without re-compiling
## - test range of "sensible" parameters (e.g., nshell=1)
## - add tauntaun shell (represented by diff. cooling constant array)
## - add constant or variable Hoth temperature
## - enable user to change geometry (i.e., enable "spherical" mode)
