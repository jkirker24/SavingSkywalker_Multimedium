import numpy as np
import matplotlib.pyplot as plt


N=10000
M=10
tlim = (0.,24.)
Ti_LS = 35.
Ti_tt = 37.
k_LS = (500*1.9/(70*4186)*60)
k_tt = (100*4.5/(290*205)*60)
Vm = 0.187249
day_H=24.
Tmin_H = -60.
r_max = 0.310037
r_max_luke = (1.9/(4*np.pi))**(1/2)

volume_max = 0.187249
volume_max_luke = ((3*r_max_luke)/(4*np.pi))**(1/2)

radius = np.linspace(0,0.310037,M+1) # +1 to make M shells
radius = radius[1:] # want radius of outer parts of shells
#above line eliminates first place 0, to avoid a division by 0
height = r_max_luke + radius

'''

##Below code works
surface_area = 2*np.pi*radius*(radius + height) # brute force would have been for-loop over M

volume = np.pi*(radius**2)*height

mass_luke = 70*(volume/volume_max_luke)

mass_tt = 290*(volume/Vm)


cooling_constant_luke = 60*(500*surface_area)/(mass_luke*4186)
cooling_constant_tt = 60*(100*surface_area)/(mass_tt*205)
##Above code works

'''


## Below is code in testing

def cooling_constant_computation(radius): ## KLC: you need to pass in a radius variable
    ## KLC: But you also nee an appropriately sized container
    nshells = len(radius)
    cooling_constant = np.empty(nshells)
    for rr in range(nshells):
        R = radius[rr] # KLC: doing this way more robust to diff in
                       # Python v2 vs v3 *and* enables rr to index
                       # cooling_constant (for R in radius would have
                       # R not applicable)
        height = 2*R
        if R <= r_max_luke: ##Below are computations for luke
            Surface_Area = 2*np.pi*R*(R+height)
            Volume = (4/3)*np.pi*(R**3) ##Volume of a sphere for Luke
            Mass = 70*(Volume/volume_max_luke)
            cooling_constant[rr] = 60*(500*Surface_Area)/(Mass*4186)
        elif R > r_max_luke and R < r_max:  ##Below are computations for tauntaun
            Surface_Area = 2*np.pi*R*(R+height)
            Volume = np.pi*(R**2)*height ## Volume of a cylinder for tauntaun
            Mass = 290*(Volume/volume_max)
            cooling_constant[rr] = 60*(100*Surface_Area)/(Mass*205)
    return cooling_constant # KLC: return the array (not within for-loop)

'''
Current Error states: TypeError: 'function' is not subscriptable
Online it seems this comes from an ambiguity objects

Intention is to have cooling_constant_computation be called by dTdt, but I 
do not know how the input variables work out.  I don't think cooling_constant_computation
needs an input variable, since it's constantly radius and height which are already defined.

Also I'm not certain how to address the [ss] loop index below for a function.
Previously I believe it was just an array or list, and easily indexed.  However,
now that it's a function I'm not sure if it becomes an array/list during the dTdt function.

'''

    
def dTdt(T_now,radius):
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
    
## Above is code in testing



    
        
##Below code works

#h*A/(m*c)*60
#k_Luke=calc_cooling_const(500,1.9,70,4186)
#k_tauntaun=calc_cooling_const(100,4.5,290,205)

time = np.linspace(tlim[0],tlim[1],N)
dt = time[1] - time[0]

'''##Below is old working code

def dTdt(T_now):
    nshells = len(T_now)
    dTdt_now = np.empty(nshells)
    dTdt_now[0] = -cooling_constant_luke[0]*(T_now[0] - T_now[1])
    for ss in range(1,nshells-1): # up to but not including outer shell
#        print('ss = {0}'.format(ss))
        dTdt_now[ss] = -cooling_constant_luke[ss]*(T_now[ss] - T_now[ss+1])
    ## Here would be another for-loop or if-statement to handle next object
    dTdt_now[-1] = -cooling_constant_luke[-1]*(T_now[-1] - Tmin_H)# - T_now[4]) error in T_now[4]: 
    return dTdt_now
'''



##Below Code to remain unaltered until functions for cooling_constant_computation and dTdt are finished 7/28/17


#T_Hoth = 0.5*Tmin_H*(1+np.sin(2*np.pi*time/day_H))

T_all = np.empty((N,M),float)

## Modeling initial temperature gradient (linear from core to surface)
T_inner = 37. # C; human core temp
T_outer = 32. # C; human freezing
slope = (T_outer - T_inner)/radius[-1] # C/m where origin at radius = 0 
T_init = slope*radius + T_inner # C

#T_all[0,:] = np.array([37,35,33,31],float) # arbitrary instantiation
T_all[0,:] = T_init 
#[0],radius[1],radius[2],radius[3]],float)
 
def answer():
    ## KLC: make all the floating "preamble" definitions part of this
    ## function
    for ii in range(1,N):
        k1 = dt * dTdt(T_all[ii-1,:],radius)
        
        k2 = dt * dTdt(T_all[ii-1,:]+0.5*k1,radius)
        
        k3 = dt * dTdt(T_all[ii-1,:]+0.5*k2,radius)
        
        k4 = dt * dTdt(T_all[ii-1,:]+k3,radius)
        
        T_all[ii,:] = T_all[ii-1,:] + (k1+ 2*(k2+k3) + k4)/6.
    print(T_all)
    return T_all




T_all = answer()


for ii in range(M):
    plt.plot(time,T_all[:,ii],label='M = {0}'.format(ii))
#plt.semilogy()
plt.legend()
plt.show()
#print(r_max_luke)


## KLC: Josh's To Do
## - fix dTdt() calling cooling_constant() function to accept k_r array
## - fix answer() to have all the "preamble" stuff and allow user to
##   change parameters without re-compiling
## - test range of "sensible" parameters (e.g., nshell=1)
## - add tauntaun shell (represented by diff. cooling constant array)
## - add constant or variable Hoth temperature
## - enable user to change geometry (i.e., enable "spherical" mode)
