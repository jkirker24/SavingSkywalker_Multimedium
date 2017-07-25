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

radius = np.linspace(0,0.310037,M+1) # +1 to make M shells
radius = radius[1:] # want radius of outer parts of shells
#above line eliminates first place 0, to avoid a division by 0
height = 2*radius


surface_area = 2*np.pi*radius*(radius + height) # brute force would have been for-loop over M

volume = np.pi*(radius**2)*height

mass_luke = 70*(volume/Vm)

mass_tt = 290*(volume/Vm)
#recreate as arrays similar to the above A(1,4) into area

#recreate cooling Constants below as a formula using a linspace or array

Ki  = (500*surface_area[0])/(mass[0]*4186)*60
Kii = (500*surface_area[1])/(mass[1]*4186)*60
Kiii = (500*surface_area[2])/(mass[2]*4186)*60
Kiv = (500*surface_area[3])/(mass[3]*4186)*60  

#not used below, intention is to eliminate the above discrete constants and use the below formula instead
cooling_constant_luke = 60*(500*surface_area)/(mass_luke*4186)
cooling_constant_tt = 60*(100*surface_area)/(mass_tt*205)


'''
    7/18/2017:
        Rewrote as single reference variables to create a list of cooling constants
    7/25/2017:
        Transformed constants Ki-Kiv into formula cooling_constant_*
        Expanded number of shells up to 10
        
    
'''

#create with same h,c values per medium
  

kconst = np.array([Ki,Kii,Kiii,Kiv])




#h*A/(m*c)*60
#k_Luke=calc_cooling_const(500,1.9,70,4186)
#k_tauntaun=calc_cooling_const(100,4.5,290,205)

time = np.linspace(tlim[0],tlim[1],N)
dt = time[1] - time[0]

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

#review above loop, 

'''
Todo:
    Create sensible file for github
    create repository on github
    decide on treating arrays with per medium constants, i.e. keep h and c the same for each array
    practice sensible variable names, comments, and units and regular sanity checks    
'''


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
    
    for ii in range(1,N):
        k1 = dt * dTdt(T_all[ii-1,:])
        
        k2 = dt * dTdt(T_all[ii-1,:]+0.5*k1)
        
        k3 = dt * dTdt(T_all[ii-1,:]+0.5*k2)
        
        k4 = dt * dTdt(T_all[ii-1,:]+k3)
        
        T_all[ii,:] = T_all[ii-1,:] + (k1+ 2*(k2+k3) + k4)/6.
    print(T_all)
    return T_all




T_all = answer()


for ii in range(M):
    plt.plot(time,T_all[:,ii],label='M = {0}'.format(ii))
#plt.semilogy()
plt.legend()
plt.show()
