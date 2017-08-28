import numpy as np
import matplotlib.pyplot as plt

def T_Hoth(t):
    global Tmin_H, day_H
    return 0.5*Tmin_H*(1+np.sin(2*np.pi*t/day_H))

def shape(geometric_solid, radius, Radius_Max, Total_Mass):
    height = Radius_Max + radius
    
    if sphere:
        Surface_Area = 2*np.pi*radius*(radius+height)
        Volume = (4/3)*np.pi*(radius**3)
        Volume_Max = (4/3)*np.pi*(Radius_Max**3)
        Mass = Total_Mass*(Volume/Volume_Max)
        return Volume, Surface_Area, Mass
    
    elif rounded_cylinder:
        Surface_Area = 2*np.pi*radius*(radius+height)
        Volume = np.pi*(radius**2)*height
        Volume_Max = np.pi*(Radius_Max**2)*height
        Mass = Total_Mass*(Volume/Volume_Max)

        
        return Volume, Surface_Area, Mass
    
    
'''
    current idea is that everything here has to be an input variable, and when function is called in
    cooling constant, that values will have to be stated there, allowing for individuality.  Question
    remains as how to assign origin/base values.
'''
    
'''
    Intention to call the above function within cooling_constant()
    problems:
        -if statements have to indicate which shape is being used, source indicated either\
        within cooling_constant() or answer()
        -for rr loop is used in shape() but stated in cooling_constant
            not certain if shape function can use a loop written in cooling_constant()
        
'''    
def cooling_constant_computation(radius):## KLC: you need to pass in a radius variable
    global r_max_tt, r_max_luke, volume_max_luke, volume_max_tt
    
    cooling_constant = np.empty(len(radius)) ## cooling constant bin
    volume_max_luke = ((3*r_max_luke)/(4*np.pi))**(1/2)
    
# turn if statement to be in regards to shape, where shape value is given in answer()
# have to incorporate radius limits


    for rr in range(len(radius)):
        height = r_max_luke + radius[rr]
        if radius[rr] <= r_max_luke:
            Surface_Area = 2*np.pi*radius[rr]*(radius[rr]+height)
            Volume = (4/3)*np.pi*(radius[rr]**3)
            Mass = 70*(Volume/volume_max_luke)
            cooling_constant[rr] = 60*(500*Surface_Area)/(Mass*4186)
        elif radius[rr] > r_max_luke and radius[rr] < r_max_tt:
            Surface_Area = 2*np.pi*radius[rr]*(radius[rr]+height)
            Volume = np.pi*(radius[rr]**2)*height
            Mass = 290*(Volume/volume_max_tt)
            cooling_constant[rr] = 60*(100*Surface_Area)/(Mass*205)
    return cooling_constant

'''
    above:
        the above has a for(radius) and an if based on the radius.
        i.e.
    for rr in range(len(radius)):
        if shape
'''

def dTdt(T_now,T_H_now,Tmin_H = -60):
    global k_LS, k_tt, radius
    
    nshells = len(T_now)
    dTdt_now = np.empty(nshells)
    cooling_const = cooling_constant_computation(radius) ## KLC: but really, pass in or make global variable e.e., k_r ****
    dTdt_now[0] = -cooling_const[0]*(T_now[0] - T_now[1])
    for ss in range(1,nshells-1): # up to but not including outer shell
#        print('ss = {0}'.format(ss))
        dTdt_now[ss] = -cooling_const[ss]*(T_now[ss] - T_now[ss+1])
    ## Here would be another for-loop or if-statement to handle next object
    dTdt_now[-1] = -cooling_const[-1]*(T_now[-1] - Tmin_H)# - T_now[4]) error in T_now[4]: 
    return dTdt_now            
    

def answer(RM_tt = 0.310037, RM_LS = (1.9/(4*np.pi))**(1/2), \
           VM_tt = 0.187249, Ti_LS = 35., Ti_tt = 37., \
           day_H = 24, tlim = (0.,24.), Tmin_H = -60, \
           N=10000, M=10):
    
    '''
    values for radii max to set global variables
    volume max of tauntaun set to global variables
    volume max of luke defined in cooling_constant_computation
        because it uses r_max_luke    
    '''

    global r_max_tt, r_max_luke, volume_max_tt, radius
    
    r_max_tt = RM_tt
    r_max_luke = RM_LS
    volume_max_tt = VM_tt

    radius = np.linspace(0,0.310037,M+2) # +1 to make M shells:
        #above changed to M+2 in order to solve index of M=1 shells
    radius = radius[1:] # 
    
    time = np.linspace(tlim[0],tlim[1],N+1)
    dt = time[1] - time[0]
    
    T_all = np.empty((N+1,M+1),float)
    T_inner = 37. # C; human core temp
    T_outer = 32. # C; human freezing
    slope = (T_outer - T_inner)/radius[-1] # C/m where origin at radius = 0 
    T_init = slope*radius + T_inner # C

#T_all[0,:] = np.array([37,35,33,31],float) # arbitrary instantiation
    T_all[0,:] = T_init 

    
    ## KLC: make all the floating "preamble" definitions part of this
    ## function
    for ii in range(1,N+1):
        k1 = dt * dTdt(T_all[ii-1,:],radius)
        
        k2 = dt * dTdt(T_all[ii-1,:]+0.5*k1,radius)
        
        k3 = dt * dTdt(T_all[ii-1,:]+0.5*k2,radius)
        
        k4 = dt * dTdt(T_all[ii-1,:]+k3,radius)
        
        T_all[ii,:] = T_all[ii-1,:] + (k1 + 2*(k2+k3) + k4)/6.
        
    for aa in range(M):
        plt.plot(time,T_all[:,aa],label='M = {0}'.format(aa+1))
#plt.semilogy()

    plt.legend()
    plt.show()
    print('shape:',T_all.shape) ##print is N+1, M+1 to allow N=M=1
    return T_all


T_all = answer()





## KLC: Josh's To Do

##8/23
  ##  Analysis of Original code
    
      
##XX - fix dTdt() calling cooling_constant() function to accept k_r array
##XX - fix answer() to have all the "preamble" stuff and allow user to
##   change parameters without re-compiling

## - test range of "sensible" parameters (e.g., nshell=1)
##XX : both N = 1 and M = 1 work

## - add tauntaun shell (represented by diff. cooling constant array) 

## - add constant or variable Hoth temperature
## : Implemented Hoth_T for constant, need to change answer() calculations and test
    
## - enable user to change geometry (i.e., enable "spherical" mode)
## : currently trying to create a function to state shape
        ## not sure how to associate radius into if statements of shape
            ## i.e. input should be:
                ##"if shape = shpere:
                    ##code for sphere"

'''
#    dTdt_now[0] = -cooling_constant_computation[0]*(T_now[0] - T_now[1])

    ## KLC: functions are accessed with parentheses, not square
    ## brackets but it seems you actually want an array of cooling
    ## constants so since the radius variable is defined and
    ## instantiated outside any function, it is available to all
    ## functions. So we'll use the function
    ## cooling_constant_computation() on the "global-ish" variable
    ## radius to define an array of cooling constants. Then use it in
    ## the calculations.

#T_Hoth = 0.5*Tmin_H*(1+np.sin(2*np.pi*time/day_H))
'''