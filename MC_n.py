import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


## Define the Hill function
def Hill(x, K, n):
    return x**n/(x**n+K**n)

def Heavi(x,K):
    if x>=K:
        return 1.0
    else:
        return 0.0

class Neurite():
    def __init__(self, beta, K, n=2, r=1, L_init=0):
        self.beta_base = beta
        self.beta = beta
        self.K = K
        self.n = n
        self.r_base = r
        self.r = r
        self.L= L_init
        self.L_list = [L_init,]
        self.T_list = [0,]
        
        
    def plot_L(self):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(self.T_list, self.L_list)
    
    def cal_r_n(self, L_list):
        # Change in the retraction rates is implemented here
        
        L_sum = np.sum(L_list)
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               
        # For retraction rate change
        alpha = 0.03
        beta = 0.0
        
        # No change
        alpha = 0.0
        beta = 0
        
        
        self.r = self.r_base*(1+alpha*L_sum/(1+beta*self.L))
        
    def cal_beta_n(self, L_list):
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # Change in the growth rates is implemented here
        
        L_sum = np.sum(L_list)
        
        
        # No change
        phi = 0
        mu = 0
        self.beta = self.beta_base/(1+phi*L_sum/(1+mu*self.L))
        
    def reset_L(self, L_init=0):
        # Reset the length to the given initial value
        # This is used in repeated MC simulation
        # since each time the length should start from a given initial value
        # not the value at the end of the last simulation.
        self.L = L_init
                
    
## Define the AW (actin wave) class
class AW():
    def __init__(self, rate, strength):
        self.rate_base = rate
        self.rate = rate
        self.strength_base = strength
        self.strength = strength   
        self.generated = False
  
    
    # Whether an AW is generated
    def cal_aw_generation(self, dt):
        prob_AW_generate = dt*self.rate
        rand_num = np.random.random()
        if rand_num<=prob_AW_generate:
            self.generated = True
            
    def reset_aw_generation(self):
        # After generating an AW, we need to reset self.generated to be 0
        # Otherwise, AWs would be generated forever.
        self.generated = False    
        
    def cal_neurite_enter(self, num_neurite):
        # Randomly choose a neurite to enter
        # 'high' is exclusive in the following random integer generator
        return np.random.randint(low=0, high=num_neurite)
        
         
    def cal_rate_n_neurites(self, L_list):
        # Calculate the aw rate for a neurite
        # Change of the aw rate with lengths is implemented here.
        L_sum = np.sum(L_list)
       
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        # Change the rate
        mu = 0.4
        self.rate = self.rate_base/(1+mu*L_sum)
        
        
        # No change
        mu = 0.0
        self.rate = self.rate_base/(1+mu*L_sum)
        
       
        
        
    def cal_strength_n_neurites(self, L_list):
        # Calculate the aw strength (the same for both neurites)
           
        
        L_sum = np.sum(L_list)
        
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                
        # 5 neurites
        phi = 0.1
        power = 1
        self.strength = self.strength_base/(1+phi*L_sum**power)
        
       
        
        # No change
        self.strength = self.strength_base*1
        
        # Change the strength
        phi = 0.4
        power = 1
        self.strength = self.strength_base/(1+phi*L_sum**power)

     
def MC_n(neurite_list, aw, pars):
    # The number of neurites
    num_neurite = len(neurite_list)
    # Use another name to simplify codes
    neu = neurite_list
    # Time constansts
    t_end, dt = pars
    
    for t in np.arange(0, t_end, dt):
        # Save a copy of the lengths
        # Length list
        L_list = []
        for neurite in neurite_list:
            L_list.append(neurite.L)
        
        # Update aw strength and rate according to the current lengths.
        aw.cal_strength_n_neurites(L_list)
        aw.cal_rate_n_neurites(L_list)
        #print("aw strength", aw.strength)
        for i in range(num_neurite): # Two neurites
            # Deterministic growth
            neu[i].cal_r_n(L_list) # Update the retraction rate
            neu[i].cal_beta_n(L_list) # Update the growth rate
            neu[i].L += (neu[i].beta * Hill(neu[i].L, neu[i].K, neu[i].n) - neu[i].r * neu[i].L) * dt
            # Jump in length
            aw.cal_aw_generation(dt)
            if aw.generated:
                neu[i].L += aw.strength
                aw.reset_aw_generation()
            
            # Random gaussian noise
            #d_gaussian = np.random.normal(loc=0.0, scale = 0.4*np.sqrt(dt))
            #neu[i].L += d_gaussian
            
            # Store the current length
            neu[i].L_list.append(neu[i].L)  
            neu[i].T_list.append(t)

def MC_n_ver2(neurite_list, aw, pars):
    # In this version, each aw randomly choose which neurite to enter
    # The number of neurites
    num_neurite = len(neurite_list)
    # Use another name to simplify codes
    neu = neurite_list
    # Time constansts
    t_end, dt = pars
    
    
    # Initialize t
    t = 0.0
    
    # Initialize time step
    steps = 0
  
    
    while t<=t_end:
        t += dt # Use while loop here to prevent creating large array with np.arange(0, t_end, dt)
        steps += 1
        if steps%10000==0:
            print("t=",t)
    #for t in np.arange(0, t_end, dt):
        #if t%10000==0:
        #   print("t=", t)
        # Save a copy of the lengths
        # Length list
        L_list = []
        for neurite in neurite_list:
            L_list.append(neurite.L)
        
        # Update aw strength and rate according to the current lengths.
        aw.cal_strength_n_neurites(L_list)
        aw.cal_rate_n_neurites(L_list)
               
        #print("aw strength", aw.strength)
        # Deterministic growth
        for i in range(num_neurite): 
            neu[i].cal_r_n(L_list) # Update the retraction rate
            neu[i].cal_beta_n(L_list) # Update the growth rate
            neu[i].L += (neu[i].beta * Hill(neu[i].L, neu[i].K, neu[i].n) - neu[i].r * neu[i].L) * dt
            
        # Jump in length
        aw.cal_aw_generation(dt)
        if aw.generated:
            neurite_to_enter = aw.cal_neurite_enter(num_neurite)
            neu[neurite_to_enter].L += aw.strength
            aw.reset_aw_generation()
            
            # Random gaussian noise
            #d_gaussian = np.random.normal(loc=0.0, scale = 0.4*np.sqrt(dt))
            #neu[i].L += d_gaussian
            
        # Store the current length 
        # For large time, store data only for some t values
        if t<10000 or (t%10000==0 and t>= 10000):
            for i in range(num_neurite):     
                neu[i].L_list.append(neu[i].L)  
                neu[i].T_list.append(t)
            

            
def plot_length_n(neurite_list):    
     fig = plt.figure()
     ax = fig.add_subplot(1,1,1)
     for neurite in neurite_list:
         ax.plot(neurite.T_list, neurite.L_list)
     ax.set_xlabel("time")
     ax.set_ylabel("lengths")
     
     ax.set_ylim([0,10])
     
def plot_length_2neurite_phaseplane(neurite_list):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    
    ax.plot(neurite_list[0].L_list, neurite_list[1].L_list)
    
    ax.set_xlabel(r"$l_1$")
    ax.set_ylabel(r"$l_2$")
    ax.set_xlim([-0.5,10])
    ax.set_ylim([-0.5,10])

def store_length_evolution(filename, neurite_list):
    with open(filename, 'w') as file_obj:
        N = len(neurite_list[0].L_list)

        for n in range(N):
            file_obj.write(str(neurite_list[0].T_list[n]))
            for neurite in neurite_list:
                file_obj.write(' '+str(neurite.L_list[n]))
            file_obj.write('\n')
            
    
def store_length_evolution_ver2(filename, sheet_name, neurite_list):
    # Generate a dictionary for pandas dataframe
    neurite_dict = {} 
    # Add time instants
    neurite_dict["Time"] = neurite_list[0].T_list
    # nth neurite
    n = 0
    for neurite in neurite_list:
         neurite_dict[n] = neurite.L_list
         n += 1 # Note that neurite index is from 0 to num_neurite-1
    
    # Create pandas dataframe
    df = pd.DataFrame(neurite_dict)
    
    # Write to excel file
    df.to_excel(filename, sheet_name=sheet_name, index=False)
            
           
     
     
     
     
def ratio_polarization(neurite_list, aw, pars):
    # calculate the ratio between #polarized case and #total simulation
    interval_length_polarization, num_simulation, t_end_MC, dt_MC = pars
    # Counter for the successful polarization cases
    count_polarization = 0
    # #neurites
    num_neurite = len(neurite_list)
    # Distribution of #neurites at the end of simulation
    distr_polarization = np.zeros(num_neurite+1)
    # Run MC simulation `num_simulation' times
    # Count the number of successful polarization
    # Calculate the ratio of successful polarization
    for n in range(num_simulation):
        print("n", n)
        # The number of long neurites (axons)
        num_axon = 0
        # Reset neurites' lengths
        for neurite in neurite_list:
            neurite.reset_L(L_init=0)
        # MC simulation
        MC_n_ver2(neurite_list, aw, pars=[t_end_MC, dt_MC])
        # Calculate the number of axons
        for neurite in neurite_list:
            if neurite.L>=interval_length_polarization[0] and neurite.L<=interval_length_polarization[1]:
                num_axon += 1
        # Successful polarization
        if num_axon == 1:
            count_polarization += 1
        # Count the number of cases of a specific number of axons
        distr_polarization[num_axon] += 1
        
    # Calculate the ratio of successful polarization        
    ratio_polarization = count_polarization/num_simulation
    print("ratio_polarization", ratio_polarization)
    # Nomalize the distribution
    distr_polarization = distr_polarization/num_simulation
    return ratio_polarization, distr_polarization


def reach_threshold(neurite_list, threshold):
    # Check whether all lengths have passed the threshold
    # Mark for success
    all_pass = 1 
    for neurite in neurite_list:
        if neurite.L<=threshold:
            all_pass *= 0
    return all_pass

def reset_lengths(neurite_list):
    # Reset all neurite length to be zero
    for neurite in neurite_list:
        neurite.L = 0
        
    
def FPT(neurite_list, aw, pars):
    # Estimating FPT by Monte Carlo
    num_simulation, dt, threshold = pars # number of simulations, time increment
    
    # The number of neurites
    num_neurite = len(neurite_list)
    
    # Mark for the stop of simulation
    all_pass = 0
    
    # Use another name to simplify codes
    neu = neurite_list
    
    # List to store the first passage times
    t_pass_list = []
    
    
    for n in range(num_simulation):
        
        print("Simulation", n)
        # Reset
        t = 0
        all_pass = 0
        reset_lengths(neurite_list)
        
        # Monte Carlo
        while all_pass==0:
            # Increase time
            t += dt
            
            # Length list
            L_list = []
            for neurite in neurite_list:
                L_list.append(neurite.L)
            
            # Update aw strength and rate according to the current lengths.
            aw.cal_strength_n_neurites(L_list)
            aw.cal_rate_n_neurites(L_list)
            
            # Deterministic growth
            for i in range(num_neurite): 
                neu[i].cal_r_n(L_list) # Update the retraction rate
                neu[i].cal_beta_n(L_list) # Update the growth rate
                neu[i].L += (neu[i].beta * Hill(neu[i].L, neu[i].K, neu[i].n) - neu[i].r * neu[i].L) * dt
            
            # Jump in length
            aw.cal_aw_generation(dt)
            if aw.generated:
                neurite_to_enter = aw.cal_neurite_enter(num_neurite)
                neu[neurite_to_enter].L += aw.strength
                aw.reset_aw_generation()
                
            all_pass = reach_threshold(neurite_list, threshold)
        
        # Store the time
        t_pass_list.append(t)
    
    FPT = np.mean(t_pass_list)
    print("FPT=",FPT)
    return FPT
    
     
     
            

##############################################################################################################
# MC for n neurites
# The number of neurites
num_neurite = 2
# Parameters for the neurites
beta = 10
K = np.sqrt(21)


neurite_list = []
# Generate n neurites with zero initial lengths
'''
for i in range(num_neurite):
    neurite = Neurite(beta, K, n=2, r=1, L_init=0) # By default, n=2, r=1, L_init=0
    neurite_list.append(neurite)
'''

# Generate n neurites with specified initial lengths
pos_init = [0.0,0.0]
for i in range(len(pos_init)):
    neurite = Neurite(beta, K, n=2, r=1, L_init=pos_init[i])
    neurite_list.append(neurite)



strength = 1
rate = 1
pulse_avg = strength * rate
aw = AW(rate=pulse_avg/strength, strength=strength)
# MC for n neurites
t_end_MC = 2000
dt_MC = 0.1
# The "correct" MC simulation
t_end_MC = 10**11 # For amp reduction
MC_n_ver2(neurite_list, aw, pars=[t_end_MC, dt_MC])


'''
strength = 1
rate = 1
aw = AW(rate=rate, strength=strength)
# MC for n neurites
t_end_MC = 2000
dt_MC = 0.1

# The "incorrect" MC, each neurite have its own actin wave source
MC_n(neurite_list, aw, pars=[t_end_MC, dt_MC])
'''

# Plot lengths
plot_length_n(neurite_list)
plot_length_2neurite_phaseplane(neurite_list)
# Store data for matlab plotting
#store_length_evolution("length_evolution.txt", neurite_list)
'''
filename = "MC_data/2neurites/pulseamp/length_evolution.xlsx"
sheet_name = "length_evolution"
store_length_evolution_ver2(filename, sheet_name, neurite_list)


filename = "MC_data/2neurites/pulserate/length_evolution.xlsx"
sheet_name = "length_evolution"
store_length_evolution_ver2(filename, sheet_name, neurite_list)


filename = "MC_data/2neurites/retraction/length_evolution.xlsx"
sheet_name = "length_evolution"
store_length_evolution_ver2(filename, sheet_name, neurite_list)
'''


'''
# Ratio of successful polarization
interval_length_polarization = [4,12]
num_simulation = 2000
ratio, distr = ratio_polarization(neurite_list, aw, pars=[interval_length_polarization, num_simulation, t_end_MC, dt_MC])

'''

# First passage time that all neurites reach the threshold
num_simulation = 2000 # Total number of simulations
threshold = 6.0 # Threshold to be passed by the neurite lengths
#FPT(neurite_list, aw, pars=[num_simulation, dt_MC, threshold])



    
    



            
        