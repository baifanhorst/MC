import numpy as np
import matplotlib.pyplot as plt


# Abbreviations
ln = np.log
exp = np.exp  



## Define the Hill function
def Hill(x, K, n):
    ratio = x/K
    if x/K>100:
        ratio = 1
    elif x/K<1e-2:
        ratio = 0
    return ratio**n/(ratio**n + 1)
    #return x**n/(x**n+K**n)

def Heavi(x,K):
    if x>=K:
        return 1.0
    else:
        return 0.0
    
class Neurite():
    def __init__(self):
        self.C = 0
        self.C_list = []
        self.L = 0
        self.L_list = []
        self.pulse_received = 0 # Number of pulses received
        self.C_diff = [] # Derivative of C
        self.C_diff_avg = [] # Derivative of C plus the avg noise term
        self.L_diff = [] # Derivative of L
    
    def initialize(self, C_init, L_init):
        self.C = C_init
        self.L = L_init
        self.C_list.append(C_init)
        self.L_list.append(L_init)
        self.C_diff.append(0)
        self.C_diff_avg.append(0)
        self.L_diff.append(0)
    
    
class Pulse():
    def __init__(self, amp, rate):
        self.amp = amp
        self.rate = rate
    
    def set_rate(self, rate):
        self.rate = rate
    def set_amp(self, amp):
        self.amp = amp
    
        
    def generate(self, dt):
        
        # Note that prob_generate is the probability of having at least one pulse
        # It is not the probability of having a single pulse
        # Therefore, the following codes work only when prob_generate is small
        prob_generate = 1 - np.exp(-dt*self.rate)
        rand_num = np.random.random()
        if rand_num<=prob_generate:
            neurite_index = np.random.randint(N)
            return [True, neurite_index]
        else:
            return [False, 0]
        
        
        '''
        prob_generate = dt*self.rate
        if prob_generate <= 1:
            rand_num = np.random.random()
            if rand_num<=prob_generate:
                neurite_index = np.random.randint(N)
                return [True, neurite_index, self.amp]
            else:
                return [False, 0, 0]
        else:
            neurite_index = np.random.randint(N)
            return [True, neurite_index, self.amp*dt*self.rate]
        '''
        
        '''
        prob_generate = dt*self.rate
        rand_num = np.random.random()
        if rand_num<=prob_generate:
            neurite_index = np.random.randint(N)
            return [True, neurite_index, self.amp]
        else:
            return [False, 0, 0]
        '''
        
class Soma():
    def __init__(self):
        self.C0 = 0
        self.C0_list = []
    
    def initialize(self, C0_init):
        self.C0 = C0_init
        self.C0_list = [C0_init,]
        
    




def time_evolution(neurites, soma, P):
    
    
    t = t_start
    num_iter = 0
    
    while t<=t_end:
        
        t += dt
        num_iter += 1
        
        C_sum = 0
        for neu in neurites:
            C_sum += neu.C
            
        
        C0 = C0max - d * C_sum
        
        
        if C0 < 0:
            C0=0
        
        soma.C0_list.append(C0)
            
        
        
        
        
        pulse_generated, pulse_neurite = P.generate(dt)
        
        for neu in neurites:
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            diff = -(neu.C-C0max)/(a*(neu.L + 1)) # No restriction on shootin1
            #diff = -(neu.C-C0)/(a*(neu.L + 1))
            neu.C += dt * diff
            neu.C_diff.append(diff)
            
        #if pulse_generated and C0>0:
        #print(pulse_generated)
        if pulse_generated:
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            #print(pulse_amp)
            
            #pulse_amp = 5 * C0
            #pulse_amp = 1.275 * C0
            pulse_amp = 18 * C0
            #pulse_amp = 1
            
            neurites[pulse_neurite].C += pulse_amp
            
        
        L_sum = 0
        for neu in neurites:
            L_sum += neu.L
        
        for neu in neurites:
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            #diff = b - ra * L_sum - rb * (1-Hill(neu.C, Kc, hc)) * Hill(neu.L, 8, 1)
            diff = b - ra * neu.L - rb * (1-Hill(neu.C, Kc, hc)) * Hill(neu.L, 8, 1)
            
            
            neu.L += dt * diff
            
            # Force a max length
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            #L_max = 75
            #if neu.L>L_max:
            #   neu.L = L_max
            
            neu.L_diff.append(diff)
           
        for neu in neurites:
            neu.C_list.append(neu.C)
            neu.L_list.append(neu.L)
        
        t_list.append(t)
            
        

        
def plot_result(neurites, t_list, soma):
    
    fig = plt.figure()
    ax = fig.add_subplot(2,1,1)
    for i in range(len(neurites)):
        ax.plot(t_list, neurites[i].C_list, label="{}".format(i))
    ax.plot(t_list, soma.C0_list, label="C0")
    ax.set_title("Concentrations")
    ax.legend(loc="upper left")
    #ax.set_ylim([-0.1, 10])
    
    
    ax = fig.add_subplot(2,1,2)
    for i in range(len(neurites)):
        ax.plot(t_list, neurites[i].L_list, label="{}".format(i))
    ax.set_title("upper left")
    ax.legend()
    #ax.set_ylim([0, 200])
    
    #ax.plot([0, np.max(t_list)], [100, 100], color="black", linestyle=":")
    
    '''
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(len(neurites)):
        ax.plot(t_list, neurites[i].C_diff_avg, label="{}".format(i))
    
    ax.set_title("C_diff_avg")
    ax.legend()
    ax.set_ylim([-100, 100])
    '''
    
def output_info(neurites, t_list, soma, P):
    for neu in neurites:
        print("Mean C", np.mean(neu.C_list))
    
   
def plot_result_onlylength(neurites, t_list):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(len(neurites)):
        ax.plot(t_list, neurites[i].L_list, label="{}".format(i))
    
    
    #ax.legend()
    
    #ax.set_xlim([0, 500])
    #ax.set_ylim([0, 120])
    
    # Hide the right and top spines
    #ax.spines['right'].set_visible(False)
    #ax.spines['top'].set_visible(False)

    # Adjust ticks
    #ax.set_yticks([0, 60, 120])
    #ax.set_xticks([0, 250, 500])
    
    # Save figure
    #filename = "./MC_figs/BMB_WTA_2D.jpg"
    #fig.savefig(filename, dpi=600)
    

# Parameter values
a = 1/8.26


pulse_rate = 1/17


b = 0.25
b = 0.4 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


ra = 0.0002

rb = 0.5

Kc = 2
hc = 5

C0max = 2
C0max = 0.56 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



d = 0.1


N = 2 # Number of neurites
neurites = [] # List of neurites
for i in range(N):
    neu = Neurite()
    neurites.append(neu)

# Initialize a long and a short neurite
neurites[0].initialize(C_init=0, L_init=20)
neurites[1].initialize(C_init=0, L_init=20)
# Generate a pulse object
# The amplitude is not important here, since it is 5*C0, updated during iteration
# Thus, amp is set to 1 for convenience
P = Pulse(amp=1, rate=pulse_rate)

# Initialize C0
soma = Soma()
C_sum = 0
L_sum = 0
for neu in neurites:
    C_sum += neu.C
    L_sum += neu.L

C0 = C0max - d * C_sum
if C0 < 0:
    C0=0
soma.initialize(C0_init=C0)




t_start = 0
dt = 0.02
t_end = 30000
t_list = [t_start,]




time_evolution(neurites, soma, P)
plot_result(neurites, t_list, soma) 
output_info(neurites, t_list, soma, P)

#plot_result_onlylength(neurites, t_list)
 
    
 


        



    
    



            
        