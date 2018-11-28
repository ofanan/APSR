# Calculate allowable values of num of schedulers and num of samples per sched', which guarantee the SLA (decline rate < epsilon).
# Input: epsilon - max allowable decline rate. n-num of bins. k-num of free bins.
# Ouput: s - num of sched', d-num of samples per sched.
#        For each (s, d) output point, d is the minimal num of samples which, in conjunction with this s, still guarantee decline rate < epsilon.

import sys
import math
import numpy as np
from scipy import special

# Calculate the minimal num of probes for a sched as a func' of the desired SLA, and the other system's param's
def calc_d (s):
	nominator = math.log (epsilon - (s-1)*(1-epsilon)/k)
	denominator = math.log ( (n-k) / n)
	return int ( math.ceil(nominator / denominator) )

def calc_s_d_vecs (epsilon, k, n):	

	if (k==0 or k==n):
		print ('illegal value of k')
		exit ()

	# Smax is the maximum possible # of sched, even if each of them samples all the bins
	Smax = int ((k*epsilon) // (1-epsilon) + 1)

	d = np.array ([0] * (Smax+1)) #will hold sample sizes
	s = np.array ([0] * (Smax+1)) #will hold # of sched's
	cur_d = 0
	conf_num = -1 #number of configuration s[conf_num] and d[conf_num] will hold the s, d for this configuration.
	for i in range (1, Smax+1):
		new_d = calc_d(i)
		if (new_d != cur_d):
			conf_num += 1
			cur_d = new_d
		d[conf_num] = new_d
		s[conf_num] = i

	# Remove 0-entries     
	non_zero_indices = [i for i in range(Smax) if s[i]>0]
	d = d[non_zero_indices]
	s = s[non_zero_indices]
	return d, s

# Inputs
n = 1000; #total number of bins
k = 500;  #num of free bins
epsilon = 0.01 #acceptable decline rate
max_s_times_d = 80
d, s = calc_s_d_vecs (epsilon, k, n)
# print ('s = ', s)
# print ('d = ', d)
s_times_d = np.multiply (s, d)
# print (s_times_d)

idx_of_highest_within_th = len ( [i for i in s_times_d if i <= max_s_times_d])-1
if (idx_of_highest_within_th >=0):
  chosen_s, chosen_d = s[idx_of_highest_within_th], d[idx_of_highest_within_th]
else:
  chosen_s, chosen_d = 1, max_s_times_d
print ('Chosen s is ', chosen_s, '. Chosen d is ', chosen_d)
