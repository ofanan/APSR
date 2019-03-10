# Input: n, max_num_of_scheds, and budget, 
# Output: a full table of recommended (s,d) for every value of k, where 1<=k<=n

# Calculate the exact epsilons, using our combinatorial balls-and-bins analysis
import math
import numpy as np
from scipy import special
from scipy.stats import binom

# Calculate P(a,b) = \sum_{1<=i<=b, 0<=X_i<=a, \sum X_i<=b} \prod (i^{X_i})
def calc_P (a, b):
	# Initialize P is a 2D-array of 1. P[a][b] will hold P(a,b)
	# As the table starts from 0, it should be actually of size (a+1)*(b+1)
	a = a+1 
	b = b+1
	P = [x[:] for x in [[0] * (b+1)] * (a+1)] 
	for cur_b in range (0, b+1):
		P[0][cur_b] = 1
	for cur_a in range (0, a+1):
		P[cur_a][1] = 1
	for cur_b in range (2, b+1):
		for cur_a in range (1, a+1):
			for j in range (0, cur_a+1):
				P[cur_a][cur_b] += pow (cur_b, j) * P[cur_a-j][cur_b-1]
	return P

def Print_P (P, Smax):
	i = 0
	for i in range (0, Smax + 1):
		print (P[i])
	
# Calculate Pr (Hs=h | Fs=f)
def calc_Pr_h_cond_f (k, h, f, P):
	if (f > k or h > f): 
		return 0
	return special.comb (k, h) * math.factorial (h) * P[f-h][h] / pow (k, f)	

# Calculate E[Hs | Fs=f]
def calc_E_H_s_cond_f (k, f, P):
	if (f > k):
		print ('Error: Cannot have more than k potentially happy agents')
		exit (0)
	E_H_s_cond_f = 0
	for h in range (1, f+1): 
		E_H_s_cond_f += h * special.comb (k, h) * math.factorial (h) * P[f-h][h];
	return E_H_s_cond_f / pow (k, f);
	
# Calculate epsilon for a sys with given n,k,s,d,P
def calc_epsilon (n, k, s, d, P):
	if (k > n or s > n):
		print ('Wrong inputs to calc_epsilon: n = ', n, ' k = ', k, ' s = ', s, ' d = ', d)
		exit (0)
	sigma = 1 - pow ( (n-k)/n, d)
	E_H_s = 0
	for f in range (1, s+1):
		E_H_s += binom.pmf(f, s, sigma) * calc_E_H_s_cond_f (k, f, P)
	return 1 - E_H_s / s

# For debug only: an accessory func'. Used for estimating the relative contribution of each item in the sum over h to E[Hs]
def check_E_H_s_cond_f ():
	k = 100
	f = 5
	P = calc_P (f, f)
	if (f > k):
		print ('Error: Cannot have more than k potentially happy agents')
		exit (0)
	E_H_s_cond_f = 0
	for h in range (1, f+1): 
		E_H_s_cond_f += (h * special.comb (k, h) * math.factorial (h) * P[f-h][h]) / pow (k, f)
		print (f'cur_P = {P[f-h][h]:}')
		print (f'E_Hs = {E_H_s_cond_f:.4f}')
	
# Main function, to be used when called from another file.
# Inputs: n, k, maximal s to check (Smax), budget and target epsilon defined by the SLA
# The function loops over 1<=s<=Smax, until finding the maximal s which still satisfies the SLA and the budget.
# The d used is floor (budget / s).
# The func' assumes that epsilon (n, k, s) is non-dec. in s.
def MaxParalForAllK (n, Smax, budget, target_epsilon):
	print ('Combinatorial analysis: n = ', n,', budget = ', budget, ', max allowed epsilon = ', target_epsilon)
	print ('******************************************************************')
	print ('Format of the TLB below is:\nk s d estimated_decline_ratio [Note]\n')
	
	s = 0 # default value (indicating a failure in finding values which satisfy the SLA & budget constraints)55
	P = calc_P(Smax, Smax)
	
	for k in range (1, n+1):
		while (s <= Smax):
			calculated_epsilon = calc_epsilon (n, k, s+1, budget // (s+1), P)
			if (calculated_epsilon > target_epsilon):
				break
			s = s + 1
		
		if (s == 0):		
			print (f'{k} {s} {budget // (s+1)} {calculated_epsilon:.4f}  Cannot provide SLA guarantees even with s=1')
			continue #No guarantees for current value of k; skip to k++
			
		# Now we know that s>0, 
		# and s holds the maximal # of sched's which can be used at the given budget.
		
		# Minimize the d, while keeping SLA guarantees
		d = budget // s
		while (1):
			calculated_epsilon = calc_epsilon (n, k, s, d-1, P)
			if (calculated_epsilon > target_epsilon):
				break
			d = d - 1
			
		calculated_epsilon = calc_epsilon (n, k, s, d, P)
		print (f'{k} {s} {d} {calculated_epsilon:.4f}')
		
	

def Debug (n, k, Smax, budget, target_epsilon)	:
	P = calc_P(Smax, Smax)
	line = 20
	print (P[line])
	# Print_P (P, Smax)
		
# Main func', executed when only this file is run.
# Defines parameters' values, and calls MaxParal.		
def main ():
	n = 1000 #total number of bins
	max_num_of_scheds = 20 #Still need some maximum, as otherwise the huge number cause problems at the binom etc.
	budget = n
	target_epsilon = 0.1

	MaxParalForAllK (n, max_num_of_scheds, budget, target_epsilon)	
	
if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
   main()	
   		
