import os, sys, argparse
import math

def shannons_alpha(p):
	# shannons = -sum(pi ln(pi))
	h = []
	for i in p:
		h.append(i * math.log(i))
	print("Shannon's diversity: %s" %(-1 *sum(h)))
	return (-1 *sum(h))

def berger_parkers_alpha(p):
	# bp is nmax/N which is equal to the max pi == max(p)
	print("Berger-parker's diversity: %s" %max(p))
	return max(p)


def fischers_alpha():	
	global np
	import numpy as np
	from scipy.optimize import fsolve
	
	zGuess = np.array([1,1])
	fish = fsolve(eqn_output,zGuess)
	# NOTE:
	# if ratio of N/S > 20 then x > 0.99 (Poole,1974)
	# x is almost always > 0.9 and never > 1.0 i.e. ~ 0.9 < x < 1.0
	print("Fischer's index: %s" %fish[1])
	return fish[1]

def eqn_output(z):
	x = z[0]
	a = z[1]

	F = np.empty((2))
	F[0] = np.exp(((x-1)*(N_f/S_f))/x) + x - 1
	F[1] = ((N_f*(1-x))/x) - a
	return F

# Main method
def main():
	# get arguments
	parser = argparse.ArgumentParser(description='Pick an alpha diversity.')
	parser.add_argument('-f','--filename',dest='filename',help='bracken file with species abundance estimates')
	parser.add_argument('-a','--alpha',dest='value',default='S',type=str,  help='type of alpha diversity to calculate S, BP or F, default = S')
	args = parser.parse_args()

	f = open(args.filename)
	
	f.readline()
	n = []
	# read in the file
	for line in f: 
		ind_abund = line.split('\t')[5] # finds the abundance estimate
		n.append(float(ind_abund)) 
	
	f.close()
	
	# calculations
	N = sum(n) # total number of individuals
	S = len(n) # total number of species
	# calculate all the pi's
	p = [] # store pi's 
	for i in n: # go through each species
		if i != 0:
			p.append(i/N) # p is the ni/N
	
	
	# find the indicated alpha
	if args.value == 'S': # calculate shannon's diversity
		shannons_alpha(p)
	elif args.value == 'BP': # calculate berger-parker's dominance index
		berger_parkers_alpha(p)
	elif args.value == 'F': # calculate fischer's alpha
		print("Fischer's alpha...loading")
		global N_f
		N_f = sum(n)
		global S_f
		S_f = len(n)
		fischers_alpha()
	else:
		print("Not a supported alpha")

if __name__ == "__main__":
    main()

