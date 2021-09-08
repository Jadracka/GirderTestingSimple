import numpy as np

#Mx is a 2D numpy array


def degeneracy_doctor(Mx, thresh = 0.005):
	#This looks through the rows and columns of matrix Mx and checks if rows or cols are nearly linearly independent. 
	#row/columns p,q are considered linearly dependent if 
	#  p dot q / |q||p| > 1.0 - thresh 
	#Inputs: 
	#  Mx : a numpy 2D array Mx
	#  thresh : the degree of similarity at which rows or columns are 
	#           considered approxomatly degenerate
	#Output: 
	# bool = true iff all rows and columns are different and matrix is full rank,

	all_diff = True
	N,M = Mx.shape

	#print("Examining matrix rank...")
	#row loop
	for n in range(N):
		r = Mx[n,:]
		for nn in range(n+1,N):
			s = Mx[nn,:]
			metric = r.dot(s)/(np.linalg.norm(r) * np.linalg.norm(s)) 
			if metric > (1.0 - thresh):
				print(f"    row {n} resembles row {nn}. (metric = {metric}, s/r scale: {np.linalg.norm(s) / np.linalg.norm(r)})")
				all_diff = False
				
	#col loop
	for n in range(M):
		r = Mx[:,n]
		for nn in range(n+1,M):
			s = Mx[:,nn]
			metric = r.dot(s)/(np.linalg.norm(r) * np.linalg.norm(s)) 
			if metric > (1.0 - thresh):
				print(f"    col {n} resembles col {nn}. (metric = {metric}, s/r scale: {np.linalg.norm(s) / np.linalg.norm(r)})")
			
				all_diff = False
	if all_diff:
		print("    Matrix appears to have full rank")
	return all_diff
		 

