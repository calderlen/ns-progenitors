import numpy as np

#P(initial binary | M1_zams;  q = M_2,zams/M_1,zams; P_zams; e_zams) = P(M_1,zams) x P(q, P, e | M1_zams) = P(M_1,zams) x P(P | M_1,zams) x P ( q | M_1,zams; P) x P( e | M_1,zams, P) --- this is the population prior, so IMF x Moe P's&Q's


def imf(M1): return	M1**-2.3

def gamma_from_q(P, q):

	logP = np.log10(P)

	gamma = 0.0
    
    if 0.1 <= q < 0.3:
    	# gamma_small: 0.1 <= q < 0.3
    	if 0.2 <= logP < 1.0: gamma = 0.1
		elif 1.0 <= logP < 3.0: gamma = 0.1 - 0.15 * (logP - 1.0)
		elif 3.0 <= logP < 5.6: gamma = -0.2 - 0.5 * (logP - 3.0)
		else: gamma = -1.5

	elif 0.3 <= q <= 1:
    	# gamma_large: 0.3 <= q <= 1
    	if 0.2 <= logP < 1.0: gamma = -0.5
		elif 1.0 <= logP < 2.0: gamma = -0.5 - 0.9 * (logP - 1.0)
		elif 2.0 <= logP < 4.0: gamma = -1.4 - 0.3 * (logP - 2.0)
		else: gamma = -2.0
        	
	return gamma

	
def period_frequency_qgeq0p3(P, M1):
	
	logP = np.log10(P)
	logM = np.log10(M1)
	
	# frequency of short-period binaries with q>0.3 to be
	f1 = 0.020 + 0.04 * logM + 0.07 * logM**2
	# frequency logP=2.7  binaries to be
    f27 = 0.039 + 0.07 * logM + 0.01 * logM**2
    # companion frequency at logP=5.5 to be
    f55 = 0.078 - 0.05 * logM + 0.04 * logM**2
	
	alpha = 0.018
    deltalogP = 0.7
	
	
	if 0.2 <= logP < 1.0:
        value = f1

    elif 1.0 <= logP < 2.7 - deltalogP:
        value = (f1 + (logP - 1.0) / (1.7 - deltalogP) * (f27 - f1 - alpha * deltalogP))

    elif 2.7 - deltalogP <= logP < 2.7 + deltalogP:
        value = f27 + alpha * (logP - 2.7)

    elif 2.7 + deltalogP <= logP < 5.5:
        value = (f27 + alpha * deltalogP + (logP - 2.7 - deltalogP) / (2.8 - deltalogP)* (f55 - f27 - alpha * deltalogP))

    elif 5.5 <= logP < 8.0:
        value = f55 * np.exp(-0.3 * (logP - 5.5))

    else:
        return 0.0
        
    # so value is the frequency of a companion mass per primary mass per unit dex in log10P, so dN/dlog10P 
        	
	return value 

def twin_fraction(M1, P):

    logM = np.log10(M1)
    logP = np.log10(P)

    F0 = 0.30 - 0.15 * logM

    if M1 <= 6.5:
        logP_twin = 8.0 - M1
    else:
        logP_twin = 1.5

    if logP < 1.0:
        return F0

    elif 1.0 <= logP < logP_twin:
        return F0 * (1.0 - (logP - 1.0) / (logP_twin - 1.0))

    else:
        return 0.0

def period_frequency_qgeq0p1(M1,P):

    logP = np.log10(P)

    if not (0.2 <= logP < 8.0):
        return 0.0

 	# the q values inserted into these functions are arbitrary, they are just to retrieve the gamma small and large needed
	gamma_small = gamma_from_q(P, 0.1)
	gamma_large = gamma_from_q(P, 1.0)
	
    # continuity at q = 0.3
    continuity = 0.3**(gamma_large - gamma_small)

    # integral from q = 0.1 to 0.3
    if gamma_small != -1.0:
    	int_small = (continuity * (0.3**(gamma_small + 1.0) - 0.1**(gamma_small + 1.0)) / (gamma_small + 1.0))
    else:
    	int_small = continuity * np.log(0.3 / 0.1)
   

    # integral from q = 0.3 to 1
    if gamma_large != -1.0:
        int_large = (1.0**(gamma_large + 1.0) - 0.3**(gamma_large + 1.0)) / (gamma_large + 1.0)
    else:
        int_large = np.log(1.0 / 0.3)

    f_twin = twin_fraction(M1, P)

    f_logP_qgeq0p3 = period_frequency_qgeq0p3(P, M1)

    return f_logP_qgeq0p3 * (1.0+ (1.0 - f_twin) * int_small / int_large)

def massratio_pdf(q, M1, P):

    logP = np.log10(P)

    if not (0.1 <= q <= 1.0):
        return 0.0

    if not (0.2 <= logP < 8.0):
        return 0.0
        
	# the q values inserted into these functions are arbitrary, they are just to retrieve the gamma small and large needed
	gamma_small = gamma_from_q(P, 0.1)
	gamma_large = gamma_from_q(P, 1.0)
	
	if gamma_large != -1.0:
		int_large = (1.0**(gamma_large + 1.0) - 0.3**(gamma_large + 1.0)) / (gamma_large + 1.0)
	else:
		int_large = np.log(1.0 / 0.3)

    # keep broken power law continuous at q = 0.3
    continuity = 0.3**(gamma_large - gamma_small)

    # integral from q = 0.1 to 0.3
    if gamma_small != -1.0:
    	int_small = (continuity * (0.3**(gamma_small + 1.0) - 0.1**(gamma_small + 1.0)) / (gamma_small + 1.0))
    else:
    	int_small = continuity * np.log(0.3 / 0.1)
    
    f_twin = twin_fraction(M1, P)	

	norm = 1+(1-f_twin)*int_small/int_large
	
    if 0.1 <= q < 0.3:
    	
        value = (1-f_twin)*continuity * q**gamma_small / (int_large*norm)
    
    elif 0.3 <= q <= 1.0:         
        if 0.95 < q < 1: p_twin = 20
     	else: p_twin = 0
        	
       	value = ((1-f_twin) * (q**gamma_large/int_large) + f_twin * p_twin)/norm

    return value
	
def eccentricity_pdf(e, P, M1):
	
	logP = np.log10(P)
	
	if P <= 2.0:
		if e == 0.0: return 1.0
		else: return 0.0
		
	e_max = 1 - (P/2)**(-2.0/3.0)
	
	if not (0.0 < e < e_max): return 0.0
		
	# late-type MS binaries
	if 0.5 < logP < 6.0 and 0.8 < M1 < 3.0:
		eta = 0.6 - 0.7/(logP - 0.5)

	# early-type MS binaries
	elif 0.5 < logP < 5.0 and M1 > 7.0:
		eta = 0.9 - 0.2/(logP - 0.5)
		
	elif 3.0 <= M1 <= 7.0:
		interpolation_weight = (M1-3)/4
		eta = (1-interpolation_weight) * (0.6 - 0.7/(logP - 0.5)) + interpolation_weight * (0.9 - 0.2/(logP - 0.5))
			
	else: return 0.0
	
	return (eta+1) * e**eta / e_max**(eta+1)


