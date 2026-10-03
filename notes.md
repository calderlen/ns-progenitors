todo

you need to compute bayes factors




(1) two main-sequence stars
(2) primary evolves off of MS / expands
(3) RLOF or CE
(4) H envelope stripped
(5) He star + MS
(6) Advanced He star burning  / pre-SN structure
(7) CCSNe
(8) NS formation + fallback
(9) SN mass loss + natal kick
(10) bound NS+MS survivor
(11) gaia detection/seleciton


Ertl catalog covers steps (5)-(8)

probability that a progenitor system would evolve into the type of system gaia would observe


P(initial binary | M1_zams;  q = M_2,zams/M_1,zams; P_zams; e_zams) = P(M_1,zams) x P(q, P, e | M1_zams) = P(M_1,zams) x P(P | M_1,zams) x P ( q | M_1,zams; P) x P( e | M_1,zams, P) --- this is the population prior, so IMF x Moe P's&Q's
x
P(M_He, M_2, a_preSN, ... | M_1,zams, M_2,zams, a_zams) --- binary evolution, maybe POSYDON or some other binary evolution grid/code
x
P(M_NS | M_He, H_explosion) --- remnant model, Ertl
x
P(bound | M_preSN, M_NS, M_2, a_preSN, v_k) --- physical selection, maybe need to do this on your own?
x
P(Gaia selection | M_NS, M_2, P, e, d, mag, ...) --- observational selection, probably don't need to do this on your own, but maybe.



their model is build from five fitted quantities: (f_{logP;q>0.3}, gamma_smallq, gamma_largeq, F_twin, eta) all functions of M_1 and P.



P(M_1,zams) \propto case 1: M^-2.3 for M > 0.5 M_sun, and case 2:M^-1.3 for 0.08 M_sun < M < 0.5 M_sun -- this is the Chabrier IMF

but obviously the IMF neq to today's galactic field mass distribution, because t_ms ~ 100Gyr*(M/M_sun)^-2.5

so the present day mass function PDMF = IMF x SF history x stellar lifetimes x dynamical evolution



then P's and Q's paper gives 

P(q,P,e | M_1,zams)

