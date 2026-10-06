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


P(initial binary | M1_zams, q_zams, P_zams, e_zams, Z)
= P(M1_zams)
x P(P_zams | M1_zams)
x P(q_zams | M1_zams, P_zams)
x P(e_zams | M1_zams, P_zams)

where:
q_zams = M2_zams / M1_zams

initial parameters:
M1_zams
M2_zams
q_zams
P_zams
a_zams
e_zams
Z

--- population prior: IMF x Moe & Di Stefano P's & Q's


x


P(preSN state |
  M1_zams, M2_zams, P_zams, e_zams, Z, H_binary)

preSN outputs:
M1_preSN
M2_preSN
M_He_init
M_preSN
M_CO
M_Fe
P_preSN
a_preSN
e_preSN
state1_preSN
state2_preSN

H_binary includes:
mass transfer
common envelope
winds
angular momentum loss
stellar evolution assumptions
etc.

--- binary evolution: POSYDON or another binary population synthesis model


x


P(CC outcome, M_rem,b, M_rem,g, M_fb, E_exp, t_exp |
  M_He_init, M_preSN, M_CO, M_Fe, M4, mu4, xi_2.5, H_engine)

CC outcome:
NS
fallback BH
direct-collapse BH

outputs:
M_rem,b = final baryonic remnant mass
M_rem,g = final gravitational remnant mass
M_fb = fallback mass
E_exp = explosion energy
t_exp = explosion time

for NS-forming systems:
M_NS = M_rem,g

H_engine:
W18
N20
W20
W15
S19.8
etc.

--- core-collapse/remnant model: e.g. Ertl+20


x


P(postSN orbit, bound/disrupted |
  M_preSN, M_rem,g, M2_preSN,
  a_preSN, e_preSN,
  v_k, theta_k, phi_k)

postSN outputs:
bound/disrupted
P_postSN
a_postSN
e_postSN
V_sys

relevant inputs:
M_preSN
M_rem,g
M2_preSN
DeltaM_SN
a_preSN
e_preSN
v_k
theta_k
phi_k

--- SN orbital dynamics + natal kick selection; POSYDON can perform this


x


P(present-day NS+MS state |
  M_NS, M2, P_postSN, e_postSN, age, H_binary)

require:
system remains bound
one component is an NS
companion is a luminous non-compact star
system has not merged
system has not evolved into an unwanted interacting state

present-day outputs:
M_NS
M2
P
a
e
R2
L2
T_eff,2

--- continued post-SN binary evolution to the present day


x


P(Galactic position and photometry |
  l, b, d, A_G, M2, R2, L2, T_eff,2)

outputs:
l
b
d
parallax
A_G
G_mag

--- Galactic placement + extinction + companion photometry


x


P(Gaia NSS selected |
  ra, dec, parallax,
  P, e,
  M_NS, M2,
  G_NS, G2)

for NS+MS systems:
G_NS is effectively negligible in the optical

--- Gaia DR3 NSS astrometric-orbit selection: Lam+25


x


P(NS-candidate selected |
  Gaia NSS selected,
  inferred dark-companion mass,
  orbital parameters,
  luminous-star properties,
  other candidate cuts)

--- Gaia NS+MS candidate selection: El-Badry+24





ok now for the Gaia NSS catalog selection function, it consists of 

Mass of neutron star, mass of main sequency star, orbital period, eccentricity, right ascension, declination, distance/parallax, Gaia G band mag of main sequence star, Gaia G band mag of NS -- effectively negligible for an ordinary NS
