from posydon.popsyn.synthetic_population import Population

pop = Population("1e+00_Zsun_population.h5")

hist = pop.history[0]

cols = [
    "time",
    "step_names",
    "state",
    "event",
    "S1_state",
    "S1_mass",
    "S1_he_core_mass",
    "S1_co_core_mass",
    "S2_state",
    "S2_mass",
    "orbital_period",
    "eccentricity",
]

print(hist[cols].to_string())

print(pop.number_of_systems)

# oneline String columns: 
    # state_i,
    # state_f, 
    # event_i, 
    # event_f, 
    # step_names_i, 
    # step_names_f, 
    # S1_state_i, 
    # S1_state_f, 
    # S2_state_i, 
    # S2_state_f, 
    # S1_SN_type, 
    # S2_SN_type, 
    # interp_class_HMS_HMS, 
    # interp_class_CO_HeMS, 
    # interp_class_CO_HMS_RLO, 
    # interp_class_CO_HeMS_RLO, 
    # mt_history_HMS_HMS, 
    # mt_history_CO_HeMS, 
    # mt_history_CO_HMS_RLO, 
    # mt_history_CO_HeMS_RLO.