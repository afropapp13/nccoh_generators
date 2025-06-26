# to setup the environment
source setup_generators.sh

# to loop over the generators
root -b loop_generators.cxx

# start plotting/overlaying just the generator predictions
root -b loop_generator_overlay.cxx
root -b generator_inte_breakdown.cxx

# include the data overlay
root -b run_data_mc_overlay.cxx