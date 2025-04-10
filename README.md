# to setup the environment
source setup_generators.sh

#to loop over the generators
root -b loop_generators.cxx

root -b loop_generator_overlay.cxx

root -b generator_inte_breakdown.cxx
