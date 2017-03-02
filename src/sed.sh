list=(
force_field/force_field.cc           force_field/potential_external.cc     force_field/potential_spring.cc             molecules/bead.cc
force_field/potential_bond.cc        force_field/potential_hard_sphere.cc  force_field/potential_truncated_lj.cc       molecules/molecule.cc
force_field/potential_ewald.cc       force_field/potential_hard_wall.cc    force_field/potential_truncated_lj_wall.cc  simulation/simulation.cc
force_field/potential_ewald_coul.cc  force_field/potential_pair.cc         force_field/potential_well_wall.cc          utilities/misc.cc
force_field/force_field.h           force_field/potential_hard_sphere.h   force_field/potential_truncated_lj_wall.h  utilities/constants.h
force_field/potential_bond.h        force_field/potential_hard_wall.h     force_field/potential_well_wall.h          utilities/misc.h
force_field/potential_ewald_coul.h  force_field/potential_pair.h          molecules/bead.h
force_field/potential_ewald.h       force_field/potential_spring.h        molecules/molecule.h
force_field/potential_external.h    force_field/potential_truncated_lj.h  simulation/simulation.h
main.cc
)

for name in ${list[@]}
do
  sed "s/$1/$2/g" $name > temp
  #sed "s/totEn(/GetTotalEnergy(/g" $name > temp
  mv temp $name

done

