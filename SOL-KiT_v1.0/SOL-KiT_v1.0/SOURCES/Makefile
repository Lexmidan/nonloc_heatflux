
#  Copyright 2020 Stefan Mijin

#  This file is part of SOL-KiT.

#      SOL-KiT is free software: you can redistribute it and/or modify
#      it under the terms of the GNU General Public License as published by
#      the Free Software Foundation, either version 3 of the License, or
#      (at your option) any later version.

#      SOL-KiT is distributed in the hope that it will be useful,
#      but WITHOUT ANY WARRANTY; without even the implied warranty of
#      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#      GNU General Public License for more details.

#      You should have received a copy of the GNU General Public License
#      along with SOL-KiT.  If not, see <https://www.gnu.org/licenses/>.

# Start of the makefile

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# Defining variables
objects = main.o input.o sparse.o printing.o normalization.o switches.o neut_and_heat.o grid.o \
					f_init.o moments.o output.o tests.o post_processing.o solver_params.o build_coulomb_coll.o \
					coll_cs_rates.o inel_grid.o build_neutral_coll.o build_crm.o build_0d.o build_x_adv.o \
					build_e_adv.o build_vlasov.o build_maxwell.o build_heating.o build_neutral_diff.o build_div_boundary.o \
					evolve.o mat_build.o ext_solver.o matrix_data.o build_ion_fluid.o build_particle_source.o variable_kind_data.o build_fluid_el.o
f90comp = mpif90
flags = -Wall  -O3 -pedantic
cdirs = -cpp
switch =
incl = -I${PETSC_DIR}/include/ -I${PETSC_DIR}/${PETSC_ARCH}/include/
lib = $(PETSC_LIB)
# Makefile
SOL-KiT: $(objects)
	$(CLINKER) $(flags) -o ../SOL-KiT $(switch)  $(objects) $(lib)

var_data.mod: var_data.o variable_kind_data.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) variable_kind_data.f90
var_data.o: variable_kind_data.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) variable_kind_data.f90

input.mod: var_data.mod input.o input.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) input.f90
input.o: var_data.mod input.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) input.f90

sparse.mod: var_data.mod sparse.o sparse.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) sparse.f90
sparse.o: var_data.mod sparse.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) sparse.f90

printing.mod: var_data.mod printing.o printing.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) printing.f90
printing.o: var_data.mod printing.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) printing.f90

normalization.mod: var_data.mod input.mod normalization.o normalization.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) normalization.f90
normalization.o: var_data.mod input.mod normalization.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) normalization.f90

switches.mod: var_data.mod input.mod switches.o switches.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) switches.f90
switches.o: var_data.mod input.mod switches.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) switches.f90

solver_params.mod: var_data.mod input.mod solver_params.o solver_params.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) solver_params.f90
solver_params.o: var_data.mod input.mod solver_params.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) solver_params.f90

neut_and_heat.mod: var_data.mod switches.mod input.mod neut_and_heat.o neut_and_heat.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) neut_and_heat.f90
neut_and_heat.o: var_data.mod switches.mod input.mod neut_and_heat.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) neut_and_heat.f90

grid.mod: var_data.mod normalization.mod switches.mod input.mod neut_and_heat.mod grid.o grid.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) grid.f90
grid.o: var_data.mod normalization.mod switches.mod input.mod neut_and_heat.mod grid.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) grid.f90

f_init.mod: var_data.mod normalization.mod switches.mod input.mod neut_and_heat.mod grid.mod f_init.o f_init.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) f_init.f90
f_init.o: var_data.mod normalization.mod switches.mod input.mod neut_and_heat.mod grid.mod f_init.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) f_init.f90

moments.mod: var_data.mod normalization.mod grid.mod moments.o moments.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) moments.f90
moments.o: var_data.mod normalization.mod grid.mod moments.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) moments.f90

output.mod: var_data.mod f_init.mod solver_params.mod switches.mod normalization.mod neut_and_heat.mod grid.mod output.o output.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) output.f90
output.o: var_data.mod f_init.mod solver_params.mod switches.mod normalization.mod neut_and_heat.mod grid.mod output.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) output.f90

post_processing.mod: var_data.mod evolve.mod build_div_boundary.mod normalization.mod switches.mod neut_and_heat.mod grid.mod moments.mod output.mod tests.mod printing.mod post_processing.o post_processing.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) post_processing.f90
post_processing.o: var_data.mod evolve.mod build_div_boundary.mod normalization.mod switches.mod neut_and_heat.mod grid.mod moments.mod output.mod tests.mod printing.mod post_processing.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) post_processing.f90

build_x_adv.mod: var_data.mod sparse.mod switches.mod grid.mod build_x_adv.o build_x_adv.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_x_adv.f90
build_x_adv.o: var_data.mod sparse.mod switches.mod grid.mod build_x_adv.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_x_adv.f90

build_maxwell.mod: var_data.mod switches.mod sparse.mod normalization.mod grid.mod build_maxwell.o build_maxwell.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_maxwell.f90
build_maxwell.o: var_data.mod switches.mod sparse.mod normalization.mod grid.mod build_maxwell.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_maxwell.f90

inel_grid.mod: var_data.mod normalization.mod grid.mod neut_and_heat.mod inel_grid.o inel_grid.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) inel_grid.f90
inel_grid.o: var_data.mod normalization.mod grid.mod neut_and_heat.mod inel_grid.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) inel_grid.f90

matrix_data.mod: var_data.mod switches.mod build_x_adv.mod neut_and_heat.mod sparse.mod grid.mod build_maxwell.mod inel_grid.mod matrix_data.o matrix_data.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) matrix_data.f90
matrix_data.o: var_data.mod switches.mod build_x_adv.mod neut_and_heat.mod sparse.mod grid.mod build_maxwell.mod inel_grid.mod matrix_data.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) matrix_data.f90

ext_solver.mod: var_data.mod grid.mod matrix_data.mod sparse.mod solver_params.mod ext_solver.o ext_solver.F90
	$(f90comp) -c $(flags) $(cdirs) $(switch) $(incl) ext_solver.F90
ext_solver.o: grid.mod matrix_data.mod sparse.mod solver_params.mod ext_solver.F90
	$(f90comp) -c $(flags) $(cdirs) $(switch) $(incl) ext_solver.F90

coll_cs_rates.mod: var_data.mod matrix_data.mod inel_grid.mod normalization.mod grid.mod neut_and_heat.mod coll_cs_rates.o coll_cs_rates.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) coll_cs_rates.f90
coll_cs_rates.o: var_data.mod matrix_data.mod inel_grid.mod normalization.mod grid.mod neut_and_heat.mod coll_cs_rates.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) coll_cs_rates.f90

tests.mod: var_data.mod coll_cs_rates.mod neut_and_heat.mod switches.mod normalization.mod grid.mod tests.o tests.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) tests.f90
tests.o: var_data.mod coll_cs_rates.mod neut_and_heat.mod switches.mod normalization.mod grid.mod tests.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) tests.f90

build_coulomb_coll.mod: var_data.mod matrix_data.mod switches.mod normalization.mod grid.mod build_coulomb_coll.o build_coulomb_coll.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_coulomb_coll.f90
build_coulomb_coll.o: var_data.mod matrix_data.mod switches.mod normalization.mod grid.mod build_coulomb_coll.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_coulomb_coll.f90

build_ion_fluid.mod: var_data.mod build_div_boundary.mod neut_and_heat.mod coll_cs_rates.mod build_coulomb_coll.mod matrix_data.mod switches.mod normalization.mod grid.mod build_ion_fluid.o build_ion_fluid.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_ion_fluid.f90
build_ion_fluid.o: var_data.mod build_div_boundary.mod neut_and_heat.mod coll_cs_rates.mod build_coulomb_coll.mod matrix_data.mod switches.mod normalization.mod grid.mod build_ion_fluid.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_ion_fluid.f90

build_fluid_el.mod:  var_data.mod build_div_boundary.mod neut_and_heat.mod coll_cs_rates.mod matrix_data.mod switches.mod normalization.mod grid.mod build_ion_fluid.o build_fluid_el.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_fluid_el.f90
build_fluid_el.o: var_data.mod build_div_boundary.mod neut_and_heat.mod coll_cs_rates.mod matrix_data.mod switches.mod normalization.mod grid.mod build_fluid_el.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_fluid_el.f90

build_neutral_coll.mod: var_data.mod matrix_data.mod inel_grid.mod coll_cs_rates.mod normalization.mod grid.mod neut_and_heat.mod build_neutral_coll.o build_neutral_coll.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_neutral_coll.f90
build_neutral_coll.o: var_data.mod matrix_data.mod inel_grid.mod coll_cs_rates.mod normalization.mod grid.mod neut_and_heat.mod build_neutral_coll.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_neutral_coll.f90

build_crm.mod: var_data.mod inel_grid.mod matrix_data.mod coll_cs_rates.mod normalization.mod grid.mod neut_and_heat.mod build_crm.o build_crm.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_crm.f90
build_crm.o: var_data.mod inel_grid.mod matrix_data.mod coll_cs_rates.mod normalization.mod grid.mod neut_and_heat.mod build_crm.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_crm.f90

build_0d.mod: var_data.mod switches.mod grid.mod build_coulomb_coll.mod build_neutral_coll.mod build_crm.mod build_0d.o build_0d.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_0d.f90
build_0d.o: var_data.mod switches.mod grid.mod build_coulomb_coll.mod build_neutral_coll.mod build_crm.mod build_0d.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_0d.f90

build_e_adv.mod: var_data.mod matrix_data.mod switches.mod grid.mod build_e_adv.o build_e_adv.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_e_adv.f90
build_e_adv.o: var_data.mod matrix_data.mod switches.mod grid.mod build_e_adv.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_e_adv.f90

build_vlasov.mod: build_maxwell.mod matrix_data.mod build_x_adv.mod var_data.mod build_e_adv.mod switches.mod grid.mod build_vlasov.o build_vlasov.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_vlasov.f90
build_vlasov.o: build_maxwell.mod matrix_data.mod build_x_adv.mod var_data.mod build_e_adv.mod switches.mod grid.mod build_vlasov.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_vlasov.f90

build_particle_source.mod: var_data.mod neut_and_heat.mod switches.mod matrix_data.mod normalization.mod grid.mod build_heating.o build_particle_source.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_particle_source.f90
build_particle_source.o: var_data.mod neut_and_heat.mod switches.mod matrix_data.mod normalization.mod grid.mod build_particle_source.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_particle_source.f90

build_heating.mod: var_data.mod neut_and_heat.mod switches.mod matrix_data.mod normalization.mod grid.mod build_heating.o build_heating.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_heating.f90
build_heating.o: var_data.mod neut_and_heat.mod switches.mod matrix_data.mod normalization.mod grid.mod build_heating.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_heating.f90

build_neutral_diff.mod: var_data.mod normalization.mod neut_and_heat.mod switches.mod matrix_data.mod grid.mod build_neutral_diff.o build_neutral_diff.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_neutral_diff.f90
build_neutral_diff.o: var_data.mod normalization.mod neut_and_heat.mod switches.mod matrix_data.mod grid.mod build_neutral_diff.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_neutral_diff.f90

build_div_boundary.mod: var_data.mod coll_cs_rates.mod normalization.mod neut_and_heat.mod switches.mod matrix_data.mod grid.mod solver_params.mod build_div_boundary.o build_div_boundary.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_div_boundary.f90
build_div_boundary.o: var_data.mod coll_cs_rates.mod normalization.mod neut_and_heat.mod switches.mod matrix_data.mod grid.mod solver_params.mod build_div_boundary.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) build_div_boundary.f90

mat_build.mod: build_fluid_el.mod var_data.mod build_particle_source.mod build_ion_fluid.mod build_coulomb_coll.mod build_div_boundary.mod build_heating.mod build_maxwell.mod build_vlasov.mod build_0d.mod printing.mod switches.mod grid.mod build_neutral_diff.mod mat_build.o mat_build.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) mat_build.f90
mat_build.o: build_fluid_el.mod var_data.mod build_particle_source.mod build_ion_fluid.mod build_coulomb_coll.mod build_div_boundary.mod build_heating.mod build_maxwell.mod build_vlasov.mod build_0d.mod printing.mod switches.mod grid.mod build_neutral_diff.mod mat_build.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) mat_build.f90

evolve.mod: var_data.mod build_div_boundary.mod coll_cs_rates.mod matrix_data.mod sparse.mod ext_solver.mod solver_params.mod mat_build.mod switches.mod neut_and_heat.mod grid.mod moments.mod printing.mod evolve.o evolve.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) evolve.f90
evolve.o: var_data.mod build_div_boundary.mod coll_cs_rates.mod matrix_data.mod sparse.mod ext_solver.mod solver_params.mod mat_build.mod switches.mod neut_and_heat.mod grid.mod moments.mod printing.mod evolve.f90
	$(f90comp) -c $(flags) $(cdirs) $(switch) evolve.f90

main.o: var_data.mod ext_solver.mod build_coulomb_coll.mod matrix_data.mod output.mod post_processing.mod normalization.mod f_init.mod evolve.mod sparse.mod coll_cs_rates.mod inel_grid.mod solver_params.mod switches.mod neut_and_heat.mod grid.mod printing.mod main.F90
	$(f90comp) -c $(flags) $(cdirs) $(switch) $(incl) main.F90

# Cleaning everything
cleanall:
	rm -f *.mod
	rm -f *.o
	rm -f ../SOL-KiT

# Cleaning only .mod and .o files, retain exec
cleancomp:
	rm -f *.mod
	rm -f *.o
# End of the makefile
