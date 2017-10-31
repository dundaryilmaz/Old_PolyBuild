# FMKMF_F90 not set: using f90
# FMKMF_SFTAG not set: using f90
# FMKMF_SPATH not set: using . 
# FMKMF_LINKOPTS not set: using no link options 
# Main program is main.f90 
# process_fsource called with arg main.f90 
# main.f90 Uses Module statistics
# main.f90 Uses Module reax_frame
# Full list of modules in main.f90: statistics reax_frame 
# Uses statistics which is in ./statistics_module.f90
# process_fsource called with arg ./statistics_module.f90 
# Full list of modules in ./statistics_module.f90:  
# Uses reax_frame which is in ./reax_frame.f90
# process_fsource called with arg ./reax_frame.f90 
# Full list of modules in ./reax_frame.f90:  

# ------------------Macro-Defs---------------------
F90=ifort -O3 -fopenmp
LIBS=-L../Libraries
# -------------------End-macro-Defs---------------------------

# Here is the link step 
OBJECTS=ceramic_type.o checkbox.o global_variables.o main.o mcmodule.o \
		polymer_type.o polymer_module.o print_files.o userinterface.o \
		vec.o nanoparticle_module.o 
polybuild:${OBJECTS}
	 $(F90) -o polybuild ${OBJECTS} ${LIBS} -lchrpak -lgeometry
	 cp polybuild ../Test/
	 cp polybuild ~/bin/

# Here are the compile steps
 
vec.o:./vec.f90  
	 $(F90) -c ./vec.f90  -o vec.o

global_variables.o:./global_variables.f90 vec.o 
	 $(F90) -c ./global_variables.f90 -o global_variables.o 

checkbox.o:./checkbox.f90 global_variables.o  
	 $(F90) -c checkbox.f90 -o checkbox.o 
polymer_type.o:./polymer_type.f90 checkbox.o global_variables.o
	 $(F90) -c polymer_type.f90 -o polymer_type.o 
mcmodule.o:./mcmodule.f90 checkbox.o global_variables.o polymer_type.o
	 $(F90) -c mcmodule.f90 -o mcmodule.o
ceramic_type.o:./ceramic_type.f90 checkbox.o
	 $(F90) -c ceramic_type.f90 -o ceramic_type.o	
polymer_module.o:./polymer_module.f90 checkbox.o global_variables.o polymer_type.o
	 $(F90) -c polymer_module.f90 -o polymer_module.o 
nanoparticle_module.o:./nanoparticle_module.f90 global_variables.o
	 $(F90) -c nanoparticle_module.f90 -o nanoparticle_module.o 
print_files.o:./print_files.f90 global_variables.o polymer_module.o ceramic_type.o polymer_type.o
	 $(F90) -c print_files.f90 -o print_files.o 
userinterface.o:./userinterface.f90 global_variables.o polymer_module.o ceramic_type.o polymer_type.o nanoparticle_module.o
	 $(F90) -c userinterface.f90 -o userinterface.o 
main.o:./main.f90 mcmodule.o vec.o global_variables.o checkbox.o polymer_type.o ceramic_type.o nanoparticle_module.o \
				polymer_module.o print_files.o userinterface.o 
	 $(F90) -c main.f90 -o main.o 
# This entry allows you to type " make clean " to get rid of
# all object and module files 
clean:
	rm -f -r f_{files,modd}* *.o *.mod *.M *.d V*.inc *.vo \
	V*.f *.dbg album F.err
  
