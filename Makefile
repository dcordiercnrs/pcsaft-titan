# ----------------------------------------------------------------------------------------------------------------------------------
#
#                                  Makefile la mise au point de la nouvelle version de l'quation PC-SAFT
#                                              avec calcul des drives par rapport  T et rho
#                                             (permet d'en dduire les Cp et la vitesse du son) 
#
# ----------------------------------------------------------------------------------------------------------------------------------
# Daniel Cordier, october 2013, Institut UTINAM, France.
#                   avril 2020, GSMA, Reims, France.
#
# - contact: daniel.cordier@univ-reims.fr
#         
# ----------------------------------------------------------------------------------------------------------------------------------
# On se renseigne sur la machine et l'OS :
MACHINE = $(shell hostname)
OS      = $(shell uname)
#
#***********************************************************************************************************************************
#
FC = gfortran-8

# Options pour 'gfortran' :
OPT= -fimplicit-none -fbacktrace -O3 -finit-real=zero -finit-integer=0 -fmax-errors=3 -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -Wunused-parameter
#
all: MESSAGE demo_pcsaft demo_binary
#
MAIN_1       = pcsaft_demo
MAIN_2       = binary_N2CH4_demo
#
# ----------------------------------------------------------------------------------------------------------------------------------
#
MESSAGE:
	@echo " "
	@echo " -----------------------------------------------------------------"
	@echo " -- > Compilation of '$(MAIN)':"
	@echo " -- "
	@echo " -- > Machine --------------------: $(MACHINE)"
	@echo " -- > Operating System -----------: $(OS)"
	@echo " -- > Compiler -------------------: $(FC)"
	@echo " --"
	@echo " ----------------------------------------------------------------"
	@echo " "
#
# ----------------------------------------------------------------------------------------------------------------------------------
# Compilation du programme de test des nouvelles routines PC-SAFT utilisant comme variables indpendantes (rho,T) ou (P,rho) :
demo_pcsaft: $(MAIN_1).f08 foul.o utils_dc.o mod_pcsaft.o 
#
	$(FC) $(OPT) -o $(MAIN_1) $(T) $(MAIN_1).f08 foul.o utils_dc.o mod_pcsaft.o 
#
demo_binary: $(MAIN_2).f08 foul.o utils_dc.o mod_pcsaft.o 
#
	$(FC) $(OPT) -o $(MAIN_2) $(T) $(MAIN_2).f08 foul.o utils_dc.o mod_pcsaft.o
#
# ----------------------------------------------------------------------------------------------------------------------------------
# Modulee :
foul.o: foul.f90
	$(FC) -fimplicit-none -fbacktrace -Ofast -finit-real=zero -finit-integer=0 -fmax-errors=3 -c foul.f90
#
utils_dc.o: utils_dc.f08
	$(FC) $(OPT) -c $(T) utils_dc.f08
#
# ----------------------------------------------------------------------------------------------------------------------------------
# PC-SAFT :
mod_pcsaft.o: mod_pcsaft.f08
	$(FC) $(OPT) -c mod_pcsaft.f08
#
# ----------------------------------------------------------------------------------------------------------------------------------
.PHONY: clean
#
# Faire le ménage :
clean:
	rm -f *.o *.mod *~ *.a $(MAIN) $(HYDRO) $(CLATH) $(TEST_LNPHI) $(SSPEED) $(CPCV) $(THERMODIFF) $(MAPHRES)
