# colscom
Cooperative Scattering Involving Optomechanics

## colscom_c_Free_v**
Classical model of cooperative scattering involving optpmechanics in 2D.

##v.02 Adds OAM to pump field - Adds single atom pump interations per documented phase and intensity gradients

##v0.3 Corrects application of single atom field interaction terms such that conjugate correctly applied either to rabi or coherence term

##v0.3 Also adds term in the RK4 calculation to keep iterated positions on fixed radius and keep iterated momenta in tangential directions

##v0.3 Correction to BetaTime Term and output file of constants created for graphing

##v0.3 07032022 Correction to cooperative force calculation to take the imaginary part of the product of coherences and spatial derivative of Greens function ##**This requires an update to the associated equations derivation document

##v0.3 01042022 No significant updates to code related python graphing file added. It should be noted that the cut-off term epsilon is 
##in dimensionless k-space in CARL OAM papers it can be referred to as a radians equivalent with a radial dependency.

##v04 16052022 calc_perAtom_XYforce function corrected to apply rabi_conj_rderiv=calc_rabi_conj_rderiv(sv)
##instead of previously incorrect calc_rabi_conj_phi_deriv(sv) which means that radial forces are 
##zero for atoms on the intensity for red detuned LG0l pump.

##24052022 Added full betaOAM calculation including collective term for each integration step
##Requires setting for Full BetaOAM (which may have alarge computational load due to inverse N*N matrix calculations per step or Simple BetaOAM 
##BetaOAM needs to be set to equal to calc_FullBetaOAM or calc_betaOAM in both calc_perAtom_XYforce(sv) and calc_coop_force(sv) functions
##26052022 State vector extended to include real and imag parts of coherence calculation for each time step

##31052022Added Beta_colscomFv04.dat output file containing na*realBeta,na*imagBeta rows with coulumns for each time iteration

##v05single pump version with exp(-ikr) in rabi term added z radiation pressure force term and free motion of atoms addded 
##in z direction
##Const_colscomF_1pumpv04

##v05 311022 Updated to apply phase Modulation to simple adiabatic coherence, (betaOAM) claculation 
##PhAmp and PhOmega added as amplitude and angular frequency of phase modulation applied in adiabatic coherence calculation
##As a result cumulative time term (tt) had to be added to rk4,dsv and betaOAM functions

##Beta e^ikr term calculated purely where k vector is 1 in z direction = e^iz (z is in units of k space)

##Modulation terms added and applied directly to rabi term calculations for AMod or PhaseMod Full coherence term used almost exclusively 
##now although simplified coherence term very good approximation in dilute Nicola/Angel regimes

