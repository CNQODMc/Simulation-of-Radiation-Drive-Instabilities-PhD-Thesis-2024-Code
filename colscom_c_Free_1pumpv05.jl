#Solves classical model for COoperative/collective SCattering including OptoMechanics
##colscom_c_Free_v03 OAM OAM added to pump and applied to rabi term. Free motion of atoms. 
##Single Atom force terms based on positional derivative of rabi term added to include OAM phase/intensity gradients
##Atoms have free motion in XY plane
##02032022 BTTimeTerm corrected - additional ouput files of constants required for graphing in Python
##07032022 Cooperative Force term updated such that the imaginary part of the whole summation*coherence products is taken
##Introduce option for randomly scattered 2D disk distribution in addition to fixed radius ring by setting option in initialise()
##Introduce option for elliptical distribution by setting additional option in initialise which can be commented out

##v04 calc_perAtom_XYforce function corrected to apply rabi_conj_rderiv=calc_rabi_conj_rderiv(sv)
##instead of previously incorrect calc_rabi_conj_phi_deriv(sv) which means that radial forces are 
##zero for atoms on the intensity for red detuned LG0l pump.

##24052022 Added full betaOAM calculation including collective term for each integration step
##Requires setting for Full BetaOAM (which may have alarge computational load due to inverse N*N matrix calculations per step 
##or Simple BetaOAM - Force terms amd output file term in RK4 function must be changed 
##Added Beta_colscomFv04.dat output file containing na*realBeta,na*imagBeta rows with coulumns for each time iteration

##single pump version with exp(-ikr) in rabi term added z radiation pressure force term and free motion of atoms addded 
##in z direction
##Const_colscomF_1pumpv04

##311022 Updated to apply phase Modulation to simple adiabatic coherence, (betaOAM) claculation 
##PhAmp and PhOmega added as amplitude and angular frequency of phase modulation applied in adiabatic coherence calculation
##As a result cumulative time term (tt) had to be added to rk4,dsv and betaOAM functions
 


using PyPlot
using DelimitedFiles
using LaTeXStrings
using LinearAlgebra

const DistRad=2.0*2.0*pi;     ##optional principal radius of random distribution

#Ellipse predefinition
#const DistRad=1.0*2.0*pi;     ##optional principal radius of random distribution
ax1=DistRad/6.0
ax2=DistRad

naOne=10450
#r=sqrt.(rand(naOne))*DistRad;
r=rand(naOne).^0.55*DistRad;
phi=rand(naOne)*2.0*pi;
x=r.*cos.(phi);
z=r.*sin.(phi);

EllipseSort=Any[]

for j in 1:naOne     
   # phi[j]=rand(1)*2.0*pi;
   # r[j]=sqrt(rand(1))*DistRad;
   # xx=r[j]*cos(phi[j]);
   # yy=r[j]*sin(phi[j]);
   if x[j]^2/ax1^2+z[j]^2/ax2^2 > 1
      push!(EllipseSort,j)
   end
end

deleteat!(x,EllipseSort)
deleteat!(z,EllipseSort)

#ELLIPSE PRE-DEFINITION DONE

#SQUARE or evenly distributed ellipse predefinition
# naSquare=20;
# Zfactor=3;
# SmallAxis=DistRad/Zfactor;
# XX=range(-SmallAxis,SmallAxis,naSquare);
# ZZ=range(-DistRad,DistRad,naSquare*Zfactor);
# #ZZ=range(-SmallAxis,SmallAxis,naSquare);
# #XX=range(-DistRad,DistRad,naSquare*Zfactor);

# x=zeros(naSquare^2*Zfactor);
# z=zeros(naSquare^2*Zfactor);

# for jj in 1:naSquare
#    for kk in 1:naSquare*Zfactor
#       x[Zfactor*naSquare*(jj-1)+(kk)]=(XX[jj]);
#       z[Zfactor*naSquare*(jj-1)+(kk)]=(ZZ[kk]);
#       #x[Zfactor*naSquare*(jj-1)+(kk)]=(XX[kk]);
#       #z[Zfactor*naSquare*(jj-1)+(kk)]=(ZZ[jj]);
#    end
# end

# EllipseSort=Any[]

# for j in 1:naSquare^2*Zfactor    
#    # phi[j]=rand(1)*2.0*pi;
#    # r[j]=sqrt(rand(1))*DistRad;
#    # xx=r[j]*cos(phi[j]);
#    # yy=r[j]*sin(phi[j]);
#    if x[j]^2/(DistRad/Zfactor)^2+z[j]^2/DistRad^2 > 1
#    #if z[j]^2/(DistRad/Zfactor)^2+x[j]^2/DistRad^2 > 1
#       push!(EllipseSort,j)
#    end
# end

#  deleteat!(x,EllipseSort)
#  deleteat!(z,EllipseSort)

#SQUARE or evenly distributed ellipse predefinition DONE


# for j in 1:na     
#    phi[j]=j*1/na*2.0*pi;
#    x[j]=DistRad*cos(phi[j]);
#    y[j]=DistRad*sin(phi[j]);
# end


#x=zeros(0)
#y=zeros(0)


const na=size(x,1)         # Number of atoms  ##for ellipse



#const na=1000;                 # disk or ring                
const ht=2.0;              # Timestep
const tmax=5.5e2;            # Max. time
const omega_r=6.3e-4;        # Dimensionless recoil frequency (Recoil time in terms of Gamma time)
const epsilon=0.22;           # Epsilon (dimensionless separation-cut-off) ***LARGE DIFFERENCE NICOLA/ROLANDO SIMULATIONS***
const Omega=12.0;             # Dimensionless Rabi frequency
const Delta=-100.0;            # Dimensionless detuning
const plots=1000;              # No. of plots to be created
const rho=1.75*2.0*pi;     #Rho where there is a fixed radius for start or whole simulation (and peak intensity of OAM)
const l_mod=0;             #Azimuthal mode number wher OAM is applied
const k0=(2*pi/780.241e-9);   ##k0 value for RD2 resonant transition

##FOR PHASE MODULATION USE PhAmp AS AMPLITUDE AND PhOmega AS FREQUENCY OF MODULATION TO PHI ;; PHI+PhAmP(SIN(PhOmega.t))

const PhAmp=1.0
const PhOmega=0.001

#const DistRad=2.0*2.0*pi;     ##optional radius of random distribution

#const Bwaist=rho/(abs(l_mod)/2)^0.5;      


#const Bwaist=rho*(abs(l_mod)/2)^0.5;
const Bwaist=rho*(2/abs(l_mod))^0.5;      #Pump beamwaist where OAM is applied

const BetaTime=omega_r^0.5*exp(-(rho^2/Bwaist^2))*(2^0.5*rho/Bwaist)^abs(l_mod)*Omega/abs(2*Delta)/rho;

###Create export files of constant names and values for export files for graphing
ConstDataVar = [na,epsilon,Omega,Delta,rho,l_mod,Bwaist,DistRad,PhAmp,PhOmega]
ConstDataList = ["Atoms";"CutOff";"RabiPeakperGamma";"DetuningperGamma";"IntRadius";"AziModeNo";"BeamWaist";"DistRad";"PhModAmp";"PhModOmega"]

ConstFile = open("Const_colscomFv04.dat", "w"); #Open new file 
close(ConstFile);
ConstFile = open("Const_colscomFv04.dat", "a");
writedlm(ConstFile,  ConstDataVar);
close(ConstFile);

ConstName = open("ConstName_colscomFv04.dat", "w"); #Open new file 
close(ConstName)
ConstName = open("ConstName_colscomFv04.dat", "a");
writedlm(ConstName, [ConstDataList]);
close(ConstName);




function open_files()
   posfile = open("positions_colscomFv04.dat", "w"); #Open new file 
   close(posfile);
   posfile = open("positions_colscomFv04.dat", "a");

   momfile = open("momenta_colscomFv04.dat", "w"); #Open new file 
   close(momfile);
   momfile = open("momenta_colscomFv04.dat", "a");

   betafile = open("Beta_colscomFv04.dat", "w"); #Open new file 
   close(betafile);
   betafile = open("Beta_colscomFv04.dat", "a");


   return [posfile ; momfile ; betafile];
end

function initialise()
   # na=naOne
   # #r=sqrt.(rand(na))*2.0*2.0*pi;
   #x=zeros(na);
   y=zeros(na);
   #z=zeros(na);

   phi=zeros(na);
   r=zeros(na);

   
   ##SAMPLE RADIUS ALLOCATION FOR ANNULAR DISTRIBUTION AROUND PEAK INTENSITY RING##
   #r=sqrt.(rand(na).*(0.1*DistRad)).+(DistRad-(0.05*DistRad));
   ##AVERAGE DENSITY RADIUS ALLOCATION FOR DISK OF MAX RADIUS SIZE DISTRAD##
   # r=sqrt.(rand(na))*DistRad;
   # phi=rand(na)*2.0*pi;
   # x=r.*cos.(phi);
   # z=r.*sin.(phi);

#    for j in 1:na     
#       phi[j]=j*1/na*2.0*pi;
#       x[j]=DistRad*cos(phi[j]);
#       #y[j]=DistRad*sin(phi[j]);
#       z[j]=DistRad*sin(phi[j]);

#    end

   # ax1=DistRad
   # ax2=DistRad/3
   # #x=zeros(0)
   # #y=zeros(0)
   # EllipseSort=Any[]
   
   # for j in 1:na     
   #    # phi[j]=rand(1)*2.0*pi;
   #    # r[j]=sqrt(rand(1))*DistRad;
   #    # xx=r[j]*cos(phi[j]);
   #    # yy=r[j]*sin(phi[j]);
   #    if x[j]^2/ax1^2+y[j]^2/ax2^2 > 1
   #       push!(EllipseSort,j)
   #    end
   # end
   
   # deleteat!(x,EllipseSort)
   # deleteat!(y,EllipseSort)
   
   # na=size(x,1)

   # end


   #x = (rand(na)-0.5*ones(na))*5.0*2.0*pi;
   #y = (rand(na)-0.5*ones(na))*5.0*2.0*pi;

   ##TEMPORARY CHANGE TO ADD INITIAL Z MOMENTA
   px=zeros(na);
   py=zeros(na);
   pz=zeros(na);
   sv=zeros(8*na);
   sv[1:na]=x;
   sv[na+1:2*na]=y;
   sv[2*na+1:3*na]=z;
   sv[3*na+1:4*na]=px;
   sv[4*na+1:5*na]=py;
   #sv[5*na+1:6*na].=0.0;
   sv[5*na+1:6*na].=15.0;
   sv[6*na+1:7*na].=0.0;
   sv[7*na+1:8*na].=0.0;
   return sv
end

# function calc_beta()
#    #Calculate beta 
#    beta=im*Omega/2.0/(im*Delta-0.5)*ones(na);
#    return beta
# end

#function calc_betaOAM(sv)
function calc_betaOAM(sv,tt)
   #Calculate beta where rabi term is position dependent adiabatically without cooperative term (same as 2 level Bloch approximation)
   radius,phi=calc_radius_phi(sv);
   betaOAM=zeros(na)im;
   x=sv[1:na];
   z=sv[2*na+1:3*na];

   ##FOR PHASE MODULATION USE PhAmp AS AMPLITUDE AND PhOmega AS FREQUENCY OF MODULATION TO PHI ;; PHI+PhAmP(SIN(PhOmega.t))
   

   for j in 1:na
      #betaOAM[j]=im*Omega/2;#*exp(im*l_mod*phi[j])*2^(abs(l_mod)/2)*(radius[j]/Bwaist)^abs(l_mod)*exp(-(radius[j]^2/Bwaist^2))/(2.0*(im*Delta-0.5));
      #betaOAM[j]=Omega/(2*Delta+im);
      DotProd=dot([0.0;k0],[x[j];z[j]]);
      betaOAM[j]= Omega/2*exp(im*DotProd)/(Delta+im/2);
      ###beta term calculated using phase modulation terms
      ##phase applied directly to rabi term
      #betaOAM[j]= Omega/2*exp(im*DotProd)*exp(im*PhAmp*sin(PhOmega*tt))/(Delta+im/2);

      ##phase applied in time derivative coherence equation to detuning term
      #betaOAM[j]= Omega/2*exp(im*DotProd)/((Delta+PhAmp*sin(PhOmega*tt))+im/2);



   end
   return betaOAM
end

function calc_NewAabsPosnArrays(sv)
   ##Calculate per Integration iteration position difference arrays needed for full coherence calculation with collective term
   NewAabsMatrix=zeros(na,na);
   NewPosnDiffXMatrix=zeros(na,na);
   NewPosnDiffZMatrix=zeros(na,na);
   x=sv[1:na];
   z=sv[2*na+1:3*na];
   radius,phi=calc_radius_phi(sv);

   for jjj in 1:na
      for ii in 1:na

         if ii==jjj
            NewAabsMatrix[jjj,ii]= 0.0;

         elseif ii!=jjj

            NewAabsX = abs(x[jjj]-x[ii]);
            NewAabsZ = abs(z[jjj]-z[ii]);
            #NewAabs = sqrt(NewAabsX^2+NewAabsY^2)
            NewAabs = sqrt(NewAabsX^2+NewAabsZ^2+epsilon^2);
            NewAabsMatrix[jjj,ii]=NewAabs;
            NewPosnDiffX = x[jjj]-x[ii];
            NewPosnDiffZ = z[jjj]-z[ii];
            NewPosnDiffXMatrix[jjj,ii]=NewPosnDiffX;
            NewPosnDiffZMatrix[jjj,ii]=NewPosnDiffZ;
         end
      end
   end
   return NewAabsMatrix#,NewPosnDiffXMatrix,NewPosnDiffYMatrix
end

function Coherence_matrix_create(Matrix,na)
   ##function to create matrix component representing coupled equations for time derivative of coherence terms
   DipoleMatrix=zeros(na,na)im
   for lll in 1:na
      for kk in 1:na
         if kk==lll
            DipoleMatrix[lll,kk]= im*Delta-0.5;
         else
            DipoleMatrix[lll,kk]= -0.5*(exp(im*Matrix[lll,kk])/(im*Matrix[lll,kk]));
         end
      end
   end
   return DipoleMatrix
end

function calc_FullBetaOAM(sv,tt)    ##CAN INSERT POSN EXP TERM FOR PUMP INTERATION IN COHERENCE OR NOT##
   #Calculate beta where rabi term is position dependent adiabatically with full cooperative term (per quantum or classical derivation 
   ## of the time derivative of coherence)
   #radius,phi=calc_radius_phi(sv);
   x=sv[1:na];
   z=sv[2*na+1:3*na];   
   FullBetaOAM=zeros(na)im;
   CollectiveVectorResultant=zeros(na)im;
   NewAbMatrix=calc_NewAabsPosnArrays(sv);
   DipoleMatrix=Coherence_matrix_create(NewAbMatrix,na);
   for hhh in 1:na
      #CollectiveVectorResultant[hhh]=im*Omega/2#*exp(im*l_mod*phi[hhh])*2^(abs(l_mod)/2)*(radius[hhh]/Bwaist)^abs(l_mod)*exp(-(radius[hhh]^2/Bwaist^2));
      DotProd=dot([0.0;k0],[x[hhh];z[hhh]]);
      #CollectiveVectorResultant[hhh]=im*Omega/2*exp(im*DotProd);
      #CollectiveVectorResultant[hhh]=im*Omega/2;

      #PhaseModulation applied
      CollectiveVectorResultant[hhh]=im*Omega/2*exp(im*DotProd)*exp(im*PhAmp*sin(PhOmega*tt));
   end
   InvDipoleMatrix=inv(DipoleMatrix);
   FullBetaOAM=InvDipoleMatrix*CollectiveVectorResultant;

   return FullBetaOAM

end



                            



function calc_radius_phi(sv)
   #Calculate per atom radius
   x=sv[1:na];
   y=sv[na+1:2*na];
   z=sv[2*na+1:3*na];
   radius=zeros(na);
   phi=zeros(na);
   for j in 1:na
      radius[j]=(x[j]^2+z[j]^2)^0.5;
      phi[j]=atan(z[j],x[j]);
   end
   #beta=im*Omega/2.0/(im*Delta-0.5)*ones(na);
   return radius,phi
end

# function calc_Bwaist(sv)
#    #Calculate Beamwaist of source Gaussian so that peak radial intensity is on atom ring
#    radius,phi=calc_radius_phi(sv);
#    Bwaist=radius./(abs(l_mod)/2).^0.5;
#    #beta=im*Omega/2.0/(im*Delta-0.5)*ones(na);
#    return Bwaist
# end

# function calc_rabi_rho(sv)
#    #Calculate OAM rabi per atom position
#    radius,phi=calc_radius_phi(sv);
#    Bwaist=calc_Bwaist(sv);
#    rabi_rho=zeros(na);
#    for j in 1:na
#       rabi_rho[j]=Omega*exp(-(radius[j]^2/Bwaist[j]^2))*(2^0.5/Bwaist[j])^abs(l_mod);
#    end
#    return rabi_rho
# end

function calc_rabi_conj_rderiv(sv)
   #Calculate per atom x radial derivative of position dependent rabi
   radius,phi=calc_radius_phi(sv)
   rabi_conj_rderiv=zeros(na)im
   for j in 1:na
      rabi_conj_rderiv[j]=-(Omega*exp(-im*l_mod*phi[j])*2^(abs(l_mod)/2)*(radius[j]/Bwaist)^abs(l_mod)*(2*radius[j]^2-abs(l_mod)*Bwaist^2)*exp(-(radius[j]^2/Bwaist^2)))/(Bwaist^2*radius[j]);
   end
   return rabi_conj_rderiv
end

function calc_rabi_conj_phi_deriv(sv)
   #Calculate per atom x azimuthal derivative of position dependent rabi convert result in radians with phi unit vector to dimensionless kr result for XY conversion
   radius,phi=calc_radius_phi(sv)
   rabi_conj_phi_deriv=zeros(na)im
   for j in 1:na
      #rabi_conj_phi_deriv[j]=-Omega*im*l_mod*exp(-im*l_mod*phi[j])*2^(abs(l_mod)/2)*(radius[j]/Bwaist)^abs(l_mod)*exp(-(radius[j]^2/Bwaist^2));#*radius[j];
      rabi_conj_phi_deriv[j]=-Omega*im*l_mod*exp(-im*l_mod*phi[j])*2^(abs(l_mod)/2)*(radius[j]/Bwaist)^abs(l_mod)*exp(-(radius[j]^2/Bwaist^2));
   end
   return rabi_conj_phi_deriv
end

function calc_perAtom_XYZforce(sv)  ##CONSIDER DIMENSIONLESS MOPMENTUM IN TERMS OF PER PHOTON AND DISTANCE/TIME
   #Calculate per atom X,Y force terms based on radial and azimuthal derivatives of position dependent rabi term
   radius,phi=calc_radius_phi(sv)
   #betaOAM=calc_FullBetaOAM(sv)
   #betaOAM=calc_betaOAM(sv)
   #rabi_conj_rderiv=calc_rabi_conj_rderiv(sv)
   #rabi_conj_phi_deriv=calc_rabi_conj_phi_deriv(sv)
   f_1x=zeros(na);
   f_1y=zeros(na);
   f_1z=zeros(na);
   for j in 1:na
      #f_1x[j]=-2*real(betaOAM[j]*(rabi_conj_rderiv[j]*cos(phi[j])))#+rabi_conj_phi_deriv[j]*-sin(phi[j])));
      #f_1y[j]=-2*real(betaOAM[j]*(rabi_conj_rderiv[j]*sin(phi[j])))#+rabi_conj_phi_deriv[j]*cos(phi[j])));
      #f_1z[j]=0.2/(Delta^2+0.25)*(Omega/2)^2
      f_1z[j]=1.0/(Delta^2+0.25)*(Omega/2)^2
   end
   return f_1x,f_1z#,f_1y
end

function calc_q(x_j,z_j,x,z)
   #Calculate q_jm
   q=((x_j*ones(na)-x).^2+(z_j*ones(na)-z).^2+epsilon^2*ones(na)).^0.5;
   return q;
end

#print(calc_radius(sv))

#function calc_coop_force(sv)
function calc_coop_force(sv,tt)
   #Calculate cooperative force on atoms
   x=sv[1:na];
   y=sv[na+1:2*na];
   z=sv[2*na+1:3*na];
   f_x=zeros(na);
   f_y=zeros(na);
   f_z=zeros(na);
   betaOAM=calc_FullBetaOAM(sv,tt);
   #betaOAM=calc_betaOAM(sv,tt)
   #betaOAM=calc_betaOAM(sv,tt);
   for j in 1:na
       q_jm=calc_q(x[j],z[j],x,z);
       #f_x[j]=-sum(conj(betaOAM[j]).*betaOAM.*(sin.(q_jm)./q_jm+cos.(q_jm)./q_jm.^2).*(x[j]*ones(na)-x)./q_jm);
       #f_y[j]=-sum(conj(betaOAM[j]).*betaOAM.*(sin.(q_jm)./q_jm+cos.(q_jm)./q_jm.^2).*(y[j]*ones(na)-y)./q_jm);
       f_x[j]=-imag(sum(conj(betaOAM[j]).*betaOAM.*(exp.(q_jm*im)./q_jm+exp.(q_jm*im).*im./q_jm.^2).*(x[j]*ones(na)-x)./q_jm));#.*exp.(q_jm*im);
       f_z[j]=-imag(sum(conj(betaOAM[j]).*betaOAM.*(exp.(q_jm*im)./q_jm+exp.(q_jm*im).*im./q_jm.^2).*(z[j]*ones(na)-z)./q_jm));#.*exp.(q_jm*im);
   end
   return f_x,f_z
end


#function dsv(sv)
function dsv(sv,tt)
   # Derivatives of state vector
    dsv=zeros(8*na);

    x=sv[1:na];
    y=sv[na+1:2*na];
    z=sv[2*na+1:3*na];
    px=sv[3*na+1:4*na];
    py=sv[4*na+1:5*na];
    pz=sv[5*na+1:6*na];
    #b1=sv[4*na+1:5*na];
    #b2=sv[5*na+1:6*na];

    #betaOAM=calc_FullBetaOAM(sv);
    #betaOAM=calc_betaOAM(sv,tt);
    #M,phi = order(theta);

    #beta=calc_beta();
    #fco_x,fco_y=calc_coop_force(sv).+calc_perAtom_XYforce(sv)
    #fcoop_x,fcoop_z=calc_coop_force(sv);
    fcoop_x,fcoop_z=calc_coop_force(sv,tt);
    fsing_x,fsing_z=calc_perAtom_XYZforce(sv);

    ##FOR TESTING APPLY ONLY SINGLE OR COOPERATIVE FORCES
    fco_x=fcoop_x+fsing_x;
    #fco_x=fsing_x;
    fco_y=zeros(na);
    #fco_y=fcoop_y+fsing_y;
    fco_z=fcoop_z+fsing_z;
    #fco_z=fsing_z;

   #  fco_x=fsing_x;   ##ENABLE IF CHECKING EFFECT OF PER ATOM FORCE
   #  fco_y=fsing_y;   ##ENABLE IF CHECKING EFFECT OF PER ATOM FORCE

    #fco_x=fcoop_x;   ##ENABLE IF CHECKING EFFECT OF COLLECTIVE FORCE
    #fco_y=fcoop_y;   ##ENABLE IF CHECKING EFFECT OF COLLECTIVE FORCE


    dsv[1:na]=2.0*omega_r*px;
    dsv[na+1:2*na]=2.0*omega_r*py;
    dsv[2*na+1:3*na]=2.0*omega_r*pz;
    dsv[3*na+1:4*na]=fco_x;
    dsv[4*na+1:5*na]=fco_y;
    dsv[5*na+1:6*na]=fco_z;
    #dsv[6*na+1:7*na]=real(betaOAM);
    #dsv[7*na+1:8*na]=imag(betaOAM);

    
    return dsv
end

function rk4(ht,tt,svt)
   
   ##add coherence terms to time iterated state vector
   #NewSV=dsv(svt);
   NewSV=dsv(svt,tt);
   RealBeta=NewSV[6*na+1:7*na];
   ImagBeta=NewSV[7*na+1:8*na];
   
   
   #4th order Runge-Kutta method for position and momenta (sv 1:8*na)

    svk1=ht*dsv(svt,tt);
    svk2=ht*dsv(svt+svk1/2.0,tt);
    svk3=ht*dsv(svt+svk2/2.0,tt);
    svk4=ht*dsv(svt+svk3,tt);

    svt=svt+1.0/6.0*(svk1 + 2.0*svk2 + 2.0*svk3  + svk4);

    ##restore calculated position dependent coherence real and imag values - which dont use RK4 to state vector
    svt[6*na+1:7*na]=RealBeta;
    svt[7*na+1:8*na]=ImagBeta;

    tt=tt+ht;
    return tt, svt
end

# function plot_xy(t,sv)
#    #Plot positions
#    @show t
#    x=sv[1]-sv[2];
#    y=sv[na+1]-sv[na+2];
#    N1N2=abs(sqrt(x^2+y^2));
#    Rabit=calc_rabi_rho(sv);
#    #plt=scatter(t,N1N2);
#    plt=scatter(t,Rabit[1]);
#    X_txt=string(L"\tau");
#    xlabel(X_txt);
#    ylabel("N1 N2 Posn Difference")
#    #title_txt=string(L"\tau =",t);
#    title("Optical Binding Simulation 3 Atoms");
#    savefig("coscom_c_PosnDiff.png")
#    draw()

#    return
# end


function plot_xz(t,sv)
   #Plot positions
   @show t
   x=sv[1:na];
   #y=sv[na+1:2*na];
   z=sv[2*na+1:3*na];
   plt=scatter(x,z);
   xlabel("X");
   #ylabel("Y");
   ylabel("Z");
   #title_txt=(L"Optical Binding Simulation $(na) Free Atoms, $(Delta) Detuning, $(Omega) Peak Rabi, \tau =" ,t);
   #title(title_txt);
   title_txt=("Optical Binding Simulation $(na) atoms, $(Delta) detuning, $(Omega) peak rabi,\n"); #Atoms, $Delta Detuning, $Omega Peak Rabi, OAM mode l=$l_mod, \tau =" ,t)
   title_txt2=string(L" \tau (recoil) =",t*omega_r);
   #title_txt3=string(", DistRadius =",DistRadius);
   #title_txt4=string("PeakI Radius =",)


   #(L"Optical Binding Simulation $(na) Atoms, $(Delta) Detuning, $(Omega) Peak Rabi, OAM mode l=$(l_mod), \tau =" ,t)
   title(title_txt*title_txt2);
   savefig("coscomv04_c_OAM_Free_1pump.png")
   draw()

   return
end


function write_dat(t,sv,file_ids)
   #Write data to files
   x=sv[1:na];
   y=sv[na+1:2*na];
   z=sv[2*na+1:3*na];

   px=sv[3*na+1:4*na];
   py=sv[4*na+1:5*na];
   pz=sv[5*na+1:6*na];
   betaR=sv[6*na+1:7*na];
   betaIm=sv[7*na+1:8*na];


   dataout=zeros(2*na+1);      # Write position data to file
   dataout[1]=t;
   dataout[2:na+1]=x;
   dataout[na+2:2*na+1]=z;
   posfile=file_ids[1];
   writedlm(posfile,  dataout', ',');

   dataout=zeros(2*na+1);      # Write momentum data to file
   dataout[1]=t;
   dataout[2:na+1]=px;
   dataout[na+2:2*na+1]=pz;
   momfile=file_ids[2];
   writedlm(momfile,  dataout', ',');

   dataout=zeros(2*na+1);      # Write beta data to file
   dataout[1]=t;
   dataout[2:na+1]=betaR;
   dataout[na+2:2*na+1]=betaIm;
   betafile=file_ids[3];
   writedlm(betafile,  dataout', ',');

   return
end

function integrate(tmax,sv,file_ids)
   #Integrate state vector in time
   tvec=LinRange(0.0,tmax,plots);
   t=0.0;
   frame=1;     
   plot_xz(t,sv);
   #plot_N1N2PosnDiff(t,sv);
   write_dat(t,sv, file_ids)
   nextplot=tvec[frame+1];
   while (t<tmax)
        t,sv=rk4(ht,t,sv);
        if t >= nextplot
            plot_xz(t,sv)
            #plot_N1N2PosnDiff(t,sv)
            write_dat(t,sv, file_ids)
            frame=frame+1;
            if frame<plots
               nextplot=tvec[frame+1];
            end
        end
    end
    return 
end

function close_files(file_ids)
   #Close data files
   for file_id in file_ids
      close(file_id);
   end
end
 
# Main loop
print("Starting...")
file_ids=open_files();
sv=initialise();       #Initialise state vector
#print(sv);
ion();                 #Interactive mode for real-time plotting
figure();              #Create figure
integrate(tmax,sv,file_ids);
close_files(file_ids);
print("Finished")