%Store dimensional parameters used in solution. (this would be better as a
%structure to avoid global variables being called)
E0      = 1*1e-2;   %entrainment coefficient
Cd      = 1*1e-3;   %drag coefficient
St      = 5.9*1e-4;   %effective stanton number 
lambda1 = -5.73*1e-2;
lambda2 = 8.32*1e-2;
lambda3 = 7.61*1e-4;  %Freezing point depth coefficient
g       = 9.81; 
rho0    = 1e3;        %reference density
alpha   = 3*1e-3;   %slope of the base at Z =0

lt      = 50;         %thickness of thermocline
bs      = 7.86*1e-4;  %Haline contraction coefficient
bt      = 3.87*1e-5;  %Thermal expansion coefficient
L       = 3.35*1e5;   %Latent heat of fusion for ice
c       = 3.974*1e3;  %specific heat capicity (water)
ci      = 2*1e3;      %specific heat capicity (ice)
Ti      = 0;          %ice temperature
Si      = 0;          %ice salinity

T0      = 0.5;          %lower layer temperature (C)
T1      = -1.5;         %upper layer temperature (C)
S0      = 34.6;         %lower layer salinity (psu)
S1      = 34.0;         %upper layer salinity (psu)
Sd      = S0 - S1;      %difference in salinity across layer
Td      = T0 - T1;      %difference in temperature across layer
zgl     = -1500;         %depth of the grounding line
tau     = T0 - (lambda1*S0 + lambda2 + lambda3*zgl);%T_{a,0} - T_{f,0}
l0      = tau/lambda3;%lengthscale of freezing pt dependence

secs_per_yr = 365*24*60*60;