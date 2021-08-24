function sol =  GetPlumeDimensional(Ti, Si, T0, T1, S0, S1, Z0, rho0, zgl, L, ci,lt,g,Cd,...
                                  c, bs, bt, E0, St, lambda1, lambda2, lambda3, tau,...
                                  Xb, Zb, dZb)
%Solve the dimensional plume equations with pycnocline with ice draft and
%ice draft slope inputted as arrays. 
%addpath('../Auxillary_functions/') %need for interpolate function
%% Anonymous functions of Z dependent variables
Ta = @(Z) (T0 + T1)/2 - (T0 - T1)/2 *tanh((Z - Z0)/lt); %ambient temperature
Sa = @(Z) (S0 + S1)/2 - (S0 - S1)/2 *tanh((Z - Z0)/lt); %ambient salinity

dTa_dZ = @(Z)  - (T0 - T1)/2/lt *sech((Z - Z0)/lt)^2; %ambient temperature gradient
dSa_dZ = @(Z)  - (S0 - S1)/2/lt *sech((Z - Z0)/lt)^2; %ambient salinity gradient
drhoa_dZ = @(Z) rho0 *(bs*dSa_dZ(Z) - bt*dTa_dZ(Z));

Tf = @(S,Z) lambda1*S + lambda2 + lambda3*Z; %liquidus equation
Delta_T_a = @(Z) Ta(Z) - Tf(Sa(Z),zgl + Z);
Delta_T_ief = @(Z) -(L + ci*(Tf(Si,Z) - Ti))/c; 
T_ief = @(Z) Tf(0,Z) + Delta_T_ief(Z); 
Delta_rho_ief = @(Z) rho0*(bs*(Sa(Z) - Si) - bt*(Ta(Z) - T_ief(Z)));

%% Initial conditions from similarity solution
xcross = 100; %specify initial condition location (dimensional)
%compute slope at 0 using 1 sided fd
dZb0 = (Zb(2) - Zb(1))/(Xb(2) - Xb(1));
E      = E0*dZb0; %notation from Hewitt 2020
Ltilde = L + ci*(lambda1*S0 + lambda2 - Ti); %notation from Hewitt 2020

Delta_T_0   = E/(E+St) *tau; %thermal driving initial condition
Delta_rho_0 = St/E * c*Delta_T_0/Ltilde *Delta_rho_ief(0);
U_0         = (E*Delta_rho_0*dZb0*g/rho0 / (2*E+(3/2)*Cd))^(1/2)*xcross^(1/2);
D_0         = 2/3*E*xcross;
Y0 = [D_0, U_0, Delta_rho_0, Delta_T_0];

%% solve ODEs
rhs = @(X,Y) forcing(X,Y,Xb,Zb, dZb, E0, g, rho0, lambda3, St, L, c, Cd,Delta_rho_ief, drhoa_dZ, Delta_T_a, Delta_T_ief);
M = @(X,Y) massmat(X,Y);
options = odeset('Mass',M, 'Events',@plumeEventsFcn,...
                 'RelTol', 1e-7);
sol  = ode45(rhs,[0,max(Xb)],Y0,options);


%% functions: mass matrix, forcing and events function
function M = massmat(X,Y)
%Return the mass matrix associated with one-d plume dynamics.
% Y = [d,u,r,q]

D = Y(1); 
U = Y(2);
delta_rho = Y(3); %buoyancy deficit
delta_T = Y(4); %dimensionless thermal driving

M = [U, D, 0, 0;
     U^2, 2*D*U, 0, 0;
     U*delta_rho, D*delta_rho, D*U, 0;
     U*delta_T, D*delta_T, 0, D*U];

end

function f = forcing(X,Y,Xb, Zb, dZb, E0, g, rho0, lambda3, St, L, c, Cd,Delta_rho_ief, drhoa_dZ, Delta_T_a, Delta_T_ief)
%return the forcing array
D = Y(1); 
U = Y(2);
delta_rho = Y(3); %buoyancy deficit
delta_T = Y(4); %dimensionless thermal driving
M0 = St/(L/c); %melt rate constant

%use interpolate to compute depth and slope
Zb_here = interpolate(Xb,Zb,X);
dZb_here = interpolate(Xb,dZb,X);

f = [E0*U*dZb_here + M0*U*delta_T;
    g*D*delta_rho*dZb_here/rho0 - Cd*U^2;
    M0*U*delta_T*Delta_rho_ief(Zb_here) + dZb_here*drhoa_dZ(Zb_here)*D*U;
    E0*U*dZb_here*Delta_T_a(Zb_here) + M0*U*delta_T*Delta_T_ief(Zb_here) - lambda3*D*U*dZb_here];

end

function [position,isterminal,direction] = plumeEventsFcn(X,Y)
position = Y(2)-1e-4; % stop integrating when speed goes through zero
isterminal = 1;  % Halt integration 
direction = -1;   %Only descrease through
end

end