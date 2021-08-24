function sol = GetPlume(eps1,eps2, eps3,eps4,delta, Pb, Pt, lambda, x0,zbF,dzbF,Xmax)
%return the solution structure associated with the plume equations with a
%non constant geometry specified by zb = zb(xb). We interpolate to the grid
%with a simple linear interpolation to nearest neighbour 


% zbF   :   Anonymous function
%           Ice shelf draft 
%
% dzbF  :   Anonymous function
%           Slope of ice shelf draft
%

%initial conditions from analytic solution (based on ignoring (a) depth
%dependence terms, (b) melting in mass equation and (c) fat plume term
xcross = 1e-4; %mu in the paper
d0 = 2/3 * xcross;
u0 = sqrt(xcross) * sqrt(lambda/(1 + eps2) * 1/(2*eps1 + 3/2));
drho0 = lambda/(1 + eps2);
dt0 = 1/(1 + eps2);
Y0 =  [d0;u0;drho0;dt0];

%solve odes:
rhs = @(x,Y)  forcing(x,Y,delta,eps3,eps4, Pb,Pt,lambda,x0, zbF,dzbF);
M = @(x,Y) massmat(x,Y, eps1,eps2);
 options = odeset('Mass',M, 'Events',@plumeEventsFcn,...
                  'RelTol', 1e-8);
%options = odeset('Mass',M, 'Events',@plumeEventsFcn);
sol  = ode15s(rhs,[0,Xmax],Y0,options);

function M = massmat(x,Y, eps1,eps2)
%Return the mass matrix associated with one-d plume dynamics.
% Y = [d,u,r,q]
D = Y(1); 
U = Y(2);
delta_rho = Y(3); %dimensionless buoyancy deficit
delta_T = Y(4); %dimensionless thermal driving

M = [U, D, 0, 0;
     eps1*U^2, eps1*2*D*U, 0, 0;
     U*delta_rho, D*delta_rho, D*U, 0;
     eps2*U*delta_T, eps2*D*delta_T, 0, eps2*D*U];
%  
end

function f = forcing(x,Y,delta,eps3,eps4, Pb,Pt, lambda,x0, zbF,dzbF)
%return the forcing array
D = Y(1);
U = Y(2);
delta_rho = Y(3);
delta_T = Y(4);

f = [U*dzbF(x) + eps3*U*delta_T;
    D*delta_rho*dzbF(x) - U^2;
    -Pb/delta * sech((x-x0)/delta)^2 * D*U*dzbF(x)+ U*delta_T*(lambda - eps4*tanh((x-x0)/delta));
    -Pt*(1 + tanh((x-x0)/delta))*dzbF(x)*U + (1-zbF(x))*dzbF(x)*U - D*U*dzbF(x) - U*delta_T];

end

function [position,isterminal,direction] = plumeEventsFcn(x,Y)
position = Y(2)-1e-4; % stop integrating when speed goes through zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction
end
end
