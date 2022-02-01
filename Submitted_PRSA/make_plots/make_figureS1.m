%Make figure supplemenrary 1: comparison between numerical solutions of the rescaled
%equations and leading order form of rescaled equations. This code produces
%two plots: the first presents the solutions in (u, p) space and the second
%has two plots of the solution in (chi,u) and (chi, p) space.

%% Preliminaries
addpath('Auxillary_functions');
figpref(4);

%uncomment for R2020a
set(0, 'defaultaxesfontsize', 12);
set(0, 'defaultaxeslinewidth', 1.5);
set(0, 'defaultlinelinewidth', 1.5);
set(0, 'defaultpatchlinewidth', 0.7); 
%% parameters
run parameters %get dimensional parameters (brings them all into global scope)

%dimensionless parameters
eps1 = E0*alpha/Cd;
eps2 = E0*alpha/St;
eps3 = tau/(L/c);
eps4 = (S0 - S1)/2/S0;
delta = lt/l0;
Pb = (L/c)/tau * (S0 - S1) /2 / S0 *( 1- bt*(T0 - T1)/bs / (S0 - S1));
Pt = (T0 - T1) / 2 / tau; %or Pt = (T0 - T1 + lambda1*(S0 - S1) / 2 / tau but makes little difference
kappa = (S0 + S1)/2 / S0 - bt*(L/c)/ bs / S0;

%order 1 quantities relevant for the transition region
k2      = eps2/eps1;
k3      = eps3/eps1;
Qc      = 0.4;      %artifically provide some Qc (in reality determined as part of solution)
Xc      = 0.4;      %as above
zbF     = @(X)  X;  %linear draft
dzbF    = @(X) 1 + 0*X; 
dzb_xc  = dzbF(Xc);  
zb_xc   = zbF(Xc);

%variable scalings
chi_scaling = kappa^(-1/4)*Qc^(1/2)*dzb_xc^(-1/2);
u_scaling = kappa^(1/4)*(dzb_xc)^(1/2)*Qc^(1/2);
p_scaling = kappa^(3/4)*(dzb_xc)^(1/2)*Qc^(1/2);

% Far field conditions
chi_ff      = -10*chi_scaling; %choose matching point (|X_ff| >> 1, X_ff < 0) Note: length scale is eps1^(3/4) \approx 0.05, so in reality we probably only need 1 length scale (i.e. no need to ramp X_ff up too massively)
U_ff        = (-chi_ff)^(1/3)*dzb_xc^(2/3)*Qc^(1/3)*kappa^(1/3) ;
D_ff        = (-chi_ff)^(-1/3)*dzb_xc^(-2/3)*Qc^(2/3)*kappa^(-1/3);
delta_P_ff  = -chi_ff*dzb_xc*kappa;
delta_T_ff  = -(-chi_ff)^(1/3)*Qc^(2/3)*dzb_xc^(1/3) * kappa^(-1/3);
Y0          = [D_ff,U_ff, delta_P_ff,delta_T_ff];

%% Plot 1:
figure(1); clf; 
ax1 = subplot(1,2,1); hold on
colmap = parula(3); %colormap for full solutions
%horizontal line where negatively buoyant
plot([0,2.5], [0,0], 'k', 'linewidth', 1.5, 'HandleVisibility', 'off')

%Full equations (regular parameters):
M1   = @(chi,Y) mass1(chi,Y, k2);
rhs1 = @(chi,Y) forcing1(chi,Y,eps1,eps4, kappa, Pt,k3,zbF, dzbF,Xc);
options = odeset('Mass',M1, 'RelTol', 1e-6,'AbsTol', 1e-6);
[chi1, Y1]  = ode15s(rhs1,[chi_ff,abs(chi_ff)],Y0,options);
U1 = Y1(:,2);
delta_P1 = Y1(:,3);

%plot nullcline as first layer
plot(U1/u_scaling,(U1/u_scaling).^3,'color',  [170,170,170]/255, 'linewidth',3)

%result with typical values
plot(U1/u_scaling, delta_P1/p_scaling,'color',  colmap(1,:), 'linewidth',3)

%repeat with small parameters scaled down:
eps1 = eps1/10; eps4 = eps4/10;
M1   = @(chi,Y) mass1(chi,Y, k2);
rhs1 = @(chi,Y) forcing1(chi,Y,eps1,eps4, kappa, Pt,k3,zbF, dzbF,Xc);
options = odeset('Mass',M1, 'RelTol', 1e-6, 'AbsTol', 1e-6);
[chi2, Y2]  = ode15s(rhs1,[chi_ff,abs(chi_ff)],Y0,options);
U2 = Y2(:,2);
delta_P2 = Y2(:,3);
plot(U2/u_scaling, delta_P2/p_scaling,'color',colmap(2,:), 'linewidth',2)


%numerical solution rescaled LO equations:
chi_ff_reduced = chi_ff/chi_scaling; %reduced far field conditions
u_tilde_ff = (-chi_ff_reduced)^(1/3);
p_tilde_ff = -chi_ff_reduced;
Y0_tilde = [u_tilde_ff, p_tilde_ff];
 
M2   = @(chi,Y) mass2(chi,Y, k2, kappa);
rhs2 = @(chi,Y) forcing2(chi,Y);
options = odeset('Mass',M2, 'RelTol', 1e-6,'AbsTol', 1e-6);
[chi_reduced, Y_reduced] = ode15s(rhs2,[chi_ff_reduced,abs(chi_ff_reduced)],Y0_tilde,options);
U_tilde = Y_reduced(:,1);
delta_P_tilde = Y_reduced(:,2);
plot(U_tilde, delta_P_tilde, 'k--','linewidth', 2)


%tidy plot
box on
ylim([-2, delta_P_ff/p_scaling])
ylabel('$\tilde{\Delta \varrho}$', 'interpreter','latex')
xlabel('$\tilde{\mathcal{U}}$', 'interpreter','latex')

% %% Plot 2: v_tilde and p_tilde versus chi
% figure(2); clf; 
% subplot(2,1,1); hold on
% plot(chi1/chi_scaling, U1/u_scaling,'color',  colmap(1,:), 'linewidth',3)
% plot(chi2/chi_scaling, U2/u_scaling,'color',  colmap(2,:), 'linewidth',3)
% plot(chi_reduced, U_tilde,'k--', 'linewidth',2.5)
% box on
% xlabel('$\tilde{\chi}$', 'interpreter','latex')
% ylabel('$\tilde{u}$', 'interpreter','latex')
% ylim([0,2.5])
% xlim([-10,5])
% xticks([-10:5:5])

ax2 = subplot(1,2,2); hold on
plot([-10,5], [0,0], 'k', 'linewidth', 1.5)
plot(chi1/chi_scaling, delta_P1/p_scaling,'color',  colmap(1,:), 'linewidth',2)
plot(chi2/chi_scaling, delta_P2/p_scaling,'color',  colmap(2,:), 'linewidth',2)
plot(chi_reduced, delta_P_tilde,'k--', 'linewidth',2.5)
box on
xlabel('$\tilde{\chi}$', 'interpreter','latex')
ylabel('$\tilde{\Delta \varrho}$', 'interpreter','latex')
ylim([-2,10])
xlim([-10,5])
xticks(-10:5:5)
fig = gcf;fig.Position(3:4) = [900,325];

txta= text(ax1,-0.6,10, '(a)','Interpreter','latex', 'FontSize', ax1.XLabel.FontSize);
txtb= text(ax2,-13.5,10, '(b)','Interpreter','latex', 'FontSize', ax1.XLabel.FontSize);

legend(ax1, {'$\tilde{\Delta \rho} = \tilde{{u}}^3$', '$\epsilon_1 = 3\times 10^{-2}$', '$\epsilon_1 = 3\times 10^{-3}$', 'Reduced equations', ''},...
    'interpreter', 'latex', 'location', 'northwest')
% saveas(gcf,'plots/figure5.png') %figure 5 in original submission
% saveas(gcf,'plots/supp_fig1.png')
%% functions
function M = mass1(chi,Y, k2)
%return the mass matrix for the solution method 1 (solving full equations
%with rescaled length and tanh == 1, sech == 0)
D = Y(1); 
U = Y(2);
delta_P = Y(3); %dimensionless buoyancy deficit
delta_T = Y(4); %dimensionless thermal driving

M = [U, D, 0, 0;
     U^2, 2*D*U, 0, 0;
     U*delta_P, D*delta_P, D*U, 0;
     k2*U*delta_T, k2*D*delta_T, 0, k2*D*U];

end

function f = forcing1(chi,Y,eps1,eps4, kappa, Pt,k3,zbF, dzbF,Xc)
%return the forcing for the solution method 1 (solving full equation
%with rescaled length and tanh == 1, sech == 0)
D = Y(1);
U = Y(2);
delta_P = Y(3);
delta_T = Y(4);

f = [eps1 * U * dzbF(Xc+eps1^(3/4)*chi) + eps1^(3/4) * k3 * U * delta_T;
    D*delta_P*dzbF(Xc+eps1^(3/4)*chi) - U^2;
    U*delta_T*(kappa - eps4);
    eps1^(1/4) * (1 - 2*Pt  - zbF(Xc+eps1^(3/4)*chi)) * dzbF(Xc+eps1^(3/4)*chi) *U + ...
        -U * delta_T - D*U*dzbF(Xc+eps1^(3/4)*chi)];

end

function M = mass2(chi,Y, k2, lambda)
%return the mass matrix for the reduced leading order equations (Note that
%this matrix is diagonal so you could easily invert it by hand!)
v       = Y(1); %rescaled speed
q       = Y(2); %rescaled buoyancy
M = [v,   0;
    0,    k2*lambda];
end

function f = forcing2(chi,Y)
%return the forcing for the reduced leading order equations
v       = Y(1); %rescaled speed
q       = Y(2); %rescaled buoyancy
f = zeros(2,1); %ode15s needs a column vector
f(1)    = q - v^3;
f(2)    = -v*(q + chi);
end

