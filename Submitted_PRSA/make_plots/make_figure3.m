%Plot example numerical solutions with and width a pyncocline in linear
%geometry. We plot u, d, \delta p and M = U\delta T
%% Preliminaries
clear; close all
numplot = 1;
addpath('Auxillary_functions')
figpref(4) %set plot defaults
cmap = [74, 67, 176]/255; %plot colours
%% Parameters and bathymetry
run Parameters/typical_parameters %get dimensional parameters (be careful with global variable names)

% pycnocline position
x0      = 0.2;      % position of the thermocline (dimensionless)
Z0      = x0*l0;    %dimensional position of the pycnocline (Z coords)
zgl/(tau/lambda3)
%dimensionless parameters
eps1 = E0*alpha/Cd;
eps2 = E0*alpha/St;
eps3 = tau/(L/c);
eps4 = (S0 - S1)/2/S0;
eps5 = E0/2;
delta = lt/l0;
Pb = (L/c)/tau * (S0 - S1) /2 / S0 *( 1- bt*(T0 - T1)/bs / (S0 - S1));
Pt = (T0 - T1) / 2 / tau; %or Pt = (T0 - T1 + lambda1*(S0 - S1) / 2 / tau; %
lambda = (S0 + S1)/2 / S0 - bt*(L/c)/ bs / S0;
%Xmax = abs(zgl)/l0; %depth corresponding to ice shelf draft
Xmax = 1.2;
k2 = eps2/eps1;

%check supercritical near base:
if  Cd/alpha > 1
    warning('plume is subcritical near the base, may experience numerical instability')
end

% problem scales
d_scale         = E0*l0;
u_scale         = sqrt(bs*S0*g*l0*tau*E0*alpha/(L/c) / Cd);
delta_rho_scale = rho0*bs*S0;

%% Linear Ice shelf draft
zbF = @(X) X;
dzbF = @(X) 1 + 0*X;
d2zbF = @(X) 0 + 0*X;
d3zbF = @(X) 0 + 0*X;

xb = linspace(0,Xmax,1000); %bathymetry grid points
zb = zbF(xb);    %ice draft at grid points
%% Numerical solution of full equations
%solve the equations
for X0 = [l0, x0] %first loop sets pycnocline very high, i.e. as if not there   
    %solve the plume equations and return a solution structure
    tic; sol = GetPlume(eps1,eps2, eps3,eps4,delta, Pb, Pt, lambda, X0,zbF,dzbF, Xmax); toc
    
    %evaluate on a regular grid
    x = sol.x;
    x = linspace(0,x(end),1000);
    z = x*alpha;
    Y = deval(sol,x);
    d = Y(1,:);
    u = Y(2,:);
    delta_rho = Y(3,:);
    delta_T = Y(4,:); %dimensionless temperature
    Q = u.*d;
    
    % Plot solutions
    figure(numplot); clf;
    subplot(1,4,1);
    plot(d,x,'color', cmap, 'linewidth', 3);
    xlim([0,1])
    ylabel('$X$', 'interpreter', 'latex')
    xlabel('$D$', 'interpreter', 'latex')
    
    subplot(1,4,2);
    plot(u,x,'color', cmap, 'linewidth', 3);
    %ylabel('x')
    xlabel('$U$', 'interpreter', 'latex')
    
    subplot(1,4,3);
    plot(delta_rho,x,'color', cmap, 'linewidth', 3);
    %ylabel('x')
    xlabel('$\Delta \rho$', 'interpreter', 'latex')
    xlim([-0.25, 1])
    
    subplot(1,4,4); hold on
    plot(delta_T.*u,x,'color', cmap, 'linewidth', 3);
    %ylabel('x')
    xlabel('$M = U \Delta T$','interpreter', 'latex')
    if X0 == l0
        xlim([-.3,.2])
    else
        xlim([-0.1, 0.2])
    end
    
    for i = 1:4
        subplot(1,4,i); hold on
        ax = gca; %get axis of i^th subplot
        xl = ax.XLim; %xlimits on this axis
        plot(xl, [max(x), max(x)], '--', 'color', [140, 140, 140]/255) %add line xstop
        if X0 == x0
            ylim([0 0.6])
            plot(xl, [x0,x0], 'color', [140, 140, 140]/255); %plot pycnocline
        else
            ylim([0, 1.2])
        end
        yticks([0,0.5, 1])
        %xlim(xl); %re-set x limits in case new plot has adjusted them
        ax = gca;
        ax.XLabel.FontSize = 18;
        ax.YLabel.FontSize = 18;
        box on
    end
    
    %figure sizing
    fig = gcf;
    fig.Position(3:4) =[1067 293];
    numplot = numplot + 1;
end
