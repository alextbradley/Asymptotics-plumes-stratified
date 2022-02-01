%Produce figure 7 in the manuscript.
%Solve equations in a flat bathymetry with ambient stratification. Compare new parametrization with
%lazeroms and numerical solutions for different positions of the
%pycnocline (a)--(c) and different values of lt (d)

%% Preliminaries
addpath('Auxillary_functions');
figpref(4);
figure(1); clf;
colmap = parula(6);
w = 0.20; 
h = 0.8;
positions = [0.07, 0.1, w, h;
             0.35, 0.1, w, h;
             0.58, 0.1, w, h;
             0.77, 0.1, w, h]; %subplot positions on page
ncols = 4;
colgap = 0.055; 
width = 0.19;
startx = (1 -width*ncols - (ncols-1)*colgap)/2;
startx = 0.06;
starty = 0.1;
height = 0.80;
positions = zeros(4,ncols);
for p = 1:ncols
positions(p,:) = [startx + (p-1)*colgap + (p-1)*width, starty , width, height];
end
%% Parameters
run parameters %get dimensional parameters 
zgl = -3000; %make artificially deeper so that transition to 0 appears
T0  = -1;
T1  = -3;
tau = T0 - (lambda1*S0 + lambda2 + lambda3*zgl);%T_{a,0} - T_{f,0}
l0  = tau/lambda3;%lengthscale of freezing pt dependence

x0s = 0.05:0.05:0.25; %dimensionless pycnocline position

%relevant variable scales:
U_scale         = sqrt(bs*S0*g*l0*tau*E0*alpha/(L/c) / Cd);
delta_T_scale   = E0 *alpha * tau/St;
X_scale         = tau/lambda3 /alpha;

%also define scales for lazeroms (uses mean of salnity and temp)
Tave = (T0 + T1)/2;
Save = (S0 + S1)/2;
tauLz = Tave - (lambda1*Save + lambda2 + lambda3*zgl);%T_{a,0} - T_{f,0}
l0Lz  = tau/lambda3;%lengthscale of freezing pt dependence
U_scaleLz = sqrt(bs*Save*g*l0Lz*tauLz*E0*alpha/(L/c) / Cd);
delta_T_scaleLz   = E0 *alpha * tauLz/St;
X_scaleLz = tauLz/lambda3/alpha;

%dimensionless parameters
eps1 = E0*alpha/Cd;
eps2 = E0*alpha/St;
eps3 = tau/(L/c);
eps4 = (S0 - S1)/2/S0;
delta = lt/l0;
Pb = (L/c)/tau * (S0 - S1) /2 / S0 *( 1- bt*(T0 - T1)/bs / (S0 - S1));
Pt = (T0 - T1) / 2 / tau; %or Pt = (T0 - T1 + lambda1*(S0 - S1) / 2 / tau; %
lambda = (S0 + S1)/2 / S0 - bt*(L/c)/ bs / S0;
Xmax = abs(zgl)/l0; %depth corresponding to ice shelf draft
M0 = St/(L/c);
k2 = eps2/eps1;

colmap = parula(length(x0s) + 1);
%% Ice shelf draft
zbF = @(X) X;
dzbF = @(X) 1 + 0*X;
d2zbF = @(X) 0*X;
d3zbF = @(X) 0*X;

xb = linspace(0,Xmax,1000); %bathymetry grid points
zb = zbF(xb);    %ice draft at grid points
%% Loop over positions of pycnocline
figure(1); clf;
subplot('Position', positions(1,:)); hold on

%add grey line indicating melt/freezing line
plot([0,0], [zgl,0], 'color', 169/255 *ones(1,3),  'HandleVisibility', 'off');

%add lazeroms/integral expression first (doesn't know anything about the pycnocline)
X_lz = linspace(0, Xmax*1.5);  %nb sol goes complex when above X = 1 (ignore)
Q_lz = zeros(1,length(X_lz));
U_lz = zeros(1,length(X_lz));  %initialize

integrand = @(x) lambda^(1/3)*dzbF(x).^(4/3) .*(1 - zbF(x)).^(1/3); %integrand used in analytic solution below thermocline
for i = 1:length(X_lz)
    Q_lz(i) =  (2/3 *integral(integrand, 0, X_lz(i)))^(3/2);
    U_lz(i) = lambda^(1/3) * dzbF(X_lz(i))^(4/3) * (1 - zbF(X_lz(i)))^(1/3) * ...
        (2/3 *integral(integrand, 0, X_lz(i)))^(1/2)/dzbF(X_lz(i)); %recall u = Q'/zb' in this region
end
delta_T_lz = (-Q_lz.*dzbF(X_lz) + U_lz.*dzbF(X_lz).*(1-zbF(X_lz)))./U_lz;
plot(M0*U_scaleLz*delta_T_scaleLz*secs_per_yr*delta_T_lz.*U_lz, X_lz*X_scaleLz*alpha + zgl, 'k--',  'HandleVisibility', 'off')

for i = 1:length(x0s)
    %solve numerically:
    sol = GetPlume(eps1,eps2, eps3,eps4,delta, Pb, Pt, lambda, x0s(i),zbF,dzbF, Xmax);
    x1 = sol.x;
    x1 = linspace(0,x1(end),1000); %regular grid to put solution on
    Y = deval(sol,x1);
    U = Y(2,:);       %dimensionless velocity
    delta_T = Y(4,:); %dimensionless temperature
    plot(M0*U_scale*delta_T_scale*secs_per_yr*U.*delta_T, x1*X_scale*alpha + zgl, 'color', colmap(i,:));
    
    %add my melt rate
    [M_AB, X_AB] = GetConstructedMeltRate(zbF, dzbF, d2zbF, d3zbF, x0s(i), Pt, Pb, delta, lambda);
    plot(M0*U_scale*delta_T_scale*secs_per_yr*M_AB, X_AB*X_scale*alpha + zgl,...
        '--', 'color', colmap(i,:),'HandleVisibility', 'off');
    
    legendinfo{i} = ['$Z_p/\ell = ' num2str(x0s(i)) '$'];
end
legend(legendinfo, 'interpreter', 'latex', 'location', 'northeast')
    
%tidy plot
xlabel('melt rate (m/yr) ', 'interpreter', 'latex', 'FontSize', 16);
ylabel('$Z + Z_{gl}$~(m)', 'interpreter', 'latex', 'FontSize', 16)
box on
ylim([zgl, 0])
yticks(-3000:1000:0)

xlim([-10,10])
ax = gca; ax.FontSize = 16;
ytickangle(90)
txta = text(-15, 300, '(a)', 'Interpreter', 'latex', 'FontSize', 16);
t1 = title({'deep grounding line', 'cold ambient'}, 'interpreter', 'latex');
%% Repeat for the more conventional parameters
run parameters %get dimensional parameters (be careful with global variable names)

%variable scales:
U_scale         = sqrt(bs*S0*g*l0*tau*E0*alpha/(L/c) / Cd);
delta_T_scale   = E0 *alpha * tau/St;
X_scale         = tau/lambda3 /alpha;

%also define scales for lazeroms (uses mean of salnity and temp)
Tave = (T0 + T1)/2;
Save = (S0 + S1)/2;
tauLz = Tave - (lambda1*Save + lambda2 + lambda3*zgl);%T_{a,0} - T_{f,0}
l0Lz  = tau/lambda3;%lengthscale of freezing pt dependence
U_scaleLz = sqrt(bs*Save*g*l0Lz*tauLz*E0*alpha/(L/c) / Cd);
delta_T_scaleLz   = E0 *alpha * tauLz/St;
X_scaleLz = tauLz/lambda3/alpha;

%add grey line indicating melt/freezing line
subplot('Position', positions(2,:)); hold on
plot([0,0], [zgl, 0], 'color', 169/255 *ones(1,3));

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
Xmax = abs(zgl)/l0; %depth corresponding to ice shelf draft
M0 = St/(L/c);    %melt rate prefactor
k2 = eps2/eps1;

%add lazeroms/integral expression first (doesn't know anything about the pycnocline)
X_lz = linspace(0, Xmax*1.5);
Q_lz = zeros(1,length(X_lz));
U_lz = zeros(1,length(X_lz));  %initialize

integrand = @(x) lambda^(1/3).*(1 - x).^(1/3); %integrand used in analytic solution below thermocline
for i = 1:length(X_lz)
    Q_lz(i) =  (2/3 *integral(integrand, 0, X_lz(i)))^(3/2);
    U_lz(i) = lambda^(1/3) * dzbF(X_lz(i))^(4/3) * (1 - zbF(X_lz(i)))^(1/3) * ...
        (2/3 *integral(integrand, 0, X_lz(i)))^(1/2)/dzbF(X_lz(i)); %recall u = Q'/zb' in this region
end
delta_T_lz = (-Q_lz.*dzbF(X_lz) + U_lz.*dzbF(X_lz).*(1-zbF(X_lz)))./U_lz;
plot(M0*U_scaleLz*delta_T_scaleLz*delta_T_lz.*U_lz*secs_per_yr, X_lz*X_scaleLz*alpha + zgl, 'k--')

for i = 1:length(x0s)
    %solve numerically:
    sol = GetPlume(eps1,eps2, eps3,eps4,delta, Pb, Pt, lambda, x0s(i),zbF,dzbF, Xmax);
    x2 = sol.x;
    x2 = linspace(0,x2(end),1000); %regular grid to put solution on
    Y = deval(sol,x2);
    U = Y(2,:);       %dimensionless velocity
    delta_T = Y(4,:); %dimensionless temperature
    plot(M0*U.*delta_T*U_scale*delta_T_scale*secs_per_yr, x2*X_scale*alpha + zgl, 'color', colmap(i,:));
    
    %add my melt rate
    [M_AB, X_AB] = GetConstructedMeltRate(zbF, dzbF, d2zbF, d3zbF, x0s(i), Pt, Pb, delta, lambda);
    plot(M0*U_scale*delta_T_scale*secs_per_yr*M_AB, X_AB*X_scale*alpha + zgl,...
        '--', 'color', colmap(i,:),'HandleVisibility', 'off');
    
end
xlabel('melt rate (m/yr) ', 'interpreter', 'latex', 'FontSize', 16);
ylabel('$Z + Z_{gl}$~(m)', 'interpreter', 'latex', 'FontSize', 16)
box on
ylim([zgl,0])
xlim([-5,12])
ax = gca; ax.FontSize = 16;
ytickangle(90)
txtb = text(-9.5, 150, '(b)', 'Interpreter', 'latex', 'FontSize', 16);
title({'typical conditions',' '}, 'interpreter', 'latex')
%% Repeat with strong stratification: force the plume to terminate
run parameters
S1 = 33.0;

%variable scales:
U_scale         = sqrt(bs*S0*g*l0*tau*E0*alpha/(L/c) / Cd);
delta_T_scale   = E0 *alpha * tau/St;
X_scale         = tau/lambda3 /alpha;

%also define scales for lazeroms (uses mean of salnity and temp)
Tave = (T0 + T1)/2;
Save = (S0 + S1)/2;
tauLz = Tave - (lambda1*Save + lambda2 + lambda3*zgl);%T_{a,0} - T_{f,0}
l0Lz  = tau/lambda3;%lengthscale of freezing pt dependence
U_scaleLz = sqrt(bs*Save*g*l0Lz*tauLz*E0*alpha/(L/c) / Cd);
delta_T_scaleLz   = E0 *alpha * tauLz/St;
X_scaleLz = tauLz/lambda3/alpha;

%add grey line indicating melt/freezing line
subplot('Position', positions(3,:)); hold on
plot([0,0], [zgl, 0], 'color', 169/255 *ones(1,3));

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
Xmax = abs(zgl)/l0; %depth corresponding to ice shelf draft
M0 = St/(L/c);    %melt rate prefactor
k2 = eps2/eps1;

%add lazeroms/integral expression first (doesn't know anything about the pycnocline)
X_lz = linspace(0, Xmax*1.5);
Q_lz = zeros(1,length(X_lz));
U_lz = zeros(1,length(X_lz));  %initialize

integrand = @(x) lambda^(1/3).*(1 - x).^(1/3); %integrand used in analytic solution below thermocline
for i = 1:length(X_lz)
    Q_lz(i) =  (2/3 *integral(integrand, 0, X_lz(i)))^(3/2);
    U_lz(i) = lambda^(1/3) * dzbF(X_lz(i))^(4/3) * (1 - zbF(X_lz(i)))^(1/3) * ...
        (2/3 *integral(integrand, 0, X_lz(i)))^(1/2)/dzbF(X_lz(i)); %recall u = Q'/zb' in this region
end
delta_T_lz = (-Q_lz.*dzbF(X_lz) + U_lz.*dzbF(X_lz).*(1-zbF(X_lz)))./U_lz;
plot(M0*U_scale*delta_T_scale*delta_T_lz.*U_lz*secs_per_yr, X_lz*X_scaleLz*alpha + zgl, 'k--')

for i = 1:length(x0s)
    %solve numerically:
    sol = GetPlume(eps1,eps2, eps3,eps4,delta, Pb, Pt, lambda, x0s(i),zbF,dzbF, Xmax);
    x2 = sol.x;
    x2 = linspace(0,x2(end),1000); %regular grid to put solution on
    Y = deval(sol,x2);
    U = Y(2,:);       %dimensionless velocity
    delta_T = Y(4,:); %dimensionless temperature
    plot(M0*U.*delta_T*U_scale*delta_T_scale*secs_per_yr, x2*X_scale*alpha + zgl, 'color', colmap(i,:));
    
    %add my melt rate
    [M_AB, X_AB] = GetConstructedMeltRate(zbF, dzbF, d2zbF, d3zbF, x0s(i), Pt, Pb, delta, lambda);
    plot(M0*U_scaleLz*delta_T_scale*secs_per_yr*M_AB, X_AB*X_scale*alpha + zgl,...
        '--', 'color', colmap(i,:),'HandleVisibility', 'off');
    
end

xlabel('melt rate (m/yr) ', 'interpreter', 'latex', 'FontSize', 16);
ylabel('$Z + Z_{gl}$~(m)', 'interpreter', 'latex', 'FontSize', 16)
box on
ylim([zgl,0])
xlim([-2, 12])
ax = gca; ax.FontSize = 16;
ytickangle(90)
txtc = text(-5.6, 150, '(c)', 'Interpreter', 'latex', 'FontSize', 16);
title({'typical grounding line', 'strong ambient salinity strat.'}, 'interpreter', 'latex')
% saveas(gcf,'plots/figure7.png')

%% Vary lp
run parameters %get dimensional parameters (be careful with global variable names)
lt = [25, 50, 100, 150, 200];
%variable scales:
U_scale         = sqrt(bs*S0*g*l0*tau*E0*alpha/(L/c) / Cd);
delta_T_scale   = E0 *alpha * tau/St;
X_scale         = tau/lambda3 /alpha;

%also define scales for lazeroms (uses mean of salnity and temp)
Tave = (T0 + T1)/2;
Save = (S0 + S1)/2;
tauLz = Tave - (lambda1*Save + lambda2 + lambda3*zgl);%T_{a,0} - T_{f,0}
l0Lz  = tau/lambda3;%lengthscale of freezing pt dependence
U_scaleLz = sqrt(bs*Save*g*l0Lz*tauLz*E0*alpha/(L/c) / Cd);
delta_T_scaleLz   = E0 *alpha * tauLz/St;
X_scaleLz = tauLz/lambda3/alpha;

%add grey line indicating melt/freezing line
subplot('Position', positions(4,:)); hold on
plot([0,0], [zgl, 0], 'color', 169/255 *ones(1,3), 'HandleVisibility', 'off');

%dimensionless parameters
eps1 = E0*alpha/Cd;
eps2 = E0*alpha/St;
eps3 = tau/(L/c);
eps4 = (S0 - S1)/2/S0;
eps5 = E0/2;
deltas = lt/l0;
Pb = (L/c)/tau * (S0 - S1) /2 / S0 *( 1- bt*(T0 - T1)/bs / (S0 - S1));
Pt = (T0 - T1) / 2 / tau; %or Pt = (T0 - T1 + lambda1*(S0 - S1) / 2 / tau; %
lambda = (S0 + S1)/2 / S0 - bt*(L/c)/ bs / S0;
Xmax = abs(zgl)/l0; %depth corresponding to ice shelf draft
M0 = St/(L/c);    %melt rate prefactor
k2 = eps2/eps1;
x0 = 0.15; %constant x0 in this case

%add lazeroms/integral expression first (doesn't know anything about the pycnocline)
X_lz = linspace(0, Xmax*1.5);
Q_lz = zeros(1,length(X_lz));
U_lz = zeros(1,length(X_lz));  %initialize

integrand = @(x) lambda^(1/3).*(1 - x).^(1/3); %integrand used in analytic solution below thermocline
for i = 1:length(X_lz)
    Q_lz(i) =  (2/3 *integral(integrand, 0, X_lz(i)))^(3/2);
    U_lz(i) = lambda^(1/3) * dzbF(X_lz(i))^(4/3) * (1 - zbF(X_lz(i)))^(1/3) * ...
        (2/3 *integral(integrand, 0, X_lz(i)))^(1/2)/dzbF(X_lz(i)); %recall u = Q'/zb' in this region
end
delta_T_lz = (-Q_lz.*dzbF(X_lz) + U_lz.*dzbF(X_lz).*(1-zbF(X_lz)))./U_lz;
plot(M0*U_scaleLz*delta_T_scaleLz*delta_T_lz.*U_lz*secs_per_yr, X_lz*X_scaleLz*alpha + zgl, 'k--','HandleVisibility', 'off')

for i = 1:length(deltas)
    %solve numerically:
    delta = deltas(i);
    sol = GetPlume(eps1,eps2, eps3,eps4,delta, Pb, Pt, lambda, x0,zbF,dzbF, Xmax);
    x2 = sol.x;
    x2 = linspace(0,x2(end),1000); %regular grid to put solution on
    Y = deval(sol,x2);
    U = Y(2,:);       %dimensionless velocity
    delta_T = Y(4,:); %dimensionless temperature
    plot(M0*U.*delta_T*U_scale*delta_T_scale*secs_per_yr, x2*X_scale*alpha + zgl, 'color', colmap(i,:));
    
    %add my melt rate
    [M_AB, X_AB] = GetConstructedMeltRate(zbF, dzbF, d2zbF, d3zbF, x0, Pt, Pb, delta, lambda);
    plot(M0*U_scale*delta_T_scale*secs_per_yr*M_AB, X_AB*X_scale*alpha + zgl,...
        '--', 'color', colmap(i,:),'HandleVisibility', 'off');
    
    legendinfo{i} = ['$\ell_p  = ' num2str(lt(i)) '~m$'];
    
end
legend(legendinfo, 'interpreter', 'latex', 'location', 'northeast')

xlabel('melt rate (m/yr) ', 'interpreter', 'latex', 'FontSize', 16);
ylabel('$Z + Z_{gl}$~(m)', 'interpreter', 'latex', 'FontSize', 16)
box on
ylim([zgl,0])
xlim([-5,12])
ax = gca; ax.FontSize = 16;
ytickangle(90)
txtd = text(-9.5, 150, '(d)', 'Interpreter', 'latex', 'FontSize', 16);
title({'typical conditions', 'vary $\ell_p$'}, 'interpreter', 'latex')
fig = gcf; fig.Position(3:4) = [1100, 482];

