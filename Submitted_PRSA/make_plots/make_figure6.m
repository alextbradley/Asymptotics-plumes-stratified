%Make figure 6 in the manuscript: produce a comparison between (1)
%numerical solutions of the full equations, (2) the Lazeroms
%parametrization, (3) the Lazeroms parametrization w/ ad hoc adjustment and
%(4) my parametrization for four different ice shelf basal geometries as
%specified in the manuscript
%% Preliminaries
clear 
addpath('Auxillary_functions')
figpref(4)
%clc

figure(1); clf;
figure(2); clf;
colmap = [74, 67, 176;
          81, 146, 246;
          244, 177, 115
           119,205, 156]/255; 
%% Parameters
run Parameters/typical_parameters.m %get dimensional parameters (be careful with global variable names)
Z0 = 1e10; %put the pycnocline super high so we never get anywhere near
T0 = (T0 + T1)/2;
S0 = (S0 + S1)/2; %take the mean of upper and loower values
%adjust tau and l0 to account for these new scales
tau     = T0 - (lambda1*S0 + lambda2 + lambda3*zgl);%T_{a,0} - T_{f,0}
l0      = tau/lambda3;%lengthscale of freezing pt dependence

%variable scales:
D_scale         = E0*l0;
U_scale         = sqrt(bs*S0*g*l0*tau*E0*alpha/(L/c) / Cd);
delta_rho_scale = rho0*bs*S0*tau/(L/c);
delta_T_scale   = E0 *alpha * tau/St;
X_scale         = tau/lambda3 /alpha;
    
M0 = St/(L/c);    %melt rate prefactor
lambda = (S0 + S1)/2 / S0 - bt*(L/c)/ bs / S0;
%% Ice Drafts:
N = 1e3; %Number of pts in the draft
Xb = zeros(4, N);
Zb = zeros(4, N); %initialize arrays. 4 different geometries
dZb= zeros(4, N);

% Geometry 1: linear
Xb(1,:) = linspace(eps, abs(zgl)/alpha, N);  %l0/alpha is X lengthscale
%Run from eps so that numerical solution (starts at 0) can be evaluated at
%first draft point
Zb(1,:) = alpha*Xb(1,:);
dZb(1,:) = alpha*ones(size(Xb(1,:)));

% Geometry 2: gently decreasing
X_scale = l0/alpha;
Xb(2,:) = linspace(eps, abs(zgl)/alpha*2, N);  %l0/alpha is X lengthscale (need to run to longer because takes longer to reach zgl with negative curvature)
p  = 0.45; %controls strength of curvature 
Zb(2,:) = l0*(Xb(2,:)/X_scale - p* (Xb(2,:)/X_scale).^2);
dZb(2,:) = l0*(1/X_scale - 2*p* (Xb(2,:)/X_scale^2));

% Geometry 3: sinusoidal
Xb(3,:) = linspace(eps, abs(zgl)/alpha, N);  %l0/alpha is X lengthscale
amp = 30; %amplitude of sin
per = 1e5; %period of sin
Zb(3,:) = alpha*Xb(3,:) + amp*sin(2*pi*Xb(3,:)/per);
dZb(3,:) = alpha*ones(size(Xb(3,:))) + (amp*2*pi/per)*cos(2*pi*Xb(3,:)/per);

% Geometry 4:  piecewise linear
%change slope in second half by a factor of slope2
Xb(4,:) = linspace(eps, abs(zgl)/alpha, N);  %l0/alpha is X lengthscale
lx = length(Xb);
Xb_cross = Xb(4,floor(lx/2)); %where the jump occurs
slope2 = 2; 
jump = slope2*alpha*Xb(4,floor(lx/2) + 1) - alpha*Xb(4,floor(lx/2)); %jump incurred if we dont shift second half
Zb(4,:) = [alpha*Xb(4,1:floor(lx/2)), slope2*alpha*Xb(4,floor(lx/2) + 1: end) - jump];
dZb(4,:) = [alpha*ones(1,floor(lx/2)), slope2*alpha*ones(1,lx -floor(lx/2))] ;

% Store as anonymous functions for the integral
zbFs = {@(x) x; ...
    @(x) x - p*x.^2;...
    @(x) x + (amp/l0)*sin(2*pi*x*X_scale/per);...
    @(x) (x.*(x <= Xb_cross/X_scale) + (slope2*x - jump/l0).*(x > Xb_cross/X_scale))};
dzbFs = {@(x) 1; ...
    @(x) 1 - 2*p*x; ...
    @(x) 1 + (amp/l0)*(2*pi*X_scale/per) *cos(2*pi*x*X_scale/per) ;...
    @(x) (1*(x <= Xb_cross/X_scale) + slope2*(x > Xb_cross/X_scale))};
%% Loop over each bathymetry
for i = 1:4
    sol =  GetPlumeDimensional(Ti, Si, T0, T1, S0, S1, Z0, rho0, zgl, L, ci,lt,g,Cd,...
                                      c, bs, bt, E0, St, lambda1, lambda2, lambda3, tau,...
                                      Xb(i,:), Zb(i,:), dZb(i,:));

    %process
    %find where draft pts within solution interval (should include zero)
    idx = ((Xb(i,:) > min(sol.x)) + (Xb(i,:) < max(sol.x)))>1;  
    X = Xb(i,idx);
    Y = deval(sol, Xb(i,idx));
    Z = X*alpha;
    U = Y(2,:);
    delta_T = Y(4,:); %dimensionless temperature
    Melt_rate = M0*U.*delta_T*secs_per_yr; %melt rate in metres per year
    figure(1);
    subplot(2,2,i);
    plot(Melt_rate, Z+zgl, 'color', colmap(1,:), 'linewidth', 3);
    ylim([zgl,0])
    
    %add Lazeroms/constant draft  
    x = X/X_scale;
    zbF = zbFs{1};
    dzbF = dzbFs{1}; %index 1 has constant slope and draft
    dzbFi = dzbFs{i}; 
    Q_lz = zeros(1,length(x)); %remember lower case x is dimensionless
    U_lz = zeros(1,length(x));  %initialize

    integrand = @(x) lambda^(1/3)*dzbF(x).^(4/3) .*(1 - zbF(x)).^(1/3); %integrand used in analytic solution below thermocline
    for j = 1:length(x)
        Q_lz(j) =  (2/3 *integral(integrand, 0, x(j)))^(3/2);
        U_lz(j) = lambda^(1/3) * dzbF(x(j))^(4/3) * (1 - zbF(x(j)))^(1/3) * ...
            (2/3 *integral(integrand, 0, x(j)))^(1/2)/dzbF(x(j)); %u = Q'/zb' in this region
    end
    delta_T_lz = (-Q_lz.*dzbF(x) + U_lz.*dzbF(x).*(1-zbF(x)))./U_lz;
    M_Lzrms_ND = delta_T_lz.*U_lz;
    M_Lzrms    = M0 * U_scale * delta_T_scale * M_Lzrms_ND * secs_per_yr;
    %M_Lzrms    = M_Lzrms * mean(dzbFi(x))^(3/2); %set the angle scaling in lazeroms to mean of shelf draft
    hold on
    plot(M_Lzrms, Z+zgl, 'k--');
    
    %add Lazeroms (ad-hoc adjustment)
    zbF = zbFs{i};
    dzbF = dzbFs{i}; 
    U_scale_adhoc         = sqrt(bs*S0*g*l0*tau*E0*(dzbF(x)*alpha)/(L/c) / Cd); %extra factor of alpha because dzbF is non-dimensional
    delta_T_scale_adhoc   = E0 *(dzbF(x) *alpha) * tau/St;
    M_Lzrms_adhoc    = M0 * U_scale_adhoc .* delta_T_scale_adhoc .* M_Lzrms_ND * secs_per_yr;
    plot(M_Lzrms_adhoc, Z+zgl, 'color', colmap(3,:), 'linewidth', 3);
    

    %add my integral solution using anonymous function
    zbF = zbFs{i};
    dzbF = dzbFs{i};
    Q_AB = zeros(1,length(x)); %remember lower case x is dimensionless
    U_AB = zeros(1,length(x));  %initialize

    integrand = @(x) lambda^(1/3)*dzbF(x).^(4/3) .*(1 - zbF(x)).^(1/3); %integrand used in analytic solution below thermocline
    for j = 1:length(x)
        Q_AB(j) =  (2/3 *integral(integrand, 0, x(j)))^(3/2);
        U_AB(j) = lambda^(1/3) * dzbF(x(j))^(4/3) * (1 - zbF(x(j)))^(1/3) * ...
            (2/3 *integral(integrand, 0, x(j)))^(1/2)/dzbF(x(j)); %recall u = Q'/zb' in this region
    end
    delta_T_AB = (-Q_AB.*dzbF(x) + U_AB.*dzbF(x).*(1-zbF(x)))./U_AB;
    M_AB = M0 * U_scale * delta_T_scale * U_AB .* delta_T_AB *secs_per_yr;
    subplot(2,2,i)
    plot(M_AB, Z+zgl, 'color', colmap(4,:),'linewidth', 3)
        
    %plot draft
    figure(2);
    subplot(2,2,i);
    fill([Xb(i,:), flip(Xb(i,:))], [zeros(size(Xb(i,:))), flip(Zb(i,:))], 'b')
    ylim([0, abs(zgl)])
    xlim([0, abs(zgl)/alpha])
        if i == 2 %quadratic case needs extended x scale
            [~, idx] = min(abs(Zb(i,:) - abs(zgl))); %find where Zb goes thru abs(zgl) 
            xlim([0, Xb(i,idx)])
        end
       
end

%% plot tidying
figure(1); 
fig = gcf;
fig.Position(3:4) =  [1074 640];
for i = 1:4
    subplot(2,2,i)
    yticks([-1500:500:0])
end
figure(2); 
for i = 1:4
    subplot(2,2,i); 
    yticks([]); 
    xticks([]);
end