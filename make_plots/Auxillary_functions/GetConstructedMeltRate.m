function [MM,XX] = GetConstructedMeltRate(zbF, dzbF, d2zbF, d3zbF, Z0, Pt, Pb, delta, lambda)
%return the dimensionless melt rate according to the algorithm described in
%the paper. Code also outputs the melt rate in each of the regions
%described in the paper.
%Inputs:
% zbF   :   Anonymous function
%           Shape of ice draft
%
% dzbF  :   Anonymous function
%           Derivative of shape of ice draft
%
% d2zbf :   Anonymous function
%           2nd derivative of shape of ice draft
%
% d3zbF :   Anonymous function
%           3rd derivative of shape of ice draft
%
% Z0    :   Scalar
%           Dimensionless center position of pycnocline
%
% Pt    :   Scalar
%           Strength of change in thermal driving across pycnocline
%
% Pb    :   Scalar
%           Strength of change in buoyancy deficit across pycnocline
%
% delta :   Scalar
%           delta = lt/l0, ratio of thermocline lengthscale and lengthscale
%           of buoyancy deficit
%
% lambda:   Scalar (= kappa in the ms)

%% Find pycnocline position in terms of X
XX = linspace(0,1,1e2); 
[~,idx] = min(abs(zbF(XX) - Z0));
X0 = XX(idx);

%% Below pycnocline
N = 2; %number of length scales to match over

X_below = linspace(eps, X0 - N*delta,300); %lower co-ordinates
Q_below = zeros(1,length(X_below));
U_below = zeros(1,length(X_below));  %initialize

integrand = @(x) lambda^(1/3)*dzbF(x).^(4/3) .*(1 - zbF(x)).^(1/3); %integrand used in analytic solution below thermocline
for i = 1:length(X_below)
    %Q_below(i) =  (2/3 *integral(integrand, 0, X_below(i)))^(3/2);
    q = [0, X_below(1:i)];
    qq = [0,integrand(X_below(1:i))];
    Q_below(i) =  (2/3 *trapz(q, qq))^(3/2);
    %U_below(i) = lambda^(1/3) * dzbF(X_below(i))^(4/3) * (1 - zbF(X_below(i)))^(1/3) * ...
    %    (2/3 *integral(integrand, 0, X_below(i)))^(1/2)/dzbF(X_below(i)); %recall u = Q'/zb' in this region
    U_below(i) = lambda^(1/3) * dzbF(X_below(i))^(4/3) * (1 - zbF(X_below(i)))^(1/3) * ...
        (2/3 *trapz(q, qq))^(1/2)/dzbF(X_below(i)); %recall u = Q'/zb' in this region
end
D_below = Q_below./U_below;
delta_rho_below = U_below.^2./D_below ./dzbF(X_below);
delta_T_below = (-Q_below.*dzbF(X_below) + U_below.*dzbF(X_below).*(1-zbF(X_below)))./U_below;

%% pycnocline
%quanties into and out of thermocline
X_in          = X_below(end);
D_in          = D_below(end);
U_in          = U_below(end);
delta_rho_in  = delta_rho_below(end);
Q_in          = Q_below(end);
delta_T_in    = delta_T_below(end);

X_out         = X0 + N*delta;
Q_out         = Q_in; %naive: Q is continuous across thermocline
Q_out         = Q_in + U_in*dzbF(X0)*(X_out - X_in); %Taylor expansion across thermocline

delta_rho_out = delta_rho_in  - 2*Pb*dzbF(X0);
U_out         = (Q_out*delta_rho_out*dzbF(X0))^(1/3);
D_out         = Q_out/U_out;
delta_T_out   = (-Q_out*dzbF(X_out) + U_out*dzbF(X_out)*(1 - 2*Pt-zbF(X_out)))/U_out;

%linearly interpolate between quantities
X_thermocline = linspace(X_in, X_out);
U_thermocline = U_in + (X_thermocline - X_in)*(U_out - U_in)/(X_out - X_in);
d_thermocline = D_in + (X_thermocline - X_in)*(D_out - D_in)/(X_out - X_in);
delta_rho_thermocline = delta_rho_in + (X_thermocline - X_in)*(delta_rho_out - delta_rho_in)/(X_out - X_in);
delta_T_thermocline = delta_T_in + (X_thermocline - X_in)*(delta_T_out - delta_T_in)/(X_out - X_in);

%% Above thermocline
% Construct quantities arising from linearization
f = 0.7;

%lambda1 = (1 - zbF(X_out))*Q_out - U_out^3 / lambda / dzbF(X_out);
%lambda2 = 0;
%K1      = lambda^(1/3) * (lambda2 + dzbF(X_out)^4*((1-zbF(X_out))*Q_out - lambda1))^(1/3);
%K2      = lambda*(4*dzbF(X_out)^3 * d2zbF(X_out)*((1 - zbF(X_out))*Q_out - lambda1) + ...
%    dzbF(X_out)^4*((1 - zbF(X_out) - 2*Pt)*K1 - dzbF(X_out)*Q_out))/3/K1^2;
%K3hat   = dzbF(X_out)^4 * ((1 - zbF(X_out) - 2*Pt)*K2/2 - dzbF(X_out)*K1 - Q_out*d2zbF(X_out)/2) + ...
%    4*d2zbF(X_out)*dzbF(X_out)^3 * ((1-zbF(X_out)- 2*Pt)*K1 - dzbF(X_out)*Q_out) +...
%    (2 * dzbF(X_out)*d3zbF(X_out) + 6*dzbF(X_out)^2*d2zbF(X_out)^2)*((1 - zbF(X_out))*Q_out - lambda1);
%K3      = (lambda*K3hat - 3*K1*K2^2)/3/K1^2;
%[K3, K2, K1 - U_out/2 *dzbF(X0)]
%r = roots([K3, K2, K1 - U_out*f *dzbF(X0)]);
%r = r(r>0);


K1 = dzbF(X_out) * U_out;
K2 = (4 * dzbF(X_out)^2 * d2zbF(X_out) * U_out^3 + lambda * dzbF(X_out)^4 *( (1 - zbF(X_out) - 2*Pt)*K1 -dzbF(X_out)*Q_out))/6/K1^2;
K3hat = lambda*dzbF(X_out)^4 * ((1- zbF(X_out) - 2*Pt)*K2 - d2zbF(X_out)*Q_out / 2 - K1 * dzbF(X_out)) + ...
        4 * lambda * dzbF(X_out)^3 * d2zbF(X_out) * ((1 - zbF(X_out) - 2*Pt)*K1 - dzbF(X_out)*Q_out) + ...
        2 * U_out^3 *( dzbF(X_out)^2 * d3zbF(X_out) + 3*dzbF(X_out)*d2zbF(X_out)^2);
K3 = (K3hat - 3*K1 *K2)/9/K1^2;


%solve the polynomial for xstar
r = roots([3*K3, 2*K2, K1 - U_out*f *dzbF(X0)]);
r = r(r>0);
if any(abs(imag(r)) <1e-10)
    r = max(r);
    xstar = X0 + r;
else %take the smallest value U acheives at x_star
    xstar = 10; %take to be large
end
%pause

% Compile terms in expansion
vareps  = xstar - X0; %parameter initicating how far beyong x_out we can go
X_upper_inner = linspace(0,1, 1e2); %scaled variable for upper
X_upper = X_out + vareps*X_upper_inner; %we'll chop this down later, when the  appropriate x_star (where u = 1/2 u_out occurs) has been determined

Q_upper_first_term = K1*X_upper_inner;
Q_upper_second_term = K2*X_upper_inner.^2 ;
Q_upper_third_term = K3*X_upper_inner.^3;
Q_upper = Q_out + vareps*Q_upper_first_term +...
    vareps^2 * Q_upper_second_term + vareps^3 * Q_upper_third_term;

U_upper         = (K1 + vareps*2*K2*X_upper_inner + vareps^2 * 3*K3 * X_upper_inner.^2)./dzbF(X_upper); %note lower case x_upper in denominator
D_upper         = Q_upper./U_upper;
delta_rho_upper = U_upper.^2 ./ D_upper ./dzbF(X_upper);
delta_T_upper  = (1 - 2*Pt - zbF(X_upper)).*dzbF(X_upper)...
    - D_upper.* dzbF(X_upper);


%% Transition region   
%use asymptotic behaviour to determine xc, where the plume speed goes
%to zero:
% we know that u ~ pref*(xc - x)^(1/3) as xc - x -> 0 (analytically we have
% const = Qc^(1/3)*zb'(xc)^(2/3) but we don't know Qc, xc a priori). Here
% we find pref (prefactor) and xc so that u is ctsly differentiable across
% x_star. 

%Compute quantities needed for matching
x_star = X_upper(end);
U_star = U_upper(end);
dU_dX_star = (U_upper(end) - U_upper(end-1))/(X_upper(end) - X_upper(end-1)); %derivative of u at x = x_star (one sided FD)
xc = x_star - U_star/3/dU_dX_star; %dividing cty equations gives xc
pref_U = U_star/(xc - x_star)^(1/3);  %continuity of u then gives prefactor


%assemble transition region quantities
X_transition = linspace(x_star, xc-1e-6, 1e3); %large number of points to get close to cusp
U_transition = pref_U*(xc - X_transition).^(1/3); %add a small quantity to prevent d going to inf
D_transition = Q_upper(end)./U_transition; %Q constant across transition layer
delta_T_transition  = (1 - 2*Pt - zbF(X_transition)).*dzbF(X_transition)...
    - D_transition.* dzbF(X_transition); %based on the exact relationship in the above thermocline region

Q_transition = D_upper(end).*U_transition;
M_transition =  (1 - 2*Pt - zbF(X_transition)).*dzbF(X_transition).*U_transition ...
    - Q_transition.* dzbF(X_transition);

%if xc < x_star, transition region unphysical. In this case, set transition
%region to be trivial and extend the upper region to the top of the draft. 
if xc < x_star
    % Compile terms in expansion
    vareps  = 2 - X0; %parameter initicating how far beyong x_out we can go
    X_upper_inner = linspace(0,1, 1e3); %scaled variable for upper
    X_upper = X_out + vareps*X_upper_inner; %we'll chop this down later, when the  appropriate x_star (where u = 1/2 u_out occurs) has been determined
    
    Q_upper_first_term = K1*X_upper_inner;
    Q_upper_second_term = K2*X_upper_inner.^2;
    Q_upper_third_term = K3*X_upper_inner.^3;
    Q_upper = Q_out + vareps*Q_upper_first_term +...
        vareps^2 * Q_upper_second_term + vareps^3 * Q_upper_third_term;
    
    U_upper         = (K1 + vareps*2*K2*X_upper_inner + vareps^2 * 3*K3 * X_upper_inner.^2)./dzbF(X_upper); %note lower case x_upper in denominator
    D_upper         = Q_upper./U_upper;
    delta_rho_upper = U_upper.^2 ./ D_upper ./dzbF(X_upper);
    delta_T_upper  = (1 - 2*Pt - zbF(X_upper)).*dzbF(X_upper)...
        - D_upper.* dzbF(X_upper);
    U_transition = [];
    D_transition = [];
    X_transition = [];
    delta_T_transition = [];
    M_transition = [];
end

%% Assemble arrays
%add small melt above transition
tidx = length([X_below, X_thermocline]); %idx where upper takes over;
X       = [X_below, X_thermocline, X_upper, X_transition,1];
delta_T = [delta_T_below, delta_T_thermocline, delta_T_upper, delta_T_transition];
delta_T(end) = delta_T(end-2); %because of how we have constructed, delta_T -> inf as U -> 0. Remove this singularity here
delta_T(end-1) = delta_T(end);
delta_T = [delta_T, 0];
U       = [U_below, U_thermocline, U_upper, U_transition, 0];

%apply the lazeroms geometric fudge factor
U(tidx + 1:end) = U(tidx + 1:end).*sqrt(dzbF(X(tidx+1:end)));
delta_T(tidx + 1:end) = delta_T(tidx + 1:end).*(dzbF(X(tidx+1:end)));

M1 = delta_T_below.*U_below;
X1 = X_below;
M2 = delta_T_thermocline.*U_thermocline;
X2 = X_thermocline;
M3l = delta_T_upper.*U_upper;
X3l = X_upper;
M3u = M_transition;
X3u = X_transition;
MM = [M1, M2, M3l, M3u];
XX = [X1, X2, X3l, X3u];
%% Add caveat for separation within pycnocline
if delta_rho_out < 0
    %find last positive velocity point in pycnocline (you could do this by
    %solving simple linear equation)
    U_out         = 0;
    delta_T_out   = 0;

    %linearly interpolate between quantities
    X_thermocline = linspace(X_in, X_out);
    U_thermocline = U_in + (X_thermocline - X_in)*(U_out - U_in)/(X_out - X_in);
    delta_T_thermocline = delta_T_in + (X_thermocline - X_in)*(delta_T_out - delta_T_in)/(X_out - X_in);
    M2 = delta_T_thermocline.*U_thermocline;

    MM = [M1, M2];
    XX = [X1, X2];
end
