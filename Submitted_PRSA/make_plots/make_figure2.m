%Make figure 2 in Bradley et al, 2021. 
%Plot of the ambient temperature and salinity profiles given by
% Sa(Z) = [(Sl + Su) - (Sl - Su)tanh((Z - Z_0)/lt)]/2
% Ta(Z) = [(Tl + Tu) - (Tl - Tu)tanh((Z - Z_0)/lt)]/2

%% Preliminaries
addpath('Auxillary_functions');
figpref(4);
Z0 = -750;

%% Parameters
run parameters %get dimensional parameters (brings them all into global scope)
Sa = @(Z)((S0 + S1) - (S0 - S1)*tanh((Z - Z0)/lt))/2; %S0 = Sl, S1 = Su in paper
Ta = @(Z)((T0 + T1) - (T0 - T1)*tanh((Z - Z0)/lt))/2; %likewise for T

%% Plot
Z = linspace(-1500, 0);
figure(1); clf; 
subplot(1,2,1); hold on
fill([33.8, 35, 35, 33.8], [-850, -850, -650, -650],  'b', 'facealpha', 0.2, 'linestyle', 'none')
plot(Sa(Z), Z, 'k', 'linewidth', 2);
ylim([-1500, 0]);
xlim([33.8, 34.8]);
box on
xlabel('$S_a(Z)$~(PSU)', 'interpreter', 'latex', 'fontsize', 16);
ylabel('$Z+Z_{gl}$~(m)', 'interpreter', 'latex', 'fontsize', 16);

subplot(1,2,2); hold on
fill([-2, 1,1,-2], [-850, -850, -650, -650],  'b', 'facealpha', 0.2, 'linestyle', 'none')
plot(Ta(Z), Z, 'k', 'linewidth', 2);
ylim([-1500, 0]);
xlim([-2, 1]);
box on
xlabel('$T_a(Z)$~(${}^\circ$C)', 'interpreter', 'latex', 'fontsize', 16);
ylabel('$Z+Z_{gl}$~(m)', 'interpreter', 'latex', 'fontsize', 16);


%% Tidy and save
fig = gcf;
fig.Position(3:4) = [480, 420];
subplot(1,2,1); txta = text(33.3, 0, '(a)', 'fontsize', 16, 'interpreter', 'latex');
subplot(1,2,2); txtb = text(-3.5, 0, '(b)', 'fontsize', 16, 'interpreter', 'latex');

%saveas(gcf,'plots/figure2.png')