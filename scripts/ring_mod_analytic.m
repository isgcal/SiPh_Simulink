%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Sen Lin <senlin@berkeley.edu>
% Integrated Systems Group, EECS, UC Berkeley
% 02/02/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ring_mod_analytic.m
clear all; close all;
% model selection
opt = 3;      % default diode
% 1--linear, 2--sqrt, 3--diode

% technology
n0 = 2;         % effective index at f0
ng = 3;         % group index
a0 = 100;       % waveguide absorption coefficient (/m) 4.3dB/cm (intrinsic) 
Ne = 5e17;      % N-doping (/cm^3)
Nh = 5e17;      % P-doping (/cm^3)
Lj = 0.5e-4;    % junction geometry factor (cm)
a1 = -0.3;      % alpha 1st order coef. for linear/sqrt model
n1 = 8e-5;      % neff 1st order coef. for linear/sqrt model
Vbi = 0.8;      % built-in voltage, only for sqrt model
f0 = 230e12;    % measurement frequency for n0, a0
wf = 1;        % waveguide factor
% device parameters microring modulator
L = 30e-6;        % ring perimeter in meter
t_in = 0.985;     % input port coupling
t_drop = 0.995;    % drop port coupling

% frequency sweep
fref = 229.9e12;  % simualtion reference frequency
fstart = fref;
fend = 230.2e12;
fstep = 5e8;
f_swp = fstart:fstep:fend;
lambda = 3e8./f_swp;

ring0 = RingModulator(opt, L, t_in, t_drop, ...
                n0, ng, a0, Ne, Nh, Lj, wf, a1, n1, Vbi);

Vdrive_0 = 0.5;
Vdrive_1 = -1.5;

[~,~,Pt0, Pd0] = ring0.tf(lambda, Vdrive_0);
[~,~,Pt1, Pd1] = ring0.tf(lambda, Vdrive_1);

close all;
figure(1); 
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 800, 400]);hold on; 
hold on; grid on;
set(gca, 'FontSize', 18, 'LineWidth', 2); 
set(gca,'Box','on');
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
plot(lambda*1e9, Pt0, 'r','LineWidth',3); 
plot(lambda*1e9, Pt1, 'b','LineWidth',3); 
axis([min(lambda*1e9) max(lambda*1e9) 0 1]);
xlabel('Wavelength (nm)');
ylabel('Power Transmission');
legend('bit0', 'bit1','Location', 'southwest');

% find the optimal laser wavelength (relative to resonance) to
% maximize OMA
[laser_opt, OMA_max, P1, P0] = ring0.modulation(lambda, Vdrive_0, Vdrive_1);
fprintf('-------Ring Modulator Analytical Model-------\n');
fprintf('Optimal laser wavelength:  %4.4f nm\n', 1e9*laser_opt);
fprintf('Optimal laser frequency:  %3.5f THz\n', 1e-12*3e8/laser_opt);
fprintf('Max. OMA: %1.3f \n', OMA_max);
fprintf('Bit-1 power P1: %1.3f \n', P1);
fprintf('Bit-0 power P0: %1.3f \n', P0);
flaser = 3e8/laser_opt; % used in simulink

figure(1);
plot([1e9*laser_opt, 1e9*laser_opt],[0, 1],'k-.','LineWidth',2);
plot([1e9*min(lambda), 1e9*max(lambda)],[P1, P1],'b:','LineWidth',2);
plot([1e9*min(lambda), 1e9*max(lambda)],[P0, P0],'r:','LineWidth',2);

% find FWHM bandwidth and quality factor of the ring
[~, ind_res] = min(Pt0);
ind_max = length(Pt0);
fres = f_swp(ind_max);
Pt0_l = Pt0(1:ind_res);
Pt0_r = Pt0((ind_res+1):ind_max);
f_l = f_swp(1:ind_res);
f_r = f_swp((ind_res+1):ind_max);
[~, ind_l] = min(abs((1-Pt0_l) - max(1-Pt0)/2));
[~, ind_r] = min(abs((1-Pt0_r) - max(1-Pt0)/2));
f_fwhm = f_r(ind_r) - f_l(ind_l);
Q = fres/f_fwhm;

fprintf('Q factor: %g \n', Q);
fprintf('FWHM bandwidth: %g GHz\n', f_fwhm/1e9);

% find the coupling state of the ring
[state, t_drop_critical] = ring0.coupling_state(laser_opt, Vdrive_0);
fprintf(state);
fprintf(', t_drop_nominal = %g\n', t_drop_critical);

% find the neff and alpha under different biases
[neff0, alpha0] = ring0.material(laser_opt, Vdrive_0);
[neff1, alpha1] = ring0.material(laser_opt, Vdrive_1);
fprintf('Vdrive_0 = %g V, ', Vdrive_0);
fprintf('alpha0 = %g/m, ', alpha0);
fprintf('neff0 = %g\n', neff0);
fprintf('Vdrive_1 = %g V, ', Vdrive_1);
fprintf('alpha1 = %g/m, ', alpha1);
fprintf('neff1 = %g\n', neff1);





