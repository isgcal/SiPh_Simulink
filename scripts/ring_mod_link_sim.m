%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Sen Lin <senlin@berkeley.edu>
% Integrated Systems Group, EECS, UC Berkeley
% 02/02/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open the library: SiPhotonics.slx
% open the testbench: ring_mod_link.slx
load_system('../SiPhotonics');
open_system('../ring_mod_link');
% recommend running ring_mod_analytic.m first
% get the static characteristic of the ring modulator
ring_mod_analytic;

% link parameters
P_laser = 1;    % laser power
Tc = 1;         % coupler loss
Rpd = 1;        % photodiode responsivity

% technology parameters 
% (if needed, uncomment to overwrite analytical model)
% n0 = 2;         % effective index at f0
% ng = 3;         % group index
% a0 = 300;       % absorption coefficient (/m)
% Ne = 5e17;      % N-doping (/cm^3)
% Nh = 5e17;      % P-doping (/cm^3)
% Lj = 0.5e-4;    % junction geometry factor (cm)
% a1 = -0.3;      % alpha 1st order coef. for linear/sqrt model
% n1 = 8e-5;      % neff 1st order coef. for linear/sqrt model
% Vbi = 0.8;      % built-in voltage, only for sqrt model
% f0 = 230e12;    % measurement frequency
% flaser = 229.95359e12;
% 
% % device parameters microring modulator
% L = 30e-6;  % ring perimeter
% t_in = 0.99;    % input port coupling
% t_drop = 1;      % drop port coupling


% fos + fref = flaser
% fos set to 0 for single optical channel without dispersion
fos = 0;  % laser offset frequency >= 0
fref = flaser - fos;  % simualtion reference frequency

%Vdrive_0 = 0;
%Vdrive_1 = -2;

VBIAS = Vdrive_0;
Va = Vdrive_0 - Vdrive_1;

Tb = 20e-12;     % bit width

fbw_tx = 20e9;   % transmitter bandwidth 
fbw_pd = 100e9;   % photodiode bandwidth
fbw_rx = 100e9;    % receiver bandwidth
% running simulink model
sim('ring_mod_link');
