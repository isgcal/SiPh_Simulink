%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Sen Lin <senlin@berkeley.edu>
% Integrated Systems Group, EECS, UC Berkeley
% 02/02/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% open the library: SiPhotonics.slx
% open the testbench: mzm_link.slx
load_system('../SiPhotonics');
open_system('../mzm_link');
% recommend running mzm_analytic.m first
mzm_analytic;

% link parameters
P_laser = 1;    % laser power
Tc = 1;         % coupler loss
Rpd = 1;        % photodiode responsivity

% technology parameters 
% (if needed, uncomment to overwrite analytical model)
% n0 = 2;         % effective index at f0
% ng = 3;         % group index
% a0 = 100;       % absorption coefficient (/m) 4.3dB/cm (intrinsic)
% Ne = 5e17;      % N-doping (/cm^3)
% Nh = 5e17;      % P-doping (/cm^3)
% Lj = 0.5e-4;    % junction geometry factor (cm)
% a1 = -0.3;      % alpha 1st order coef. for linear/sqrt model
% n1 = 1e-4;      % neff 1st order coef. for linear/sqrt model
% Vbi = 0.8;      % built-in voltage, only for sqrt model
% f0 = 230e12;    % measurement frequency for n0, a0
% 
% % device parameters for MZM
% L = 1e-3;             % arm length in meter
% phase_os = -0.5*pi;   % phase offset between two arms (quarature point)

% fos + fref = flaser
% fos set to 0 for single optical channel without dispersion
fos = 0;  % laser offset frequency >= 0
fref = flaser - fos;  % simualtion reference frequency

VBIAS0 = Vdrive_0;
Va = Vdrive_0-Vdrive_1;
VBIAS1 = Vdrive_1;

Tb = 20e-12;     % bit width

fbw_tx = 100e9;   % transmitter bandwidth 
fbw_pd = 100e9;   % photodiode bandwidth
fbw_rx = 30e9;    % receiver bandwidth
% running simulink model
sim('mzm_link');
