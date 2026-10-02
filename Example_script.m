%% Reset-control toolbox example
% This script accompanies the toolbox chapter of the PhD thesis and shows
% how the analysis and design tools developed in the thesis can be applied
% to one simple closed-loop example.
%
% The example uses a mass plant with time delay. A conventional PID
% controller is first designed as a linear baseline. A PI+CgLp controller is
% then added and the following tools are demonstrated:
%
%   1) CgLp design using calc_omega_f
%   2) Reset HOSIDFs using computeResetHOSIDF_SF
%   3) Closed-loop pseudo-sensitivity using computePseudoSens_SF
%   4) Stability assessment using theta_deg
%   5) Nonlinearity quantification using sigma_2
%   6) Nonlinearity shaping using complementary pre/post filters
%   7) Reliability assessment using sigma_t and sigma_d
%   8) Shaping-filter design using design_Cs_maxPhase
%
% The script is intentionally kept as one sequence of MATLAB commands so
% that each design step can be followed and modified directly.
%
% Required MATLAB functions:
%   calc_omega_f.m
%   computeResetHOSIDF_SF.m
%   computePseudoSens_SF.m
%   computePseudoSens_SF_withMetrics.m
%   convertToLure_SF.m
%   theta_deg.m
%   design_Cs_maxPhase.m
%
% Required MATLAB toolboxes:
%   Control System Toolbox
%   Optimization Toolbox (for design_Cs_maxPhase)

clear;
clc;
close all;

%% Plot settings
% LaTeX formatting is used for all axes and legends. The same frequency
% range is used in every frequency-domain figure.
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',13);
set(groot,'defaultLineLineWidth',1.8);

% Colors used throughout the example.
cPID   = [0.15 0.15 0.15];
cReset = [0.00 0.40 0.65];
cThird = [0.60 0.10 0.10];
cShape = [0.50 0.30 0.65];
cGray  = [0.60 0.60 0.60];

% Fix the random seed because design_Cs_maxPhase uses random initial
% guesses in its multi-start optimization.
rng(1);

s = tf('s');

%% Frequency grid
% All frequency-response calculations and all frequency-domain figures use
% this same grid. The grid is linear because the reset-analysis functions
% use harmonic indices based on the frequency-vector spacing.
f = 0.1:0.1:2000;              % Frequency [Hz]
w = 2*pi*f;                    % Angular frequency [rad/s]

%% 1. Mass plant with time delay
% A mass is selected as a transparent benchmark because its dynamics are
% dominated by the double integrator, while the delay introduces a clear
% bandwidth and phase-margin limitation. This makes the influence of the
% reset controller easy to interpret without introducing unnecessary plant
% dynamics.
m  = 1;                        % Mass [kg]
Td = 0.5e-3;                   % Time delay [s]

P = 1/(m*s^2);
P.InputDelay = Td;

f_BW = 100;                    % Desired crossover frequency [Hz]
w_BW = 2*pi*f_BW;              % Desired crossover frequency [rad/s]
PM   = 30;                     % Desired phase margin [deg]

%% 2. Baseline PID controller
% A filtered PID is designed by loop shaping. The PI zero is placed one
% decade below the desired crossover frequency. A lead section is then
% centered around the crossover and supplies the phase required to obtain
% the selected phase margin.

% PI corner.
w_i = w_BW/10;

% Plant phase at the desired crossover frequency.
P_BW = reshape(squeeze(freqresp(P,w_BW)),1,[]);
phi_P = rad2deg(angle(P_BW));

% MATLAB can represent -180 deg as +180 deg. Shift to the branch used for
% the loop-shaping calculation.
if phi_P > 0
    phi_P = phi_P - 360;
end

% Phase contribution of the PI term and required lead contribution.
phi_PI   = -atand(w_i/w_BW);
phi_lead = (-180 + PM) - phi_P - phi_PI;

% Lead zero and pole, with the maximum phase lead centered at w_BW.
alpha = (1-sind(phi_lead))/(1+sind(phi_lead));
w_z = w_BW*sqrt(alpha);
w_p = w_BW/sqrt(alpha);

% Controller gain is selected so that |P(jw_BW)C(jw_BW)| = 1.
C_PID_0 = (1 + w_i/s)*(1 + s/w_z)/(1 + s/w_p);
K_PID = 1/abs(reshape(squeeze(freqresp(P*C_PID_0,w_BW)),1,[]));
C_PID = minreal(K_PID*C_PID_0);

L_PID = P*C_PID;

[~,PM_PID,~,Wcp_PID] = margin(L_PID);

fprintf('\nPID design\n');
fprintf('----------\n');
fprintf('Crossover frequency : %.2f Hz\n',Wcp_PID/(2*pi));
fprintf('Phase margin        : %.2f deg\n',PM_PID);
fprintf('PI corner           : %.2f Hz\n',w_i/(2*pi));
fprintf('Lead zero           : %.2f Hz\n',w_z/(2*pi));
fprintf('Lead pole           : %.2f Hz\n',w_p/(2*pi));

% Open-loop frequency response.
L_PID_f = reshape(squeeze(freqresp(L_PID,w)),1,[]);
phase_PID = rad2deg(unwrap(angle(L_PID_f)));

if phase_PID(1) > 0
    phase_PID = phase_PID - 360;
end

figure('Color','w');
subplot(211)
semilogx(f,20*log10(abs(L_PID_f)),'Color',cPID);
hold on;
yline(0,'--','Color',cGray);
xline(f_BW,'--','Color',cGray);
ylabel('Magnitude (dB)');
xlim([0.1 2000]);
grid on;

subplot(212)
semilogx(f,phase_PID,'Color',cPID);
hold on;
xline(f_BW,'--','Color',cGray);
xlabel('Frequency (Hz)');
ylabel('Phase (degree)');
xlim([0.1 2000]);
grid on;

%% Reset-analysis settings
% nrHOSIDFsMax limits the number of harmonics used in the closed-loop
% reconstruction. Fs and overSampleRatio determine the time resolution used
% by computePseudoSens_SF and computePseudoSens_SF_withMetrics.
nrHOSIDFsMax = 31;
Fs = 2000;                     % Sampling frequency [Hz]
overSampleRatio = 10;

%% 3. PI+CgLp add-on controller
% The additional PI increases the low-frequency loop gain. The CgLp element
% is designed to compensate the phase lag introduced by this PI at the
% desired crossover frequency.

A_rho = 0;                     % Full reset

% Additional PI controller.
f_i_add = 15;                  % PI corner [Hz]
w_i_add = 2*pi*f_i_add;
C_PI_add = 1 + w_i_add/s;

% CgLp lead corner.
f_l = 70;                      % Lead corner [Hz]
w_l = 2*pi*f_l;

% The CgLp phase contribution required at the crossover frequency is chosen
% to compensate the phase lag introduced by the additional PI.
theta_CgLp = atand(w_i_add/w_BW);

% No shaping filter is used for the initial CgLp design.
frfCs = ones(size(f));

% Calculate the CgLp filter frequency w_f. The function uses the frequency
% vector f [Hz], the target bandwidth w_BW [rad/s], the requested CgLp phase
% [deg], the lead corner w_l [rad/s], the reset coefficient A_rho, and the
% frequency response of the shaping filter C_s.
w_f = calc_omega_f(f,w_BW,theta_CgLp,w_l,A_rho,frfCs);

% CgLp parameters.
w_r = w_l/sqrt(1 + (4*(1-A_rho)/(pi*(1+A_rho)))^2);
D_r = w_l/(w_f-w_l);
k_c = (w_f-w_l)/w_f;
C_c = (1+s/w_l)/(1+s/w_f);

% Strictly proper reset element. In computePseudoSens_SF, the direct term
% D_r is represented by the parallel linear path C_par.
sysR = ss(-w_r,1,w_r,0);

% Blocks used in the closed-loop representation.
C_pos = k_c*C_c*C_PI_add*C_PID;

fprintf('\nPI+CgLp design\n');
fprintf('--------------\n');
fprintf('Added PI corner       : %.2f Hz\n',w_i_add/(2*pi));
fprintf('Required CgLp phase   : %.2f deg\n',theta_CgLp);
fprintf('w_l/(2*pi)            : %.2f Hz\n',w_l/(2*pi));
fprintf('w_r/(2*pi)            : %.2f Hz\n',w_r/(2*pi));
fprintf('w_f/(2*pi)            : %.2f Hz\n',w_f/(2*pi));
fprintf('D_r                   : %.4f\n',D_r);

%% 3.1 Reset HOSIDFs and CgLp HOSIDFs
% The first- and third-order HOSIDFs are sufficient here to visualize the
% fundamental response and the dominant higher-order nonlinear contribution.

% Reset-element HOSIDFs, including the direct term D_r.
H1 = computeResetHOSIDF_SF( ...
    sysR.A,sysR.B,sysR.C,D_r,A_rho,f,1,frfCs);
H3 = computeResetHOSIDF_SF( ...
    sysR.A,sysR.B,sysR.C,D_r,A_rho,f,3,frfCs);

% Linear CgLp filter evaluated at the corresponding harmonic frequencies.
frfCc_1 = reshape(squeeze(freqresp(C_c,w)),1,[]);
frfCc_3 = reshape(squeeze(freqresp(C_c,3*w)),1,[]);

CgLp_1 = k_c*frfCc_1.*H1;
CgLp_3 = k_c*frfCc_3.*H3;

phase_CgLp = rad2deg(unwrap(angle(CgLp_1)));
phase_CgLp = phase_CgLp - 360*round(phase_CgLp(1)/360);

figure('Color','w');
subplot(211)
semilogx(f,20*log10(abs(CgLp_1)),'Color',cReset);
hold on;
semilogx(f,20*log10(abs(CgLp_3)),'--','Color',cThird);
yline(0,':','Color',cGray);
ylabel('Magnitude (dB)');
legend('$|\mathfrak{C}_1|$','$|\mathfrak{C}_3|$','Location','best');
xlim([0.1 2000]);
grid on;

subplot(212)
semilogx(f,phase_CgLp,'Color',cReset);
xlabel('Frequency (Hz)');
ylabel('Phase (degree)');
xlim([0.1 2000]);
grid on;

frfC_PI_add = reshape(squeeze(freqresp(C_PI_add,w)),1,[]);

frfC_PI_add_cgLp=frfC_PI_add.*CgLp_1;

figure('Color','w');

subplot(211)
semilogx(f,20*log10(abs(frfC_PI_add)), ...
    'Color',cPID,'LineWidth',1.5);
hold on;
semilogx(f,20*log10(abs(frfC_PI_add_cgLp)), ...
    'Color',cReset,'LineWidth',1.5);

ylabel('Magnitude (dB)','Interpreter','latex');

legend('PI', 'PI+CgLp', ...
    'Interpreter','latex', ...
    'Location','best');

xlim([0.1 2000]);
grid on;


subplot(212)

phase_PI_add = rad2deg(unwrap(angle(frfC_PI_add)));
phase_PI_add_cglp  = rad2deg(unwrap(angle(frfC_PI_add_cgLp)));

semilogx(f,phase_PI_add, ...
    'Color',cPID,'LineWidth',1.5);
hold on;
semilogx(f,phase_PI_add_cglp, ...
    'Color',cReset,'LineWidth',1.5);

xlabel('Frequency (Hz)','Interpreter','latex');
ylabel('Phase (degree)','Interpreter','latex');

xlim([0.1 2000]);
grid on;

%% 3.2 Closed-loop pseudo-sensitivity
% Convert the LTI blocks to frequency-response arrays on the same frequency
% grid used by the reset-analysis functions.
frfPlant = reshape(squeeze(freqresp(P,w)),1,[]);
frfCpre  = ones(size(f));
frfCpar  = D_r*ones(size(f));
frfCpos  = reshape(squeeze(freqresp(C_pos,w)),1,[]);

[pseudoMag,~,Swz,~,~] = computePseudoSens_SF( ...
    f,sysR,A_rho,frfCs,frfCpre,frfCpar,frfCpos,frfPlant, ...
    Fs,nrHOSIDFsMax,overSampleRatio);

% Linear PID sensitivity for comparison.
S_PID = 1./(1+L_PID_f);

figure('Color','w');
semilogx(f,20*log10(abs(S_PID)),'Color',cPID);
hold on;
semilogx(f,20*log10(pseudoMag),'Color',cReset);
semilogx(f,20*log10(abs(Swz(1,:))),'--','Color',cShape);
semilogx(f,20*log10(abs(Swz(3,:))),'Color',cThird);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('$|S_{\mathrm{PID}}|$','$|S_{\infty}|$','$|S_{r,e}^{1}|$','$|S_{r,e}^{3}|$', ...
    'Location','best');
xlim([0.1 2000]);
grid on;

%% 4. Stability assessment
% theta_deg evaluates the angle theta_N used by the frequency-domain
% stability condition developed in the thesis. The condition requires the
% total angular spread to remain below 180 deg and theta_N to remain inside
% one of the admissible angular sectors.

% Base-linear reset element, including D_r.
R_bl = w_r/(s+w_r)+D_r;
frfRbl = reshape(squeeze(freqresp(R_bl,w)),1,[]);

% No shaping filter is used in the reset-trigger path for the controller
% designed in Section 3, hence C_s = 1.
theta_N = theta_deg( ...
    frfPlant,ones(size(f)),frfCpos,zeros(size(f)), ...
    ones(size(f)),frfRbl,D_r);

% Express theta_N on the branch used by the stability condition,
% theta_N in [-90,270) deg.
theta_N(theta_N < -90) = theta_N(theta_N < -90) + 360;
theta_N(theta_N >= 270) = theta_N(theta_N >= 270) - 360;

theta_min = min(theta_N);
theta_max = max(theta_N);
theta_spread = theta_max-theta_min;

% Sufficient frequency-domain stability condition for the first-order
% reset element: the angular spread must be below 180 deg and theta_N must
% remain inside one of the admissible angular sectors.
stable_reset = (theta_spread < 180) && ...
    (all(theta_N>-90 & theta_N<180) || ...
     all(theta_N>0 & theta_N<270));

fprintf('\nStability assessment\n');
fprintf('--------------------\n');
fprintf('Certified           : %d\n',stable_reset);
fprintf('Angular spread      : %.2f deg\n',theta_spread);
fprintf('Minimum theta_N     : %.2f deg\n',theta_min);
fprintf('Maximum theta_N     : %.2f deg\n',theta_max);

figure('Color','w');
semilogx(f,theta_N,'Color',cReset);
hold on;
yline(-90,':','Color',cGray);
yline(0,':','Color',cGray);
yline(180,':','Color',cGray);
yline(270,':','Color',cGray);
xlabel('Frequency (Hz)');
ylabel('$\theta_{\mathcal{N}}$ (degree)');
legend('$C_s=1$','Location','best');
xlim([0.1 2000]);
grid on;

%% 5. Robust performance and sigma_2
% To make the influence of higher-order harmonics visible, a more aggressive
% additional PI is used. The same analysis procedure is then repeated.

f_i_aggr = 30;                 % Added PI corner [Hz]
w_i_aggr = 2*pi*f_i_aggr;
C_PI_aggr = 1 + w_i_aggr/s;

f_l_aggr = 50;                 % CgLp lead corner [Hz]
w_l_aggr = 2*pi*f_l_aggr;
theta_CgLp_aggr = atand(w_i_aggr/w_BW);

% Start the aggressive design without a reset-trigger shaping filter.
% This keeps C_s = 1 until the reliability-shaping step in Section 6.
frfCs_aggr = ones(size(f));

w_f_aggr = calc_omega_f( ...
    f,w_BW,theta_CgLp_aggr,w_l_aggr,A_rho,frfCs_aggr);

w_r_aggr = w_l_aggr/sqrt( ...
    1 + (4*(1-A_rho)/(pi*(1+A_rho)))^2);
D_r_aggr = w_l_aggr/(w_f_aggr-w_l_aggr);
k_c_aggr = (w_f_aggr-w_l_aggr)/w_f_aggr;
C_c_aggr = (1+s/w_l_aggr)/(1+s/w_f_aggr);

sysR_aggr = ss(-w_r_aggr,1,w_r_aggr,0);
C_pos_aggr = k_c_aggr*C_c_aggr*C_PI_aggr*C_PID;

frfCpre_aggr = ones(size(f));
frfCpar_aggr = D_r_aggr*ones(size(f));
frfCpos_aggr = reshape(squeeze(freqresp(C_pos_aggr,w)),1,[]);

[pseudoMag_aggr,~,Swz_aggr,~,~] = computePseudoSens_SF( ...
    f,sysR_aggr,A_rho,frfCs_aggr, ...
    frfCpre_aggr,frfCpar_aggr,frfCpos_aggr,frfPlant, ...
    Fs,nrHOSIDFsMax,overSampleRatio);

% sigma_2 measures the contribution of higher-order harmonics relative to
% the first-order closed-loop response. The result is expressed in percent.
sigma_2 = (sqrt(sum(abs(Swz_aggr).^2,1,'omitnan')) ...
    - abs(Swz_aggr(1,:)))./abs(Swz_aggr(1,:))*100;

figure('Color','w');
subplot(211)
semilogx(f,20*log10(abs(S_PID)),'Color',cPID);
hold on;
semilogx(f,20*log10(pseudoMag_aggr),'Color',cReset);
semilogx(f,20*log10(abs(Swz_aggr(1,:))),'--','Color',cShape);
ylabel('Magnitude (dB)');
legend('$|S_{\mathrm{PID}}|$','$|S_{\infty}|$','$|S_{r,e}^{1}|$', ...
    'Location','best');
xlim([0.1 2000]);
grid on;

subplot(212)
semilogx(f,20*log10(abs(Swz_aggr(3,:))),'--','Color',cThird);
hold on;
semilogx(f,20*log10(abs(Swz_aggr(1,:))),'Color',cShape);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('$|S_{r,e}^{3}|$','$|S_{r,e}^{1}|$','Location','best');
xlim([0.1 2000]);
grid on;

figure('Color','w');
semilogx(f,sigma_2,'Color',cPID);
xlabel('Frequency (Hz)');
ylabel('$\sigma_2$ (\%)');
xlim([0.1 2000]);
grid on;

%% 5.1 Nonlinearity shaping using Psi
% The desired upper bound for sigma_2 is written as a ratio. For example,
% sigma_2_max = 0.15 corresponds to 15 percent.
sigma_2_max = 0.15;

% Psi gives the frequency-dependent bound used for the complementary
% pre/post-filter design.
sum_higher = sum(abs(Swz_aggr(3:end,:)).^2,1,'omitnan');
Psi = abs(Swz_aggr(1,:)).*sqrt( ...
    (sigma_2_max^2 + 2*sigma_2_max)./sum_higher);
Psi(sum_higher==0) = Inf;

% Notch filter used in the Chapter 4 example of the thesis:
%
%            s^2/w_n^2 + s/(Q_1 w_n) + 1
%   F(s) = ---------------------------------
%            s^2/w_n^2 + s/(Q_2 w_n) + 1
%
% The parameters below are used for the present toolbox example.
Q1 = 3;
Q2 = 2;
w_n = 2*pi*27.5;

F = (s^2/w_n^2 + s/(Q1*w_n) + 1) / ...
    (s^2/w_n^2 + s/(Q2*w_n) + 1);
F_inv = 1/F;

frfF = reshape(squeeze(freqresp(F,w)),1,[]);

% Approximate max_{odd n>=3}|F^{-1}(jnw)| using the odd harmonics up to 41.
Finv_max = zeros(size(f));
for n = 3:2:41
    frfFinv_n = reshape(squeeze(freqresp(F_inv,n*w)),1,[]);
    Finv_max = max(Finv_max,abs(frfFinv_n));
end

filter_bound = abs(frfF).*Finv_max;

figure('Color','w');

semilogx(f,20*log10(abs(Psi)),'Color',cPID);
hold on;
semilogx(f,20*log10(abs(filter_bound)),'Color',cReset);

xlabel('Frequency (Hz)','Interpreter','latex');
ylabel('Magnitude (dB)','Interpreter','latex');

legend('$\Psi(\omega)$', ...
    '$|F(j\omega)|\max_{n\geq3}|F^{-1}(jn\omega)|$', ...
    'Interpreter','latex','Location','best');

xlim([0.1 2000]);
grid on;

% Apply F before the reset element and F^{-1} after the reset element. The
% complementary placement preserves the first-order loop while shaping the
% higher-order harmonic contribution.
frfCpre_F = reshape(squeeze(freqresp(F,w)),1,[]);
frfCpar_F = D_r_aggr*ones(size(f));
frfCpos_F = reshape(squeeze(freqresp(C_pos_aggr*F_inv,w)),1,[]);

[~,~,Swz_F,~,~] = computePseudoSens_SF( ...
    f,sysR_aggr,A_rho,frfCs_aggr, ...
    frfCpre_F,frfCpar_F,frfCpos_F,frfPlant, ...
    Fs,nrHOSIDFsMax,overSampleRatio);

sigma_2_F = (sqrt(sum(abs(Swz_F).^2,1,'omitnan')) ...
    - abs(Swz_F(1,:)))./abs(Swz_F(1,:))*100;

figure('Color','w');
semilogx(f,sigma_2,'Color',cPID);
hold on;
semilogx(f,sigma_2_F,'Color',cReset);
yline(100*sigma_2_max,':','Color',cGray);
xlabel('Frequency (Hz)');
ylabel('$\sigma_2$ (\%)');
legend('Before shaping','After shaping','$\sigma_{2,\max}$', ...
    'Location','best');
xlim([0.1 2000]);
grid on;

%% 6. Reliability metrics sigma_t and sigma_d
% The reliability metrics are first evaluated with C_s = 1. The same
% closed-loop structure is then evaluated after designing C_s with
% design_Cs_maxPhase.

[~,nrResets_0,~,~,~,~,sigma_d_0,sigma_t_0] = ...
    computePseudoSens_SF_withMetrics( ...
    f,sysR_aggr,A_rho,frfCs_aggr, ...
    frfCpre_F,frfCpar_F,frfCpos_F,frfPlant, ...
    Fs,nrHOSIDFsMax,overSampleRatio);

%% 6.1 Shaping-filter design
% The shaping filter C_s is designed to improve the reliability of the
% reset action. The filter is designed subject to the shaping-filter
% constraints introduced in the thesis. The design parameters specify the
% filter order, high-frequency gain, desired phase-peak frequency, and
% maximum phase contribution.

p = 2;                         % Shaping-filter order
k_inf = 0.30;                  % Desired high-frequency gain
w_pk = w_BW;                   % Desired phase-peak frequency [rad/s]
phi_max = 15;                  % Target maximum phase [deg]

[Cs_rel,~,~,info] = design_Cs_maxPhase( ...
    k_inf,w_r_aggr,p,w_pk,phi_max);

% Frequency response of the designed shaping filter
frfCs_rel = reshape(squeeze(freqresp(Cs_rel,w)),1,[]);


%% Shaping filter and reset base-linear system
% The frequency response of C_s is shown together with the base-linear
% dynamics of the reset element. This illustrates the frequency range in
% which the shaping filter modifies the signal entering the reset element.

R_bl = ss(sysR_aggr.A,sysR_aggr.B,sysR_aggr.C,0);

frfRbl = reshape(squeeze(freqresp(R_bl,w)),1,[]);

figure('Color','w');

subplot(211)
semilogx(f,20*log10(abs(frfRbl)), ...
    'Color',cPID,'LineWidth',1.5);
hold on;
semilogx(f,20*log10(abs(frfCs_rel)), ...
    'Color',cReset,'LineWidth',1.5);

ylabel('Magnitude (dB)','Interpreter','latex');

legend('$R_{\mathrm{bl}}$', '$C_s$', ...
    'Interpreter','latex', ...
    'Location','best');

xlim([0.1 2000]);
grid on;


subplot(212)

phase_Rbl = rad2deg(unwrap(angle(frfRbl)));
phase_Cs  = rad2deg(unwrap(angle(frfCs_rel)));

semilogx(f,phase_Rbl, ...
    'Color',cPID,'LineWidth',1.5);
hold on;
semilogx(f,phase_Cs, ...
    'Color',cReset,'LineWidth',1.5);

xlabel('Frequency (Hz)','Interpreter','latex');
ylabel('Phase (degree)','Interpreter','latex');

xlim([0.1 2000]);
grid on;


%% Recalculate reliability metrics with the designed shaping filter

[~,nrResets_1,~,~,~,~,sigma_d_1,sigma_t_1] = ...
    computePseudoSens_SF_withMetrics( ...
    f,sysR_aggr,A_rho,frfCs_rel, ...
    frfCpre_F,frfCpar_F,frfCpos_F,frfPlant, ...
    Fs,nrHOSIDFsMax,overSampleRatio);

fprintf('\nReliability shaping filter\n');
fprintf('--------------------------\n');
fprintf('Order                  : %d\n',p);
fprintf('Phase-peak frequency   : %.2f Hz\n',w_pk/(2*pi));
fprintf('Achieved phase         : %.2f deg\n',info.bestPhase_deg);
fprintf('High-frequency gain    : %.4f\n',info.k_inf_achieved);


%% Reliability metrics before and after shaping-filter design
% Compare the number of reset events and the two reliability metrics before
% and after redesigning the shaping filter C_s.

figure('Color','w');

subplot(311)
semilogx(f,nrResets_0, ...
    'Color',cPID,'LineWidth',1.5);
hold on;
semilogx(f,nrResets_1, ...
    'Color',cReset,'LineWidth',1.5);

yline(2,':','Color',cGray);

ylabel('$N_r$', ...
    'Interpreter','latex');

legend('Before $C_s$ design','After $C_s$ design', ...
    'Interpreter','latex', ...
    'Location','best');

xlim([0.1 2000]);
grid on;


subplot(312)
semilogx(f,sigma_t_0, ...
    'Color',cPID,'LineWidth',1.5);
hold on;
semilogx(f,sigma_t_1, ...
    'Color',cReset,'LineWidth',1.5);

ylabel('$\sigma_t$ (\%)', ...
    'Interpreter','latex');

xlim([0.1 2000]);
grid on;


subplot(313)
semilogx(f,sigma_d_0, ...
    'Color',cPID,'LineWidth',1.5);
hold on;
semilogx(f,sigma_d_1, ...
    'Color',cReset,'LineWidth',1.5);

yline(100,':','Color',cGray);

xlabel('Frequency (Hz)', ...
    'Interpreter','latex');

ylabel('$\sigma_d$ (\%)', ...
    'Interpreter','latex');

xlim([0.1 2000]);
grid on;
%% 6.2 Effect of C_s on the reset HOSIDFs
% Compare the first- and third-order reset HOSIDFs before and after the
% shaping-filter design.
H1_before = computeResetHOSIDF_SF( ...
    sysR_aggr.A,sysR_aggr.B,sysR_aggr.C,D_r_aggr, ...
    A_rho,f,1,frfCs_aggr);
H3_before = computeResetHOSIDF_SF( ...
    sysR_aggr.A,sysR_aggr.B,sysR_aggr.C,D_r_aggr, ...
    A_rho,f,3,frfCs_aggr);

H1_after = computeResetHOSIDF_SF( ...
    sysR_aggr.A,sysR_aggr.B,sysR_aggr.C,D_r_aggr, ...
    A_rho,f,1,frfCs_rel);
H3_after = computeResetHOSIDF_SF( ...
    sysR_aggr.A,sysR_aggr.B,sysR_aggr.C,D_r_aggr, ...
    A_rho,f,3,frfCs_rel);

phase_H1_before = rad2deg(unwrap(angle(H1_before)));
phase_H1_after  = rad2deg(unwrap(angle(H1_after)));
phase_H1_before = phase_H1_before - 360*round(phase_H1_before(1)/360);
phase_H1_after  = phase_H1_after  - 360*round(phase_H1_after(1)/360);

figure('Color','w');
subplot(211)
semilogx(f,20*log10(abs(H1_before)),'Color',cPID);
hold on;
semilogx(f,20*log10(abs(H3_before)),'--','Color',cThird);
semilogx(f,20*log10(abs(H1_after)),'Color',cReset);
semilogx(f,20*log10(abs(H3_after)),'--','Color',cShape);
ylabel('Magnitude (dB)');
legend('$H_1$ before','$H_3$ before','$H_1$ after','$H_3$ after', ...
    'Location','best');
xlim([0.1 2000]);
grid on;

subplot(212)
semilogx(f,phase_H1_before,'Color',cPID);
hold on;
semilogx(f,phase_H1_after,'Color',cReset);
xlabel('Frequency (Hz)');
ylabel('Phase (degree)');
xlim([0.1 2000]);
grid on;
