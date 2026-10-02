function [pseudoMag, nrResets, Swz, Swy, Ln] = computePseudoSens_SF( ...
    freqs, sysR, A_rho, frfCs, frfCpre, frfCpar, frfCpos, ...
    frfPlant, Fs, nrHOSIDFsMax, overSampleRatioUser)

% ---------------------------------------------------------------------
% This function can be utilized to compute the pseudo-sensitivity
% magnitude (from w to z) of a general zero-crossing reset control system:
%
%             z   _____    w
%           <----|     |<-----
%             y  |  G  |   u
%           -----|_____|<-----
%           |     _____      |
%           |    |     |     |
%           |--->|  R  |-----|
%                |_____|
%
%               ---------> R --------
%   U_in        |          |         |               U_out
%    ---->Cpre--| ----Cs---          + --> Cpos --> P -->
%               |                    |
%               ---------> Cpar ------
%
% Signal z and the underlying higher-order sinusoidal-input sensitivity
% functions (HOSISFs) are computed using [1, Theorem 3.1]. The pseudo-
% sensitivity is computed based on the maximum of z, as in (37) of [1].
%
% v00 - Ali Hosseini  (10-09-2026)
%
% Code based on:
% v00 - Luke van Eijk (19/05/2025)
%
% [1] L.F. van Eijk, D. Kostić, S.H. HosseinNia,
%     "Frequency Response Analysis of General Zero-Crossing Reset
%     Control Systems," submitted to IEEE Control Systems Letters
%
% ---------------------------------------------------------------------
%
% Input definition:
%
% freqs - 1-by-n linearly-spaced frequency array [Hz], such that
%         freqs = [f_1, f_2, ..., f_n] with f_k = k*f_1
%
% sysR  - struct containing variables A_R, B_R, and C_R of system R
%
% A_rho - variable A_rho of system R
%
% frfCpre - 1-by-n complex-valued array with FRF of SISO LTI
%           controller Cpre at frequencies 'freqs'
%
% frfCpar - 1-by-n complex-valued array with FRF of SISO LTI
%           controller Cpar at frequencies 'freqs'
%
% frfCpos - 1-by-n complex-valued array with FRF of SISO LTI
%           controller Cpos at frequencies 'freqs'
%
% frfCs   - 1-by-n complex-valued array with FRF of SISO LTI
%           shaping filter Cs at frequencies 'freqs'
%
% frfPlant - 1-by-n complex-valued array with FRF of SISO LTI
%            plant P at frequencies 'freqs'
%
% Fs - sampling frequency [Hz]
%
% nrHOSIDFsMax (optional) - largest HOSIDF to take into account
%
% overSampleRatioUser (optional) - factor by which Fs is increased
%
% ---------------------------------------------------------------------
%
% Output definition:
%
% pseudoMag - 1-by-n array with pseudo-sensitivity magnitudes
%
% nrResets  - 1-by-n array with number of resets per period
%
% Swz       - n-by-n (or nrHOSIDFsMax-by-n) array with HOSISFs
%             from w to z
%
% Swy       - n-by-n (or nrHOSIDFsMax-by-n) array with HOSISFs
%             from w to y
%
% Ln        - n-by-n (or nrHOSIDFsMax-by-n) array with HOSISFs
%             from U_in to U_out
%
% ---------------------------------------------------------------------

%% Input handling

nrFreqs = length(freqs);

if nargin == 9

    nrHOSIDFs = nrFreqs;
    overSampleRatio = 100;

elseif nargin == 10

    nrHOSIDFs = min(nrHOSIDFsMax, nrFreqs);
    overSampleRatio = 100;

elseif nargin == 11

    nrHOSIDFs = min(nrHOSIDFsMax, nrFreqs);
    overSampleRatio = overSampleRatioUser;

else

    error('An unexpected error occurred.');

end


%% Compute HOSISFs

% Convert system to Lure form
[frfGwz, frfGuz, frfGwy, frfGuy] = ...
    convertToLure_SF(frfCpre, frfCpar, frfCpos, frfPlant);

% Reset-element state-space matrices
A_R = sysR.A;
B_R = sysR.B;
C_R = sysR.C;


%% Compute HOSIDFs of reset element

H_n = NaN(nrHOSIDFs, nrFreqs);

for nn = 1:nrHOSIDFs

    H_n(nn,:) = computeResetHOSIDF_SF( ...
        A_R, B_R, C_R, 0, A_rho, freqs, nn, frfCs);

end


% FRF of reset-element base-linear system (BLS)
frfRbl = computeResetHOSIDF_SF( ...
    A_R, B_R, C_R, 0, eye(length(B_R)), freqs, 1, frfCs);


%% Preallocate HOSISFs

Swy     = NaN(nrHOSIDFs, nrFreqs);
Swz     = NaN(nrHOSIDFs, nrFreqs);
Ln      = NaN(nrHOSIDFs, nrFreqs);
frfCsnn = NaN(nrHOSIDFs, nrFreqs);


%% Compute HOSISFs

for nn = 1:nrHOSIDFs

    if nn == 1

        % Fundamental harmonic
        Swy(1,:) = frfGwy ./ ...
            (1 - frfGuy .* H_n(1,:));

        Swz(1,:) = frfGwz + ...
            frfGuz .* H_n(1,:) .* Swy(1,:);

        Ln(1,:) = ...
            frfCpre .* ...
            (H_n(1,:) + frfCpar) .* ...
            frfCpos .* ...
            frfPlant;

        frfCsnn(1,:) = frfCs;

    else

        % Valid input-frequency indices for harmonic nn
        omegaIdxs = 1:floor(nrFreqs / nn);

        % Corresponding output-frequency indices
        nnOmegaIdxs = nn * omegaIdxs;

        % Auxiliary term
        frfDummy = ...
            H_n(nn,omegaIdxs) .* ...
            Swy(1,omegaIdxs) .* ...
            exp(1i*(nn-1).*angle(Swy(1,omegaIdxs))) ./ ...
            (1 - frfGuy(nnOmegaIdxs) .* frfRbl(nnOmegaIdxs));

        % HOSISFs
        Swy(nn,omegaIdxs) = ...
            frfGuy(nnOmegaIdxs) .* frfDummy;

        Swz(nn,omegaIdxs) = ...
            frfGuz(nnOmegaIdxs) .* frfDummy;

        Ln(nn,omegaIdxs) = ...
            frfCpre(omegaIdxs) .* ...
            H_n(nn,omegaIdxs) .* ...
            frfCpos(nnOmegaIdxs) .* ...
            frfPlant(nnOmegaIdxs);

        frfCsnn(nn,omegaIdxs) = ...
            frfCs(nnOmegaIdxs);

    end

end


%% Compute pseudo-sensitivity magnitude and number of resets

% Sample time with oversampling
Ts = 1 / (Fs * overSampleRatio);

pseudoMag = zeros(1, nrFreqs);
nrResets  = zeros(1, nrFreqs);


for ff = 1:nrFreqs

    %% Input frequency

    freqInput = freqs(ff);

    % Period of the input signal
    Tperiod = 1 / freqInput;

    % Time vector over one period
    time = 0:Ts:Tperiod;

    % Largest HOSIDF available at this frequency
    finalHOSIDF = min( ...
        nrHOSIDFs, ...
        floor(freqs(end) / freqInput));

    nrSamples = length(time);


    %% Reconstruct performance output and reset-element input

    perfOutput = zeros(1, nrSamples);
    resetInput = zeros(1, nrSamples);

    for nn = 1:finalHOSIDF

        % Performance output
        perfOutput = perfOutput + ...
            abs(Swz(nn,ff)) .* ...
            sin( ...
                nn*2*pi*freqInput*time + ...
                angle(Swz(nn,ff)));

        % Reset-element input
        resetFRF = frfCsnn(nn,ff) .* Swy(nn,ff);

        resetInput = resetInput + ...
            abs(resetFRF) .* ...
            sin( ...
                nn*2*pi*freqInput*time + ...
                angle(resetFRF));

    end


    %% Pseudo-sensitivity magnitude

    pseudoMag(ff) = max(abs(perfOutput));


    %% Number of reset instants per period

    [~, counted_crossings_with_half] = zerocrossrate(resetInput);

    % First nonzero value is counted as a half crossing by zerocrossrate.
    % Taking floor removes this contribution.
    counted_crossings_unchecked = ...
        floor(counted_crossings_with_half);

    % For a periodic signal, the number of zero crossings should be even.
    % If an odd number is detected, one crossing was missed at the period
    % boundary.
    if mod(counted_crossings_unchecked, 2) == 1

        counted_crossings = ...
            counted_crossings_unchecked + 1;

    else

        counted_crossings = ...
            counted_crossings_unchecked;

    end

    nrResets(ff) = counted_crossings;

end

end