function [pseudoMag, nrResets, Swz, Swy, Ln, sigma_2, sigma_d, sigma_t] = ...
    computePseudoSens_SF_withMetrics( ...
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
% Input definition:
%
% freqs - 1-by-n linearly-spaced frequency array [Hz], such that
%         freqs = [f_1, f_2, ..., f_n] with f_k = k*f_1
%
% sysR  - struct containing variables A_R, B_R, and C_R of system R
%         as in (3) of [2]
%
% A_rho - variable A_rho of system R as in (3) of [2]
%
% frfCpre - 1-by-n complex-valued array with frequency-response function
%           (FRF) of SISO LTI controller Cpre at frequencies 'freqs'
%
% frfCpar - 1-by-n complex-valued array with FRF of SISO LTI controller
%           Cpar at frequencies 'freqs'
%
% frfCpos - 1-by-n complex-valued array with FRF of SISO LTI controller
%           Cpos at frequencies 'freqs'
%
% frfCs   - 1-by-n complex-valued array with FRF of SISO LTI shaping
%           filter Cs at frequencies 'freqs'
%
% frfPlant - 1-by-n complex-valued array with FRF of SISO LTI plant P
%            at frequencies 'freqs'
%
% Fs - sampling frequency [Hz]
%
% nrHOSIDFsMax (optional) - Largest HOSIDF that should be taken into
%                           account
%
% overSampleRatioUser (optional) - Factor by which 'Fs' is increased
%
% ---------------------------------------------------------------------
% Output definition:
%
% pseudoMag - 1-by-n array with pseudo-sensitivity magnitudes at
%             frequencies 'freqs'
%
% nrResets  - 1-by-n array with number of resets per period of w at
%             frequencies 'freqs'
%
% Swz       - n-by-n (or nrHOSIDFsMax-by-n) array with higher-order
%             sinusoidal-input sensitivity functions (HOSISFs) at
%             frequencies 'freqs' (from w to z)
%
% Swy       - n-by-n (or nrHOSIDFsMax-by-n) array with HOSISFs at
%             frequencies 'freqs' (from w to y)
%
% Ln        - n-by-n (or nrHOSIDFsMax-by-n) array with HOSISFs at
%             frequencies 'freqs' (from U_in to U_out)
%
% sigma_2   - 1-by-n array containing the robustness metric sigma_2
%             (in percent), as defined in [2]
%
% sigma_d   - 1-by-n array containing the reset signal value reliability
%             metric sigma_d (in percent), as defined in [3]
%
% sigma_t   - 1-by-n array containing the reset event time reliability
%             metric sigma_t (in percent), as defined in [3]
%
% [2] Hosseini, Ali, Dragan Kostić, and Hassan HosseinNia.
%     "Robust performance analysis and nonlinearity shaping for
%     closed-loop reset control systems."
%     IEEE Transactions on Control Systems Technology (2026).
%
% [3] Hosseini, Ali, Dragan Kostić, and Hassan HosseinNia.
%     "Reliability Assessment and Performance Enhancement of Reset
%     Control Systems." arXiv:2606.21431 (2026).
%
% ---------------------------------------------------------------------


%% Input handling

nrFreqs = length(freqs);

if nargin == 9

    % Maximum number of HOSIDFs is equal to the number of frequencies
    nrHOSIDFs = nrFreqs;

    % Due to sampling, the maximum error might occur between two samples.
    % To approximate the worst-case pseudo-sensitivity, sample at a
    % 100-times higher rate than the sampling frequency.
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

% Compute LTI element G in Lure form
[frfGwz, frfGuz, frfGwy, frfGuy] = ...
    convertToLure_SF(frfCpre, frfCpar, frfCpos, frfPlant);

% Reset-element state-space matrices
A_R = sysR.A;
B_R = sysR.B;
C_R = sysR.C;


%% Compute HOSIDFs of the reset element

H_n = NaN(nrHOSIDFs, nrFreqs);

for nn = 1:nrHOSIDFs

    H_n(nn,:) = computeResetHOSIDF_SF( ...
        A_R, B_R, C_R, 0, A_rho, freqs, nn, frfCs);

end


% FRF of the reset-element base-linear system (BLS), i.e. A_rho = I_m
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
        Swy(1,:) = ...
            frfGwy ./ (1 - frfGuy .* H_n(1,:));               % (16) in [2]

        Swz(1,:) = ...
            frfGwz + frfGuz .* H_n(1,:) .* Swy(1,:);          % (19) in [2]

        Ln(1,:) = ...
            frfCpre .* (H_n(1,:) + frfCpar) .* ...
            frfCpos .* frfPlant;

        frfCsnn(1,:) = frfCs;

    else

        % Valid input-frequency indices for harmonic nn
        omegaIdxs = 1:floor(nrFreqs / nn);

        % Corresponding higher-harmonic frequency indices
        nnOmegaIdxs = nn * omegaIdxs;

        % Auxiliary frequency-response term
        frfDummy = ...
            H_n(nn,omegaIdxs) .* ...
            Swy(1,omegaIdxs) .* ...
            exp(1i*(nn-1).*angle(Swy(1,omegaIdxs))) ./ ...
            (1 - frfGuy(nnOmegaIdxs) .* frfRbl(nnOmegaIdxs));

        Swy(nn,omegaIdxs) = ...
            frfGuy(nnOmegaIdxs) .* frfDummy;                   % (17) in [2]

        Swz(nn,omegaIdxs) = ...
            frfGuz(nnOmegaIdxs) .* frfDummy;                   % (20) in [2]

        Ln(nn,omegaIdxs) = ...
            frfCpre(omegaIdxs) .* ...
            H_n(nn,omegaIdxs) .* ...
            frfCpos(nnOmegaIdxs) .* ...
            frfPlant(nnOmegaIdxs);

        frfCsnn(nn,omegaIdxs) = frfCs(nnOmegaIdxs);

    end

end


%% Compute pseudo-sensitivity and reliability metrics

% Sample time with oversampling
Ts = 1 / (Fs * overSampleRatio);

% Preallocate outputs and intermediate metric variables
pseudoMag      = zeros(1, nrFreqs);
nrResets       = zeros(1, nrFreqs);
Sigma_tr       = zeros(1, nrFreqs);
delta_distance = zeros(1, nrFreqs);


for ff = 1:nrFreqs

    %% Reconstruct performance output and reset-element input

    freqInput = freqs(ff);

    % Period corresponding to the current input frequency
    Tperiod = 1 / freqInput;

    % Time vector over one period
    time = 0:Ts:Tperiod;

    % Largest HOSIDF available for this input frequency
    finalHOSIDF = min( ...
        nrHOSIDFs, ...
        floor(freqs(end) / freqInput));

    nrSamples = length(time);

    perfOutput  = zeros(1, nrSamples);
    resetInput  = zeros(1, nrSamples);
    resetInput_1 = zeros(1, nrSamples);


    for nn = 1:finalHOSIDF

        % Performance output based on (18) in [2]
        % for \hat{w} = 1 and \varphi_w = 0
        perfOutput = perfOutput + ...
            abs(Swz(nn,ff)) .* ...
            sin( ...
                nn*2*pi*freqInput*time + ...
                angle(Swz(nn,ff)));

        % Reset-element input based on (15) in [2]
        resetFRF = frfCsnn(nn,ff) .* Swy(nn,ff);

        resetInput = resetInput + ...
            abs(resetFRF) .* ...
            sin( ...
                nn*2*pi*freqInput*time + ...
                angle(resetFRF));

    end


    % First-harmonic approximation of the reset-triggering signal
    resetFRF_1 = frfCsnn(1,ff) .* Swy(1,ff);

    resetInput_1 = ...
        abs(resetFRF_1) .* ...
        sin( ...
            2*pi*freqInput*time + ...
            angle(resetFRF_1));


    %% Pseudo-sensitivity magnitude

    pseudoMag(ff) = max(abs(perfOutput));


    %% Reliability metrics for reset-triggering signal

    % q  = actual reset-triggering signal
    % q1 = first-harmonic approximation
    q  = resetInput(:).';
    q1 = resetInput_1(:).';

    N = numel(time);


    %% Find one nominal zero crossing of q1(t)

    [~, k1_max_q1] = max(q1);
    [~, k1_min_q1] = min(q1);

    idx_min_q1 = min(k1_min_q1, k1_max_q1);
    idx_max_q1 = max(k1_min_q1, k1_max_q1);

    % Zero crossings of q1(t)
    nr_1_all = find(q1 .* circshift(q1,1) <= 0);

    % Keep the nominal crossing between the maximum and minimum of q1(t)
    nr_1 = nr_1_all( ...
        nr_1_all > idx_min_q1 & ...
        nr_1_all < idx_max_q1);

    % Safety check
    if isempty(nr_1)

        Sigma_tr(ff)       = NaN;
        delta_distance(ff) = NaN;

        return

    end

    % If multiple neighboring indices are detected because of an exact
    % numerical zero, keep one representative crossing.
    nr_1 = nr_1(1);


    %% Find actual zero crossings of q(t) in the same half-period window

    nr_ss_all = find(q .* circshift(q,1) <= 0);

    nr_ss = nr_ss_all( ...
        nr_ss_all > idx_min_q1 & ...
        nr_ss_all < idx_max_q1);

    % Remove duplicate neighboring detections caused by exact zeros
    if ~isempty(nr_ss)

        nr_ss = nr_ss([true, diff(nr_ss) > 1]);

    end


    if isempty(nr_ss)

        Sigma_tr(ff)       = NaN;
        delta_distance(ff) = NaN;

        return

    end


    %% sigma_t: reset-time deviation metric

    % For Nr = 2, this represents the shift between the nominal crossing
    % nr_1 and the actual crossing.
    %
    % For Nr > 2, it represents the spread between nr_1 and all local
    % actual crossings.
    delta_sample = ...
        max([nr_1, nr_ss]) - min([nr_1, nr_ss]);

    % Normalize with respect to half a period
    Sigma_tr(ff) = delta_sample / (N/2) * 100;


    %% sigma_d: zero-crossing reliability metric

    q_inf = max(abs(q));

    if q_inf == 0

        delta_distance(ff) = NaN;

    else

        % Derivatives with respect to sample index.
        % Scaling with the sample time is unnecessary for the sign
        % conditions used below.
        dq  = gradient(q);
        ddq = gradient(dq);
        dq1 = gradient(q1);


        if numel(nr_ss) > 1

            %% Case Nr > 2

            % Search only between the nominal crossing and the group of
            % actual crossings.
            idx_I_min = min([nr_1, nr_ss]);
            idx_I_max = max([nr_1, nr_ss]);

            idx_I_min = max(idx_I_min, 2);
            idx_I_max = min(idx_I_max, N-1);

            % Stationary points of q(t)
            k_d_all = find(dq .* circshift(dq,1) <= 0);

            k_d_all = k_d_all( ...
                k_d_all >= idx_I_min & ...
                k_d_all <= idx_I_max);

            % Remove duplicate neighboring detections
            if ~isempty(k_d_all)

                k_d_all = k_d_all([true, diff(k_d_all) > 1]);

            end

            % Conditions for E_{>2,k}:
            %
            %   dot(q) = 0
            %   q*ddot(q) < 0
            %   dot(q1)*ddot(q) > 0
            valid_idx = ...
                q(k_d_all) .* ddq(k_d_all) < 0 & ...
                dq1(k_d_all) .* ddq(k_d_all) > 0;

            k_d_0 = k_d_all(valid_idx);

            if ~isempty(k_d_0)

                Delta_d = max(abs(q(k_d_0)));

                delta_distance(ff) = ...
                    (1 + Delta_d/q_inf) * 100;

            else

                % Metric cannot be evaluated reliably when no valid
                % stationary point is detected.
                delta_distance(ff) = NaN;

            end


        else

            %% Case Nr = 2

            % Search for critical stationary points between the maximum
            % and minimum of q(t).
            [~, k_max_q] = max(q);
            [~, k_min_q] = min(q);

            idx_min_q = min(k_min_q, k_max_q);
            idx_max_q = max(k_min_q, k_max_q);

            idx_min_q = max(idx_min_q, 2);
            idx_max_q = min(idx_max_q, N-1);

            % Stationary points of q(t)
            k_d_all = find(dq .* circshift(dq,1) <= 0);

            k_d_all = k_d_all( ...
                k_d_all >= idx_min_q & ...
                k_d_all <= idx_max_q);

            % Remove duplicate neighboring detections
            if ~isempty(k_d_all)

                k_d_all = k_d_all([true, diff(k_d_all) > 1]);

            end

            % Conditions for E_2:
            %
            %   dot(q) = 0
            %   q*ddot(q) > 0
            valid_idx = ...
                q(k_d_all) .* ddq(k_d_all) > 0;

            k_d_0 = k_d_all(valid_idx);

            if ~isempty(k_d_0)

                Delta_d = min(abs(q(k_d_0)));

                delta_distance(ff) = ...
                    (1 - Delta_d/q_inf) * 100;

            else

                % No critical stationary point corresponds to the most
                % reliable case:
                %
                % Delta_d = ||q||_inf  -->  sigma_d = 0
                delta_distance(ff) = 0;

            end

        end

    end


    %% Compute number of reset instants per period

    [~, counted_crossings_with_half] = zerocrossrate(resetInput);

    % The first nonzero value is counted as half a crossing by
    % zerocrossrate, so take the integer part.
    counted_crossings_unchecked = ...
        floor(counted_crossings_with_half);

    % The number of crossings over one complete period should be even.
    % If an odd number is obtained, one crossing was missed at the
    % boundary of the sampled period.
    if mod(counted_crossings_unchecked, 2) == 1

        counted_crossings = ...
            counted_crossings_unchecked + 1;

    else

        counted_crossings = ...
            counted_crossings_unchecked;

    end

    nrResets(ff) = counted_crossings;

end


%% Metric outputs

sigma_2 = ...
    (sqrt(nansum(abs(Swz).^2)) - abs(Swz(1,:))) ./ ...
    abs(Swz(1,:)) * 100;

sigma_d = delta_distance;
sigma_t = Sigma_tr;

end