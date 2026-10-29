function [Hn_sf] = computeResetHOSIDF_SF(A_R,B_R,C_R,D_R,A_rho,freqs,n,frfCs)
    % ---------------------------------------------------------------------
    % This function computes a higher-order sinusoidal input describing
    % function (HOSIDF) of a reset controller as in (1) of [1]
    % v00 - Ali Hosseini (10-09-2026)
    % Code based on:
    % v03 - Luke van Eijk (17/06/2024)
    % [1] Karbasizadeh & HosseinNia, "Complex-order Reset Control System"
    % ---------------------------------------------------------------------
    % Input definition:
    % A_R, B_R, C_R, D_R, A_rho - variables of reset controller as in (1) of [1]
    % freqs - 1-by-m frequency vector [Hz]
    % n     - HOSIDF order [-]
    % frfCs - 1-by-n complex-valued array with FRF of SISO LTI shaping
    %           filter Cs at frequencies 'freqs'
    % Output definition:
    % Hn    - 1-by-m complex-valued array with the reset controller's 
    %           n^th-order HOSIDF, as in (2) of [1]
    % ---------------------------------------------------------------------

    if min(freqs) <= 0
        error('Only positive frequencies are allowed')
    end
    comR = zeros(1,length(freqs));
    Hn_sf = zeros(1,length(freqs));
    % Formulas below taken from (2) in [1]
    nrFreqs = length(freqs);
    if mod(n,2) == 0 % Even orders
        Hn = zeros(1,nrFreqs);
    elseif mod(n,2) == 1 % Odd orders
        Hn = zeros(1,nrFreqs);

        omega = 2*pi*freqs;
        nrStates = length(B_R);   
        In = eye(nrStates); On = zeros(nrStates);
        for ff = 1:nrFreqs
            Lambda = omega(ff)^2 * In + A_R^2;
            Delta = In + expm(pi*A_R/omega(ff));
            Delta_rho = In + A_rho * expm(pi*A_R/omega(ff));

            Lambda = omega(ff)*omega(ff)*eye(size(A_R)) + A_R^2;
            Delta = eye(size(A_R)) + expm(A_R*pi/omega(ff));
            Delta_rho = eye(size(A_R)) + A_rho*expm(A_R*pi/omega(ff));

            %Calculate the Phase of shaping filter
            phase_Cs=phase(frfCs(ff));
            phi = -1*phase_Cs;
            Zeta0 = (omega(ff)*cos(phi)*eye(size(A_R)) + A_R*sin(phi)) / Lambda * B_R; 
            Theta_phi= ( -2*1i*omega(ff)*exp(-1i*phi) / pi) * Delta * ( eye(size(A_R)) - Delta_rho \ (A_rho*Delta) ) * Zeta0;



            % Calculate Delta_r & Gamma_r
            % if A_rho == On
            %     Delta_r = In;
            %     Gamma_r = On;
            % else
                
                Gamma = inv(Delta_rho) * A_rho * Delta * inv(Lambda);
            % end
            
            Theta_D = -2*omega(ff)^2 / pi * Delta * (Gamma - inv(Lambda));
            if n == 1
                Hn(ff) = C_R * inv(1i*omega(ff)*In - A_R) * (In + 1i*Theta_D) * B_R + D_R;
                comR(ff) = C_R / (A_R - 1i*omega(ff)*eye(size(A_R))) * Theta_phi - C_R * (1i*omega(ff)*eye(size(A_R)) + A_R)/Lambda * B_R;
                comR(ff) = comR(ff) + D_R;
            else
                Theta_phi= ( -2*1i*omega(ff)*exp(-1i*phi*n) / pi) * Delta * ( eye(size(A_R)) - Delta_rho \ (A_rho*Delta) ) * Zeta0;
                Hn(ff) = C_R * inv(1i*n*omega(ff)*In - A_R) * 1i*Theta_D * B_R;
                comR(ff) = C_R / (A_R - 1i*omega(ff)*n*eye(size(A_R))) * Theta_phi; 
            end
            Hn_sf(ff)=comR(ff);
        end
    else
        error('Only natural numbers are allowed for the HOSIDF order')
    end
end

