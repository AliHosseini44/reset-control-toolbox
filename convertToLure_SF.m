function [frfGwz, frfGuz, frfGwy, frfGuy] = convertToLure_SF(frfCpre, frfCpar, frfCpos, frfPlant)
    % ---------------------------------------------------------------------
    % This function converts the reset control system below into the Lure
    % form (as in Fig. 1 of [1]) with reference 'r' as external input and 
    % error 'e' as performance output
    %
    %          --> C2 --> R --> C3 --
    %          |                    |                
    % --> C1 --|                    + --> C5 --> P -->
    %          |                    |                
    %          ---------> C4 --------
    %
    % v03 - Luke van Eijk (17/06/2024)
    % Code based on:
    % [1] L. F. van Eijk, "Loop-shaping for reset control systems with
    %       linear control elements," v00, internal report.
    % ---------------------------------------------------------------------
    % Input definition:
    % frfC1 - 1-by-n complex-valued array with frequency-response function
    %               (FRF) of SISO LTI controller C1
    % frfC2 - 1-by-n complex-valued array with FRF of SISO LTI controller C2
    % frfC3 - 1-by-n complex-valued array with FRF of SISO LTI controller C3
    % frfC4 - 1-by-n complex-valued array with FRF of SISO LTI controller C4
    % frfC5 - 1-by-n complex-valued array with FRF of SISO LTI controller C5
    % frfPlant - 1-by-n complex-valued array with FRF of SISO LTI plant P
    %
    % Output definition:
    % frfGwz    - 1-by-n complex-valued array with FRF of SISO LTI element
    %               from external input 'w' to performance output 'z'
    % frfGuz    - 1-by-n complex-valued array with FRF of SISO LTI element
    %               from 'u' to 'z'
    % frfGwy    - 1-by-n complex-valued array with FRF of SISO LTI element
    %               from 'w' to 'y'
    % frfGuy    - 1-by-n complex-valued array with FRF of SISO LTI element
    %               from 'u' to 'y'
    % ---------------------------------------------------------------------

    %% Convert to form in Fig. 4 of [1]
    %            ---> R -----
    %            |          |                
    % --> Cpre --|          + --> Cpos --> P -->
    %            |          |                
    %            --> Cpar ---
    % 
    % frfCpre = frfC1 .* frfC2;
    % frfCpar = frfC4 ./ frfC2 ./ frfC3;
    % frfCpos = frfC3 .* frfC5;


    %% Convert to Lure form in Fig. 1 of [1]
    frfGwz = 1 ./ (1 + frfPlant .* frfCpos .* frfCpar .* frfCpre);  % (39) in [1]
    frfGuz = -frfPlant .* frfCpos .* frfGwz;                        % (40) in [1]
    frfGwy = frfCpre .* frfGwz;                                     % (41) in [1]
    frfGuy = frfCpre .* frfGuz;                                     % (42) in [1]
end