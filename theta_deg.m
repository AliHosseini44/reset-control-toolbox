function theta=theta_deg(frfPlant,frfCpre,frfCpos,frfCpar,frfCs,CR,Dr)

    % closed-loop system with SISO open-loop as follows:
    %
    %           - -------^---> R -- -----
    %          | ----    |             |                
    % --> C1 --|    |_Cs_|             + --> C2 --> G -->
    %          |                       |               
    %           ------------> C3 -------
    % 

C1=frfCpre;
C3=frfCpar;
C2=frfCpos;
Cs=frfCs;
G=frfPlant;


%%
L=C1.*C2.*G;

M1=1+L.*(CR+C3);
M2=L.*Cs.*(CR-Dr);
M3=(1+L.*(C3+Dr)).*(CR-Dr);
 
Nx=real(conj(M1).*M2);
Ny=real(conj(M1).*M3);
theta=rad2deg(atan2(Ny,Nx));