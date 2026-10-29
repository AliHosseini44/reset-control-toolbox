function [omega_f] = calc_omega_f(f,w_BW,theta_cglp,omega_l,A_rho,frfCs)







Theta=theta_cglp;
w=w_BW;
kw=find(abs(2*pi*f-w)<1);
kw=min(kw);

%%

g=sqrt(1+16*(1-A_rho)^2/pi^2/(1+A_rho)^2);
wr_r=omega_l/g;
Ar=-wr_r;Br=1;Cr=wr_r;
sys=ss(Ar,Br,Cr,0);
%FORE1=hosidfcalc(sys,A_rho,1,2*pi*f);
FORE1=computeResetHOSIDF_SF(Ar,Br,Cr,0,A_rho,f,1,frfCs);

%%
wr=omega_l;
a=real(FORE1(kw));
b=imag(FORE1(kw));
Q=tan(deg2rad(Theta)-atan(w/wr));


k1=a*Q-b;
k2=b*w*Q+b*wr+a*w-(a-1)*wr*Q;
k3=-w*wr*(b*Q+a-1);

wf1=(-k2+sqrt(k2^2-4*k1*k3))/(2*k1);
wf2=(-k2-sqrt(k2^2-4*k1*k3))/(2*k1);

wf=max(wf1,wf2);

if wf>0 && phase(wf)==0
omega_f=wf;
else
    error('The requested phase exceeds the maximum achievable phase')
end


end

