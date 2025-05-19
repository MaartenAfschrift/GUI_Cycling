function [fisomrel] = Muscle_F_L_VU(lcerel)
%Muscle_F_L_VU Force-length relation muscle model VU
%   Detailed explanation goes here
Fmin=1e-2;
c1_gaus = 0.407742573856005;
c2_gaus = 0.447996880165858; % w = 0.56
fisomrel = ((1-Fmin)/2) *( exp(-((lcerel-1)./c1_gaus).^(2)) +...
    exp(-((lcerel-1)./c2_gaus).^(4))) +Fmin;
end