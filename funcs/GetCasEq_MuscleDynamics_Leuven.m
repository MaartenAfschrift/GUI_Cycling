function [CasF] = GetCasEq_MuscleDynamics_Leuven(parms)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

import casadi.*

% create some generic casadi variables used for the functions
nseg = 4; % number of segments in the model
nmus = 9;
lcerel = SX.sym('lcerel', nmus);
a      = SX.sym('a',nmus,1);
fi     = SX.sym('fi', nseg,1);
stim   = SX.sym('stim',nmus,1);
vcerel_help = SX.sym('fcerel_help',nmus,1);
vcerel_in = SX.sym('vcerel_in', nmus);

% function for smooth transition
sigma = @(x,x0,w) 1./(1+exp(-(x-x0)/w));

%% Muscle geometry (this should be in another function)

fijo = [fi(3)-fi(2) fi(4)-fi(3) 1.1-fi(4)];
lmtc = zeros(1,length(parms.nmus)) + parms.A0(:,1)' + parms.A1(:,1)'.*fijo(1) + parms.A2(:,1)'.*fijo(1).^2;
lmtc = lmtc +  parms.A0(:,2)' + parms.A1(:,2)'.*fijo(2) + parms.A2(:,2)'.*fijo(2).^2;
lmtc = lmtc +  parms.A0(:,3)' + parms.A1(:,3)'.*fijo(3) + parms.A2(:,3)'.*fijo(3).^2;
momarm = parms.A1 + 2*repmat(fijo,parms.nmus,1).*parms.A2;
momarm = -momarm;
CasF.f_lMT_dM = Function('f_lMT_dM',{fi},{lmtc,momarm},...
    {'fi'},{'lmtc','momarm'});

%% simple activation dynamics
a_dot = (stim-a)/0.03;
CasF.f_adot = Function('f_adot',{a,stim},{a_dot});

%% Force velocity model Leuven

[Fpe,FMltilde,FMvtilde] = getForceLengthVelocityProperties(lcerel,vcerel_help,10);

CasF.f_FLV_Leuven = Function('f_FLV_Leuven',{lcerel,vcerel_help},...
    {Fpe,FMltilde,FMvtilde},{'lcerel','vcerel'},...
    {'Fpe','FMltilde','FMvtilde'});

%% FL tendon leuven
lse=lmtc'-lcerel.*parms.lce_opt'; %tendon length
lTtilde = lse./parms.lse0';
fse = (exp(35.*(lTtilde - 0.995)))/5-0.25+0;
fse_N= fse.*parms.fmax';
CasF.f_FL_Tendon_Leuven = Function('f_FL_Tendon_Leuven',{lcerel,fi},...
    {fse_N},{'lcerel','fi'},{'fse'});


end