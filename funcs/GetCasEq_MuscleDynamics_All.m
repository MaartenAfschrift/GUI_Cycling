function [CasF] = GetCasEq_MuscleDynamics_All(parms)
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
fcerel_help = SX.sym('fcerel_help',nmus,1);

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

%% force-length tendon (VU)
lse=lmtc'-lcerel.*parms.lce_opt'; %tendon length
ese=lse-parms.lse0'; % tendon elongation
sig_fse = sigma(ese,0,0.001);
fse = (1-sig_fse).*((parms.b_Fse/.1)*ese+parms.b_Fse) + sig_fse.*(parms.kse'.*ese.^2+parms.b_Fse);
CasF.f_FL_tendon = Function('f_FL_tendon',{lcerel,fi},{fse},...
    {'lcerel','fi'},{'fse'});

% fiber length assuming rigid tendon
lcerel_rigidTendon = (lmtc-parms.lse0)./parms.lce_opt;
CasF.f_lce_Rigid_tendon = Function('f_lce_Rigid_tendon',{fi},{lcerel_rigidTendon'},...
    {'LMT'},{'lcerel'});

% muscle force length relation
fisomrel = Muscle_F_L_VU(lcerel);
CasF.f_FL_muscle = Function('f_FL_muscle',{lcerel},{fisomrel},...
    {'lcerel'},{'fisomrel'});

% relative fiber force
fce=fse;
fcerel=fce./parms.fmax'; % actual relative muscle force
CasF.f_fcerel = Function('f_fcerel',{lcerel,fi},{fcerel},...
    {'lcerel','fi'},{'fcerel'});

%% force-velocity
fcerel = fcerel_help;

% bool_ActivationDependent = false; % F/v profile independent of activation
% bool_LengthDependent = true;
[vcerel] = Muscle_F_V_VU_v2(fcerel,a,lcerel,parms);
% function for state derivative muscle 
CasF.f_FV_muscle = Function('f_FV_muscle',{lcerel,a,fcerel_help},{vcerel},...
    {'lcerel','a','fcerel_help'},{'vcerel'});

%% very simple force velocity
fcerel_out = 0.1*vcerel_in+1;
CasF.f_FV_muscle_simple = Function('f_FV_muscle_simple',{vcerel_in},{fcerel_out},...
    {'vcerel'},{'fcerel'});

% idea is that 
% Fm = Fl*Fv*a
% Fm = Fl * (0.1*vcerel_in+1)*a;
% vc = 10* Fm/(Fl*a) - 10
%    --> main problem is that you divide by a (which can be very small)
%     => high velocities


end