function [Fmetab] = CasF_MetabolicEnergy(parms)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


import casadi.*

%% get muscle properties

SpecficTension = parms.tension;
NMuscle = parms.nmus;
FMo = parms.fmax';
volM = FMo.*parms.lce_opt';
musclemass = volM.*(1059.7)./(SpecficTension*1e4);

%% Metabolic energy models: Barghava 2004
act_SX          = SX.sym('act_SX',NMuscle,1); % Muscle activations
exc_SX          = SX.sym('exc_SX',NMuscle,1); % Muscle excitations
lMtilde_SX      = SX.sym('lMtilde_SX',NMuscle,1); % N muscle fiber lengths
vM_SX           = SX.sym('vM_SX',NMuscle,1); % Muscle fiber velocities
Fce_SX          = SX.sym('FT_SX',NMuscle,1); % Contractile element forces
Fpass_SX        = SX.sym('FT_SX',NMuscle,1); % Passive element forces
Fiso_SX         = SX.sym('Fiso_SX',NMuscle,1); % N forces (F-L curve)
pctst_SX        = SX.sym('pctst_SX',NMuscle,1); % Slow twitch ratio
modelmass_SX    = SX.sym('modelmass_SX',1); % Model mass
b_SX            = SX.sym('b_SX',1); % Parameter determining tanh smoothness
ScaleRate       = SX.sym('ScaleRate',1); % scale heath rate components
% Bhargava et al. (2004)
[energy_total_sm_SX,Adot_sm_SX,Mdot_sm_SX,Sdot_sm_SX,Wdot_sm_SX,...
    energy_model_sm_SX] = getMetabolicEnergySmooth2004all(exc_SX,act_SX,...
    lMtilde_SX,vM_SX,Fce_SX,Fpass_SX,musclemass,pctst_SX,Fiso_SX,...
    FMo,modelmass_SX,b_SX,ScaleRate);
fgetMetabolicEnergySmooth2004all = ...
    Function('fgetMetabolicEnergySmooth2004all',...
    {exc_SX,act_SX,lMtilde_SX,vM_SX,Fce_SX,Fpass_SX,...
    pctst_SX,Fiso_SX,modelmass_SX,b_SX,ScaleRate},{energy_total_sm_SX,...
    Adot_sm_SX,Mdot_sm_SX,Sdot_sm_SX,Wdot_sm_SX,energy_model_sm_SX});

%% Umberger metabolic energy models
vMtilde = vM_SX./parms.lce_opt';

% 2003 energy model
[energy_total_SX,energy_am_SX,energy_sl_SX,energy_mech_SX,energy_model_SX] = ...
    getMetabolicEnergySmooth2003all(exc_SX,act_SX,lMtilde_SX,vMtilde,vM_SX,Fce_SX,...
        musclemass,pctst_SX,10,Fiso_SX,modelmass_SX,b_SX);
fgetMetabolicEnergySmooth2003all = ...
    Function('fgetMetabolicEnergySmooth2003all',...
    {exc_SX,act_SX,lMtilde_SX,vM_SX,Fce_SX,pctst_SX,Fiso_SX,modelmass_SX,b_SX},...
        {energy_total_SX,energy_am_SX,energy_sl_SX,energy_mech_SX,energy_model_SX});

% 2010 energy model
[energy_total_SX,energy_am_SX,energy_sl_SX,energy_mech_SX,energy_model_SX] = ...
    getMetabolicEnergySmooth2010all(exc_SX,act_SX,lMtilde_SX,vMtilde,vM_SX,Fce_SX,...
        musclemass,pctst_SX,10,Fiso_SX,modelmass_SX,b_SX);
fgetMetabolicEnergySmooth2010all = ...
    Function('fgetMetabolicEnergySmooth2010all',...
    {exc_SX,act_SX,lMtilde_SX,vM_SX,Fce_SX,pctst_SX,Fiso_SX,modelmass_SX,b_SX},...
        {energy_total_SX,energy_am_SX,energy_sl_SX,energy_mech_SX,energy_model_SX});

% 2016 energy model
[energy_total_SX,energy_am_SX,energy_sl_SX,energy_mech_SX,energy_model_SX] = ...
    getMetabolicEnergySmooth2016all(exc_SX,act_SX,lMtilde_SX,vMtilde,vM_SX,Fce_SX,...
        musclemass,pctst_SX,10,Fiso_SX,modelmass_SX,b_SX);
fgetMetabolicEnergySmooth2016all = ...
    Function('fgetMetabolicEnergySmooth2016all',...
    {exc_SX,act_SX,lMtilde_SX,vM_SX,Fce_SX,pctst_SX,Fiso_SX,modelmass_SX,b_SX},...
        {energy_total_SX,energy_am_SX,energy_sl_SX,energy_mech_SX,energy_model_SX});

%% store functions in structure
Fmetab.fgetMetabolicEnergySmooth2004all = fgetMetabolicEnergySmooth2004all;
Fmetab.fgetMetabolicEnergySmooth2003all = fgetMetabolicEnergySmooth2003all;
Fmetab.fgetMetabolicEnergySmooth2010all = fgetMetabolicEnergySmooth2010all;
Fmetab.fgetMetabolicEnergySmooth2016all = fgetMetabolicEnergySmooth2016all;

end