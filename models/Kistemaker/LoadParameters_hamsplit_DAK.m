function [parms] = LoadParameters_hamsplit_DAK()

%% load matlab structure with all the default parameters

mus = 1:9;
tmp = load('kvs_parms.mat');

%% Muscle parameters
% General
arel_c = 0.41;
brel_c = 5.2;
fasymp = 1.5;
slopfac = 2;
vfactmin = 0.1;
width = 0.56;
C_fl = -1/width^2;
b_Fse = 1;
sloplin = vfactmin.*brel_c./(slopfac.*0.005.*0.0975.*(1+arel_c));

% Muscle specific
tmp.A0 = [tmp.A0(:,3:5); 0 0.3388 0]; %  0.338839879800784
tmp.A1 = [tmp.A1(:,3:5); 0 -0.0260 0];
tmp.A2 = [tmp.A2(:,3:5); 0 0 0];

% Nieuwe fmax voor bi en mono ham
tmp.fmax = [tmp.ce_F_max tmp.ce_F_max(7)*(5/40)];
tmp.fmax(7) = tmp.fmax(7)*(35/40);
tmp.muscle_names = [tmp.muscle_names; '9 biceps femoris 2 '];

tmp.lce_opt = [tmp.ce_len_opt 0.11]; % Lce_opt_mono_ham = 0.11
tmp.lse0 = [tmp.se_len_slack 0.20]; % LSE lengte rust mono_ham

A0 = tmp.A0(mus,:);
A1 = tmp.A1(mus,:);
A2 = tmp.A2(mus,:);

lce_opt = tmp.lce_opt;

lse0 = tmp.lse0;

fmax = tmp.fmax;

kse = fmax./((lse0.*0.04).^2);

%% Activation dynamics parameters
q0 = tmp.hafo_q0(1);
rm = tmp.hafo_m(1);
gamma_0 = 0.00001;
kCa = 0.8e-5;
a_act = -4.587;
a1_act=log10(exp(a_act));
b_act = [5.168 1.081 -0.1909];

%% add parameters to output structure
parms.arel_c = arel_c;
parms.brel_c = brel_c; % coeff for force velocity
parms.fasymp = fasymp;
parms.slopfac = slopfac;
parms.vfactmin = vfactmin; % coeff for force velocity
parms.width = width;
parms.C_fl = C_fl;
parms.b_Fse = b_Fse; % coeff for F/L tendon
parms.sloplin = sloplin;

% muscle geometry
parms.A0 = A0;
parms.A1 = A1;
parms.A2 = A2;

% muscle properties
parms.fmax = fmax; % max force
parms.lce_opt = lce_opt; % optimal fiber length
parms.lse0 = lse0; % tendon slack length
parms.kse = kse; % coeff for F/L tendon

parms.q0 = q0;
parms.rm = rm;
parms.gamma_0 = gamma_0;
parms.kCa = kCa;
parms.a_act = a_act;
parms.a1_act= a1_act;
parms.b_act = b_act;
parms.nmus = length(mus);
parms.muscle_names = tmp.muscle_names;



end