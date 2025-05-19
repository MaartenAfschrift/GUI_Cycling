function [Out] = PredSim_Cycling(S)
%PredSim_Cycling Predictive simulations of cycling
%   Detailed explanation goes here

import casadi.*
scaling = S.scaling;
SolverSetup = S.SolverSetup;



%% Load function to evaluate muscle and sekeleton dynamics

% To Do;
%   (1): add opensim model (using PredSim functions)
%
if strcmp(S.MSK_Model,'Kistemaker')
    % get muscle parameters
    MuscleParams = LoadParameters_hamsplit_DAK();
    % set specific tension
    MuscleParams.tension = S.specific_tension;
    % get casadi equations for rigid body dynamics
    [CasF_RBD,RBD_parms] = RigidBodyDyn_VU_CyclingModel(S);
    % equations metabolic energy
    CasF_metab = CasF_MetabolicEnergy(MuscleParams);
    if (strcmp(S.MuscleModel,'VU') || strcmp(S.MuscleModel,'VU-FV-simple'))
        % get casadi equations for muscle dynamics
        % CasF_Muscles = GetCasEq_MuscleDynamics_VU(MuscleParams);
        CasF_Muscles = GetCasEq_MuscleDynamics_All(MuscleParams);
    elseif strcmp(S.MuscleModel,'Leuven')
        % get casadi equations for muscle dynamics
        CasF_Muscles = GetCasEq_MuscleDynamics_Leuven(MuscleParams);
    elseif strcmp(S.MuscleModel,'hybrid')
        % get casadi equations for muscle dynamics -- hybrid
        CasF_Muscles = GetCasEq_MuscleDynamics_All(MuscleParams);
    end
end

%% check some inputs
[nr,nc] = size(S.Cycling.Saddle);
if nr == 1 && nc == 2
    S.Cycling.Saddle = S.Cycling.Saddle';
end


%% To Do':

% find feasible saddle position


% search for initial configuration that satisfies the kinematic constraints
% % create opti instance
% opti_Init = Opti();
% fi = opti_Init.variable(nseg,1);
% fid = opti_Init.variable(nseg,1);
% 
% % constraint on hip position and velocity
% [rhc,vhc] = eval_hip(fi,fid);
% opti_Init.subject_to(rhc == [-0.0549; 0.8743]);
% opti_Init.subject_to(vhc == [0; 0]);
% 
% % initial position cranck
% opti_Init.subject_to(fi(1) == 0.5*pi); % fi(1) at top-dead center
% opti_Init.subject_to(fid(1) == v_crank) % fid(1);
% 
% % path constraints
% opti_Init.subject_to(-pi <= fi(3)-fi(2) <= 0);
% opti_Init.subject_to(0 <= fi(4)-fi(3) <= pi);
% opti_Init.subject_to(-pi <= 1.1-fi(4,:) <= 0);
% 
% opti_Init.solver(SolverSetup.nlp.solver,SolverSetup.optionssol);
% 
% % some random objective function
% J = sumsqr(fi) + sumsqr(fid);
% opti_Init.minimize(J);
% % solve the optimal control problem
% solInit = opti_Init.solve();
% % extract solution
% IG.fi = solInit.value(fi);
% IG.fid = solInit.value(fid);




%% Optimal control Torque driven model
% we do a torque driven simulation and use this as an initial guess
opti = Opti();

% unpack settings
N = S.Coll.N;
cf =  S.Cycling.cf;
if S.Cycling.isokin
    v_crank = -2*pi*cf;
    time = 1/cf; % duration time [s]
    h = time/(N-1); % time interval between nodes
end
tVect = linspace(0,time,N);
nseg    = 4;
Topt    = 300;

% create optimization variables
Sfi = opti.variable(nseg,N);
Sfid = opti.variable(nseg,N);
Sfidd = opti.variable(nseg,N);
stim = opti.variable(3,N);
SFr5x = opti.variable(1,N);
SFr5y = opti.variable(1,N);
SM1 = opti.variable(1,N);

% unscale opt variables
fi = Sfi.*scaling.fi;
fid = Sfid.*scaling.fid;
fidd = Sfidd.*scaling.fidd;
Fr5x= SFr5x.*scaling.Fr5x;
Fr5y= SFr5y.*scaling.Fr5y;
M1= SM1.*scaling.M1;

% set some bounds on variables
opti.subject_to(pi/2+0.3<fi(2,:)<pi+0.3)
opti.subject_to(-pi <= fi(3,:)-fi(2,:) <= 0);
opti.subject_to(0 <= fi(4,:)-fi(3,:) <= pi);
opti.subject_to(-pi <= 1.1-fi(4,:) <= 0);
opti.subject_to(-1 < stim < 1);

% initial guess (based on something that satisfies )
fi_guess = [1.5708;    0.9117;    0.9117;    2.6898];
fid_guess = [v_crank; 1.2534; 3.4792; 2.5262];
Fr5x_guess = -200;
Fr5y_guess = -200;
M1_guess = -50;
opti.set_initial(Sfi,repmat(fi_guess,1,N)./scaling.fi);
opti.set_initial(Sfid,repmat(fid_guess,1,N)./scaling.fid);
opti.set_initial(Sfidd,0);
opti.set_initial(SFr5x,Fr5x_guess./scaling.Fr5x);
opti.set_initial(SFr5y,Fr5y_guess./scaling.Fr5y);
opti.set_initial(SM1,M1_guess./scaling.M1);

% constraints on initial state (angle crank)
opti.subject_to(fi(1,1) == 0.5*pi); % fi(1) at top-dead center

% constant crank angular velocity during motion
opti.subject_to(fid(1,:) == v_crank) % fid(1);d

% try first with an trapezoidal integration scheme
x = [fi(:,1:N-1); fid(:,1:N-1) ];
xt1 = [fi(:,2:N); fid(:,2:N) ];
xd = [fid(:,1:N-1); fidd(:,1:N-1) ];
xdt1 = [fid(:,2:N); fidd(:,2:N) ];
errInt = TrapezoidalIntegrator(x,xt1,xd,xdt1,h);
opti.subject_to(errInt ==  0);

% pre allocte some matrices for constraints
Terr =  MX(nseg,N);
rhc = MX(2,N);
% loop over mesh points
for k = 1:N
    % evaluate skeleton dynamics
    Tm = stim(:,k).*Topt;
    Terr(:,k)  = CasF_RBD.f_skeldyn(fi(:,k),fid(:,k), fidd(:,k),...
        Fr5x(1,k),Fr5y(1,k),M1(1,k),Tm);
    % Constraints on hip and hip velocity
    [rhcTemp,vhcTemp] = CasF_RBD.f_Hip_kin(fi(:,k),fid(:,k));
    rhc(:,k) = rhcTemp;
end

% impose equality constraints
opti.subject_to(Terr == 0)
opti.subject_to(rhc == repmat(S.Cycling.Saddle,1,N));

% % impose constraint on periodic motion
opti.subject_to(fi(2,1) == fi(2,end));
opti.subject_to(fid(2,1) == fid(2,end));

% equality constraint on a desired power
AveragePower = sum(-fid(1,:).*M1)/N;
opti.subject_to(AveragePower == 200)

% objective function (here minimize stim mainly)
J = sumsqr(stim) + 0.001.*sumsqr(fidd); % minimize power by moment on cranck
opti.solver(SolverSetup.nlp.solver,SolverSetup.optionssol);
opti.minimize(J);

% solve the optimal control problem
sol = opti.solve();

% unpack solution
Tdriven.fi = sol.value(fi);
Tdriven.fid = sol.value(fid);
Tdriven.fidd = sol.value(fidd);
Tdriven.T = [sol.value(M1); sol.value(stim)*Topt];
Tdriven.Fr5x = sol.value(Fr5x);
Tdriven.Fr5y = sol.value(Fr5y);
Tdriven.M1 = sol.value(M1);
Tdriven.AvPower = sol.value(AveragePower);
Tdriven.time = tVect;
Tdriven.JointPower = Tdriven.T.*Tdriven.fid;
Tdriven.JointWork = trapz(tVect',Tdriven.JointPower');
Tdriven.rhc = sol.value(rhc);

clear opti sol

%% compute initial guess fiber lengths based on rigid tendon assumption
% (or based on fiber lengths ?)


%% Formulate optimization problem (Muscle driven simulation)

% create opti object
opti = Opti();

% unpack settings
N = S.Coll.N;
if ~S.Cycling.Opt_cf
    cf =  S.Cycling.cf;
else
    cf = opti.variable(1,1);
    opti.set_initial(cf,S.Cycling.cf);
    opti.subject_to(30/60 < cf < 200/60);
end
if S.Cycling.isokin
    v_crank = -2*pi*cf;
    time = 1/cf; % duration time [s]
    h = time/(N-1); % time interval between nodes
end

% create optimization variables
nseg    = 4;
nmus    = MuscleParams.nmus;
Sfi     = opti.variable(nseg,N);
Sfid    = opti.variable(nseg,N);
a       = opti.variable(nmus,N);
lcerel  = opti.variable(nmus,N);

Sfidd   = opti.variable(nseg,N);
stim    = opti.variable(nmus,N);
SFr5x   = opti.variable(1,N);
SFr5y   = opti.variable(1,N);
SM1     = opti.variable(1,N);
vcerel_Helper = opti.variable(nmus,N);
if strcmp(S.MuscleModel,'VU')
    fcerel_Helper = opti.variable(nmus,N);
end

% unscale opt variables
fi = Sfi.*scaling.fi;
fid = Sfid.*scaling.fid;
fidd = Sfidd.*scaling.fidd;
Fr5x= SFr5x.*scaling.Fr5x;
Fr5y= SFr5y.*scaling.Fr5y;
M1= SM1.*scaling.M1;

% set some bounds on variables
opti.subject_to(0 < stim < 1);
opti.subject_to(0 < a < 1);
opti.subject_to(0.1 < lcerel < 1.7);
opti.subject_to(-10< vcerel_Helper < 10);
if strcmp(S.MuscleModel,'VU')
    opti.subject_to(0 < fcerel_Helper < 3);
end

% set initial guess
opti.set_initial(Sfi,Tdriven.fi./scaling.fi);
opti.set_initial(Sfid,Tdriven.fid./scaling.fid);
opti.set_initial(Sfidd,Tdriven.fidd./scaling.fidd);
opti.set_initial(SFr5x,Tdriven.Fr5x./scaling.Fr5x);
opti.set_initial(SFr5y,Tdriven.Fr5y./scaling.Fr5y);
opti.set_initial(SM1,Tdriven.M1./scaling.M1);
opti.set_initial(a,1);
opti.set_initial(lcerel,1);
opti.set_initial(stim,1);
opti.set_initial(vcerel_Helper,0);
if strcmp(S.MuscleModel,'VU')
    opti.set_initial(fcerel_Helper,0.2);
end

% constraints on initial state (angle crank)
opti.subject_to(fi(1,1) == 0.5*pi); % fi(1) at top-dead center

% constant crank angular velocity during motion
opti.subject_to(fid(1,:) == v_crank) % fid(1);

% Pre allocate some variables
adot    = MX(nmus,N);
T_id_err = MX(nseg,N);
rhc     = MX(2,N);
vhc     = MX(2,N);
errFEq  = MX(nmus,N); % equilibrium forces
fseM    = MX(nmus,N);
TMusM   = MX(3,N);
FceV    = MX(nmus,N);
FpasV   = MX(nmus,N);
FMltildeV = MX(nmus,N);
FMvtildeV = MX(nmus,N);
FMV     = MX(nmus,N);
energy_total = MX(nmus,N);
energy_am = MX(nmus,N);
energy_a = MX(nmus,N);
energy_m = MX(nmus,N);
energy_sl = MX(nmus,N);
energy_mech = MX(nmus,N);

if strcmp(S.MuscleModel,'VU')
    vcerel = MX(nmus,N);
    fcerel = MX(nmus,N);
end
% loop over mesh points
for k = 1:N
    if strcmp(S.MSK_Model,'Kistemaker')
        % muscle geometry
        [lmt, momarm] = CasF_Muscles.f_lMT_dM(fi(:,k));
        if strcmp(S.MuscleModel,'Leuven')
            % torque generated by muscles
            fse =CasF_Muscles.f_FL_Tendon_Leuven(lcerel(:,k),fi(:,k));
        elseif (strcmp(S.MuscleModel,'VU') || strcmp(S.MuscleModel,'VU-FV-simple'))
            % torque generated by muscles
            fse =CasF_Muscles.f_FL_Tendon_Leuven(lcerel(:,k),fi(:,k));
            % fse =  CasF_Muscles.f_FL_tendon(lcerel(:,k),fi(:,k));
        elseif strcmp(S.MuscleModel,'hybrid')
            % torque generated by muscles
            fse =CasF_Muscles.f_FL_Tendon_Leuven(lcerel(:,k),fi(:,k));
            % fse =  CasF_Muscles.f_FL_tendon(lcerel(:,k),fi(:,k));            
        end
        TMus = fse'*momarm;
        fseM(:,k) = fse;
        TMusM(:,k) = TMus;
    end

    % evaluate skeletal dynamics at current mesh point
    T_id_err(:,k)  = CasF_RBD.f_skeldyn(fi(:,k),fid(:,k), fidd(:,k),...
        Fr5x(1,k),Fr5y(1,k),M1(1,k),TMus');

    % Constraints on hip and hip velocity
    [rhcTemp,vhcTemp] = CasF_RBD.f_Hip_kin(fi(:,k),fid(:,k));
    rhc(:,k) = rhcTemp;
    vhc(:,k) = vhcTemp;

    % get state derivative act
    adot(:,k) = CasF_Muscles.f_adot(a(:,k),stim(:,k));

    % get force length - velocity properties
    if strcmp(S.MuscleModel,'Leuven')
        % force length velocity
        [Fpe,FMltilde,FMvtilde] = CasF_Muscles.f_FLV_Leuven(lcerel(:,k),...
            vcerel_Helper(:,k));
        % active muscle force
        Fce = a(:,k).*FMltilde.*FMvtilde;
        % total muscle force
        FM = Fce+Fpe;
        % error force equilibrium (force muscle = force tendon)
        errFEq(:,k) =  FM-(fse./MuscleParams.fmax');
    elseif strcmp(S.MuscleModel,'VU')
        % get state derivate fiber length
        vcerel(:,k) = CasF_Muscles.f_FV_muscle(lcerel(:,k),a(:,k),fcerel_Helper(:,k));
        % constraint on help force muscle
        fcerel(:,k) = CasF_Muscles.f_fcerel(lcerel(:,k),fi(:,k));
        Fpe = zeros(nmus,1);
        Fce = fcerel(:,k);
        FMltilde = ones(nmus,1); % this is wrong: ToDo adapt (only for energy equations)
        FMvtilde = nan(nmus,1); % only to store results (ToDo: compute this)
        FM = Fce+Fpe;
    elseif strcmp(S.MuscleModel,'VU-FV-simple')
        % get state derivate fiber length
        FMltilde = CasF_Muscles.f_FL_muscle(lcerel(:,k));
        FMvtilde = CasF_Muscles.f_FV_muscle_simple(vcerel_Helper(:,k));
        Fce = a(:,k).*FMltilde.*FMvtilde;
        FM = Fce;
        Fpe = zeros(size(FM));
        % error force equilibrium (force muscle = force tendon)
        errFEq(:,k) =  FM-(fse./MuscleParams.fmax');
    elseif strcmp(S.MuscleModel,'hybrid')
        % force length velocity
        % [Fpe,~,FMvtilde] = CasF_Muscles.f_FLV_Leuven(lcerel(:,k),...
        %     vcerel_Helper(:,k));
        FMltilde = CasF_Muscles.f_FL_muscle(lcerel(:,k));
        FMvtilde = CasF_Muscles.f_FV_muscle_simple(vcerel_Helper(:,k));
        % active muscle force
        Fce = a(:,k).*FMltilde.*FMvtilde;
        % total muscle force
        Fpe = zeros(size(Fce));
        FM = Fce; % removed passive forces
        % error force equilibrium (force muscle = force tendon)
        errFEq(:,k) =  FM-(fse./MuscleParams.fmax');
    end
    % metabolic power
    vM = MuscleParams.lce_opt'.*vcerel_Helper(:,k);
    modelmass = 75;
    b = 100;
    pctst = 0.5;
    Fm = Fce'.*MuscleParams.fmax;
    if strcmp(S.Metab.Model,'Bhargava2004')
        [energy_total(:,k), energy_a(:,k), energy_m(:,k), energy_sl(:,k),...
            energy_mech(:,k)] = CasF_metab.fgetMetabolicEnergySmooth2004all(stim(:,k)',...
            a(:,k)',lcerel(:,k)',vM',Fm',Fpe',pctst,FMltilde',modelmass,b,...
            S.Metab.scaleRate);
    elseif strcmp(S.Metab.Model,'Umberger2003')
        [energy_total(:,k), energy_am(:,k), energy_sl(:,k),...
            energy_mech(:,k)]  = CasF_metab.fgetMetabolicEnergySmooth2003all(stim(:,k)',...
            a(:,k)',lcerel(:,k)',vM',Fm',pctst,FMltilde',modelmass,b);
    elseif strcmp(S.Metab.Model,'Umberger2010')
        [energy_total(:,k), energy_am(:,k), energy_sl(:,k),...
            energy_mech(:,k)] = CasF_metab.fgetMetabolicEnergySmooth2010all(stim(:,k)',...
            a(:,k)',lcerel(:,k)',vM',Fm',pctst,FMltilde',modelmass,b);
    elseif strcmp(S.Metab.Model,'Umberger2016')
        [energy_total(:,k), energy_am(:,k), energy_sl(:,k),...
            energy_mech(:,k)] = CasF_metab.fgetMetabolicEnergySmooth2016all(stim(:,k)',...
            a(:,k)',lcerel(:,k)',vM',Fm',pctst,FMltilde',modelmass,b);
    end

    % store some variables
    FceV(:,k)   = Fce;
    FpasV(:,k)  = Fpe;
    FMvtildeV(:,k) = FMvtilde;
    FMltildeV(:,k) = FMltilde;
    FMV(:,k) = FM;
end

% equality constraints
if strcmp(S.MuscleModel,'Leuven') || strcmp(S.MuscleModel,'VU-FV-simple')
    opti.subject_to(errFEq == 0); % equilibrium tendon and muscle force
elseif strcmp(S.MuscleModel,'VU')
    % constraint on helper variables
    opti.subject_to(vcerel_Helper == vcerel);
    opti.subject_to(fcerel_Helper == fcerel);
elseif strcmp(S.MuscleModel,'hybrid')
    opti.subject_to(errFEq == 0);
end
opti.subject_to(rhc == repmat(S.Cycling.Saddle,1,N)); % fixed position pelvis
opti.subject_to(T_id_err == 0); % inverse dynamics eq

% trapezoidal integration
x = [   fi(:,1:N-1);    fid(:,1:N-1);   a(:,1:N-1);     lcerel(:,1:N-1)];
xt1 = [ fi(:,2:N);      fid(:,2:N);     a(:,2:N);       lcerel(:,2:N)];
xd = [  fid(:,1:N-1);   fidd(:,1:N-1);  adot(:,1:N-1);  vcerel_Helper(:,1:N-1)];
xdt1 = [fid(:,2:N);     fidd(:,2:N);    adot(:,2:N);    vcerel_Helper(:,2:N)];
errInt = TrapezoidalIntegrator(x,xt1,xd,xdt1,h);
opti.subject_to(errInt ==  0);

% power cranck
CranckPower = -fid(1,:).*M1;
AveragePower = sum(CranckPower)./N;
if ~isnan(S.Cycling.FixedPower) && ~strcmp(S.Objective.type,'maxpower')
    opti.subject_to(AveragePower == S.Cycling.FixedPower);
end

% impose constraint on periodic motion for states (skeleton and muscle).
% note that we have only 1dof in the kinematic model, so imposing this at
% 1 joint is sufficient
opti.subject_to(fi(2,1) == fi(2,end));
opti.subject_to(fid(2,1) == fid(2,end));
opti.subject_to(lcerel(:,1) == lcerel(:,end));
opti.subject_to(vcerel_Helper(:,1) == vcerel_Helper(:,end));

% solver
opti.solver(SolverSetup.nlp.solver,SolverSetup.optionssol);

% objective function
if strcmp(S.Objective.type,'Multi_a_E')
    Jmetab = S.Objective.w_metab*sum(sum(energy_total))./N./nmus;
    Jstim = S.Objective.w_stim*sumsqr(stim)./N./nmus;
    Jfidd = S.Objective.w_qdd.*sumsqr(fidd)./N./nseg;
    Jvce = S.Objective.w_vMtilde.*sumsqr(vcerel_Helper)./N./nmus;
    JM1 = S.Objective.w_M1.*sumsqr(M1)./N./nmus;
    J_ConF = S.Objective.C_forces.*sum(Fr5x.^2 + Fr5y.^2)./N;
    J = S.Objective.scale.*(Jmetab + Jstim + Jfidd + Jvce + JM1 + J_ConF);
elseif strcmp(S.Objective.type,'stim')
    Jstim = S.Objective.w_stim*sumsqr(stim)./N./nmus;
    Jvce = S.Objective.w_vMtilde.*sumsqr(vcerel_Helper);
    JM1 = S.Objective.w_M1.*sumsqr(M1)./N./nmus;
    J_ConF = S.Objective.C_forces.*sum(Fr5x.^2 + Fr5y.^2)./N;
    J = S.Objective.scale.*(Jstim + Jvce + JM1 + J_ConF);
elseif strcmp(S.Objective.type,'stim_echt')
    Jstim = S.Objective.w_stim*sumsqr(stim)./N./nmus;
    J_ConF = S.Objective.C_forces.*sum(Fr5x.^2 + Fr5y.^2)./N;
    J = S.Objective.scale.*(Jstim + J_ConF) ;
elseif strcmp(S.Objective.type,'maxpower')
    % we want to maximize power and we don't want negative values (must be possible to improve this)
    JCranckPower = S.Objective.w_cranckP.*(AveragePower.*-1 + 3000);
    Jstim = S.Objective.w_stim*sumsqr(stim)./N./nmus;
    Jvce = S.Objective.w_vMtilde.*sumsqr(vcerel_Helper);
    JM1 = S.Objective.w_M1.*sumsqr(M1)./N./nmus;
    J_ConF = S.Objective.C_forces.*sum(Fr5x.^2 + Fr5y.^2)./N;
    J = S.Objective.scale.*(Jstim + Jvce + JM1 + JCranckPower + J_ConF);
elseif strcmp(S.Objective.type,'MinNegFiberWork')
    Jstim = S.Objective.w_stim*sumsqr(stim)./N./nmus;
    Jvce = S.Objective.w_vMtilde.*sumsqr(vcerel_Helper);
    JM1 = S.Objective.w_M1.*sumsqr(M1)./N./nmus;
    J_ConF = S.Objective.C_forces.*sum(Fr5x.^2 + Fr5y.^2)./N;
    Pmuscle = (FceV.*MuscleParams.fmax') .* (MuscleParams.lce_opt'.*vcerel_Helper);
    PmuscleNeg = Pmuscle.*(0.5.*tanh(0.1.*-Pmuscle)+0.5);
    J_NegFiberWork = S.Objective.w_minNegWork.*(-sum(sum(PmuscleNeg))./N);
    % J = S.Objective.scale.*(Jstim + Jvce + JM1 + J_ConF + J_NegFiberWork);
    J = S.Objective.scale.*(J_ConF + J_NegFiberWork);

else
    disp('Not yet supported, minimizing stim');
    Jstim = S.Objective.w_stim*sumsqr(stim)./N./nmus;
    J = S.Objective.scale.*(Jstim);
end
opti.minimize(J);

Pmuscle = (FceV.*MuscleParams.fmax') .* (MuscleParams.lce_opt'.*vcerel_Helper);

% solve the optimal control problem
sol = opti.solve();

% unpack solution
Out.fi = sol.value(fi);
Out.fid = sol.value(fid);
Out.fidd = sol.value(fidd);
Out.Fr5x = sol.value(Fr5x);
Out.Fr5y = sol.value(Fr5y);
Out.M1 = sol.value(M1);
Out.stim = sol.value(stim);
Out.TMus = sol.value(TMusM);
Out.fse = sol.value(fseM);
Out.a = sol.value(a);
Out.lMtilde = sol.value(lcerel);
Out.Fce = sol.value(FceV);
Out.vMtilde = sol.value(vcerel_Helper);
Out.FMltilde = sol.value(FMltildeV);
Out.FMvtilde = sol.value(FMvtildeV);
Out.Fpas = sol.value(FpasV);
Out.FM = sol.value(FMV);
% energy equations
Out.energy_total = sol.value(energy_total);
Out.energy_a = sol.value(energy_a);
Out.energy_m = sol.value(energy_m);
Out.energy_am = sol.value(energy_am);
Out.energy_mech = sol.value(energy_mech);
Out.Pmuscle = sol.value(Pmuscle);

Out.J = sol.value(J);

if ~S.Cycling.Opt_cf
    Out.t = tVect;
    Out.cf = cf;
else
    Out.cf = sol.value(cf);
    Out.t = linspace(0,1/Out.cf,N);
end

if strcmp(S.Objective.type,'Multi_a_E')
    Out.Jfidd = S.Objective.scale.*sol.value(Jfidd);
    Out.Jmetab = S.Objective.scale.*sol.value(Jmetab);
    Out.J_ConF = S.Objective.scale.*sol.value(J_ConF);
    Out.Jstim = S.Objective.scale.*sol.value(Jstim);
    Out.Jvce = S.Objective.scale.*sol.value(Jvce);
    Out.JM1 = S.Objective.scale.*sol.value(JM1);
    Out.JCranckPower = NaN;
    Out.J_NegFiberWork = NaN;
end
if strcmp(S.Objective.type,'maxpower')
    Out.JCranckPower = S.Objective.scale.*sol.value(JCranckPower);
    Out.Jfidd = NaN;
    Out.Jmetab = NaN;
    Out.J_ConF = S.Objective.scale.*sol.value(J_ConF);
    Out.Jstim = S.Objective.scale.*sol.value(Jstim);
    Out.Jvce = S.Objective.scale.*sol.value(Jvce);
    Out.JM1 = S.Objective.scale.*sol.value(JM1);
    Out.J_NegFiberWork = NaN;
end
if strcmp(S.Objective.type,'stim')
    Out.JCranckPower = NaN;
    Out.Jfidd = NaN;
    Out.Jmetab = NaN;
    Out.J_ConF = S.Objective.scale.*sol.value(J_ConF);
    Out.Jstim = S.Objective.scale.*sol.value(Jstim);
    Out.Jvce = S.Objective.scale.*sol.value(Jvce);
    Out.JM1 = S.Objective.scale.*sol.value(JM1);
    Out.J_NegFiberWork = NaN;
end

if strcmp(S.Objective.type,'MinNegFiberWork')
    Out.JCranckPower = NaN;
    Out.Jfidd = NaN;
    Out.Jmetab = NaN;
    Out.J_ConF = S.Objective.scale.*sol.value(J_ConF);
    Out.Jstim = S.Objective.scale.*sol.value(Jstim);
    Out.Jvce = S.Objective.scale.*sol.value(Jvce);
    Out.JM1 = S.Objective.scale.*sol.value(JM1);
    Out.J_NegFiberWork = S.Objective.scale.*sol.value(J_NegFiberWork);
end

Out.CranckPower  = sol.value(CranckPower);

% display some results
disp(['average power crank ' num2str(mean(Out.CranckPower)) ' W']);
Acranck = trapz(Out.t',Out.CranckPower');
Ametab = sum(trapz(Out.t',Out.energy_total'));
Pmuscle = -(Out.Fce.*MuscleParams.fmax').* (Out.vMtilde.*MuscleParams.lce_opt');
Amech = sum(trapz(Out.t',Pmuscle'));
disp(['metabolic work ' num2str(Ametab) ' J']);
disp(['Mechanical work cranck ' num2str(Acranck) ' J']);
disp(['Muscle mechanical work ' num2str(Amech) ' J']);
disp(['average metabolic effiency ' num2str(Amech./Ametab) ' %']);
Out.Acranck = Acranck;
Out.Ametab = Ametab;
Out.Pmuscle = Pmuscle;
Out.Amech= Amech;
Out.stats = sol.stats;
Out.S = S;
Out.MuscleParams = MuscleParams;
Out.RBD_parms = RBD_parms;

% joint locations
Out.rJoint = getJointPositions(RBD_parms,Out.fi);


%% compute metabolic energy for all models during postprocessing

for k = 1:N
    if strcmp(S.MuscleModel,'Leuven') && strcmp(S.MSK_Model,'Kistemaker')
        % torque generated by muscles
        fse =CasF_Muscles.f_FL_Tendon_Leuven(Out.lMtilde(:,k),Out.fi(:,k));
        [lmt, momarm] = CasF_Muscles.f_lMT_dM(Out.fi(:,k));
        TMus = fse'*momarm;
    elseif (strcmp(S.MuscleModel,'VU') || strcmp(S.MuscleModel,'VU-FV-simple')) ...
            && strcmp(S.MSK_Model,'Kistemaker')
        % torque generated by muscles
        fse =  CasF_Muscles.f_FL_tendon(Out.lMtilde(:,k),Out.fi(:,k));
        [lmt, momarm] = CasF_Muscles.f_lMT_dM(fi(:,k));
        TMus = fse'*momarm;
    end

    % get force length - velocity properties
    if strcmp(S.MuscleModel,'Leuven')
        % force length velocity
        [Fpe,FMltilde,FMvtilde] = CasF_Muscles.f_FLV_Leuven(Out.lMtilde(:,k),...
            Out.vMtilde(:,k));
        % active muscle force
        Fce = Out.a(:,k).*FMltilde.*FMvtilde;
        % total muscle force
        FM = Fce+Fpe;
        % elseif strcmp(S.MuscleModel,'VU')
        %     % get state derivate fiber length
        %     vcerel(:,k) = CasF_Muscles.f_FV_muscle(Out.lMtilde(:,k),a(:,k),fcerel_Helper(:,k));
        %     % constraint on help force muscle
        %     fcerel(:,k) = CasF_Muscles.f_fcerel(lcerel(:,k),fi(:,k));
        %     Fpe = zeros(nmus,1);
        %     Fce = fcerel(:,k);
        %     FMltilde = ones(nmus,1); % this is wrong: ToDo adapt (only for energy equations)
        %     FMvtilde = nan(nmus,1); % only to store results (ToDo: compute this)
        %     FM = Fce+Fpe;
    elseif strcmp(S.MuscleModel,'VU-FV-simple')
        % get state derivate fiber length
        FMltilde = CasF_Muscles.f_FL_muscle(Out.lMtilde(:,k));
        FMvtilde = CasF_Muscles.f_FL_muscle(Out.vMtilde(:,k));
        Fce = Out.a(:,k).*FMltilde.*FMvtilde;
        FM = Fce;
    end

    % metabolic power
    vM = MuscleParams.lce_opt'.*Out.vMtilde(:,k);
    modelmass = 75;
    b = 100000;
    pctst = 0.5;
    Fm = Fce'.*MuscleParams.fmax;
    % Bhargava 2004
    [a, bb, c, d, e] = CasF_metab.fgetMetabolicEnergySmooth2004all(Out.stim(:,k)',...
        Out.a(:,k)',Out.lMtilde(:,k)',vM',Fm',Fpe',pctst,FMltilde',modelmass,b,...
        S.Metab.scaleRate);
    Bharg.energy_total(:,k)= full(a);
    Bharg.energy_a(:,k) = full(bb);
    Bharg.energy_m(:,k) = full(c);
    Bharg.energy_sl(:,k)= full(d);
    Bharg.energy_mech(:,k) = full(e);
    % Umberger 2003
    [a, bb, c, d]  = CasF_metab.fgetMetabolicEnergySmooth2003all(Out.stim(:,k)',...
        Out.a(:,k)',Out.lMtilde(:,k)',vM',Fm',pctst,FMltilde',modelmass,b);
    Umb2003.energy_total(:,k) = full(a);
    Umb2003.energy_am(:,k)= full(bb);
    Umb2003.energy_sl(:,k) = full(c);
    Umb2003.energy_mech(:,k) = full(d);
    % Umberger 2010
    [a, bb, c, d] = CasF_metab.fgetMetabolicEnergySmooth2010all(Out.stim(:,k)',...
        Out.a(:,k)',Out.lMtilde(:,k)',vM',Fm',pctst,FMltilde',modelmass,b);
    Umb2010.energy_total(:,k) = full(a);
    Umb2010.energy_am(:,k)= full(bb);
    Umb2010.energy_sl(:,k) = full(c);
    Umb2010.energy_mech(:,k) = full(d);
    % Umberger 2016
    [a, bb, c, d] = CasF_metab.fgetMetabolicEnergySmooth2016all(Out.stim(:,k)',...
        Out.a(:,k)',Out.lMtilde(:,k)',vM',Fm',pctst,FMltilde',modelmass,b);
    Umb2016.energy_total(:,k) = full(a);
    Umb2016.energy_am(:,k)= full(bb);
    Umb2016.energy_sl(:,k) = full(c);
    Umb2016.energy_mech(:,k) = full(d);
end

Out.Metab.Bharg = Bharg;
Out.Metab.Umb2003 = Umb2003;
Out.Metab.Umb2010 = Umb2010;
Out.Metab.Umb2016 = Umb2016;

%% plot solution
if S.BoolPlot
    SetFigureDefaults();
    figure('Color',[1 1 1]);
    nr = 4;
    for i=1:4
        subplot(nr,4,i)
        plot(Out.fi(i,:));
        ylabel(['seg orientation ' num2str(i)]);
        subplot(nr,4,i+4)
        if i<4
            plot(Out.TMus(i,:));
        end
        ylabel(['joint moment ' num2str(i)]);
        subplot(nr,4,8)
        plot(Out.M1);
        ylabel('M1 moment');
    end
    subplot(nr,4,9)
    plot(Out.CranckPower);
    ylabel('Power Cranck')
    subplot(nr,4,10)
    bar(1,Out.Jmetab./Out.J); hold on;
    bar(2,Out.Jfidd ./Out.J);
    bar(3,Out.Jstim ./Out.J);
    bar(4,Out.JM1 ./Out.J);
    bar(5,Out.Jvce ./Out.J);
    bar(6,Out.J_ConF ./Out.J);
    bar(7,Out.JCranckPower ./Out.J);
    bar(8,Out.J_NegFiberWork./Out.J)
    set(gca,'XTickLabelRotation',45);
    set(gca,'XTick',1:6);
    set(gca,'XTickLabel',{'E','fidd','a','Ma','vM','Fc','Pmax','NegPm'});
    ylabel('% objective')
    subplot(nr,4,11:12);
    plot(Pmuscle');
    ylabel('power muscle');
    set(gca,'box','off');
    for i= 1:10
        subplot(nr,4,i)
        set(gca,'box','off');
    end
end


%% plot figure for F/l/v

if strcmp(S.MuscleModel,'VU') && strcmp(S.MSK_Model,'Kistemaker')
    % kracht lengte relatie
    lcerel = 0.3:0.01:1.5;
    Fl = [];
    for i=1:length(lcerel)
        x = repmat(lcerel(i),9,1);
        Fl_temp = CasF_Muscles.f_FL_muscle(x);
        Fl(i) = full(Fl_temp(1));
    end

    % kracht/snelheidsrelatie
    fcerel = 0:0.01:1.3;
    aV = [0 0.01 0.2 0.5 1];
    aV = [0.01];
    lcerel_V = 0.8;
    Fv = [];
    for i=1:length(fcerel)
        for j = 1:length(aV)
            x = repmat(fcerel(i),9,1);
            a = repmat(aV(j),9,1);
            Fv_temp =  CasF_Muscles.f_FV_muscle(lcerel_V,a,x);
            Fv(i,j) = full(Fv_temp(1));
        end
    end
    figure();
    subplot(1,2,1);
    plot(lcerel,Fl);
    xlabel('lrel');
    ylabel('FV');
    subplot(1,2,2);
    for j = 1:length(aV)
        plot(Fv(:,j),fcerel); hold on;
    end
    % legend('a0','a001','a02','á05','a1')
    set(gca,'ylim',[-0.2 2]);
    set(gca,'xlim',[-10 10]);
    xlabel('vrel');
    ylabel('FV');
end

end