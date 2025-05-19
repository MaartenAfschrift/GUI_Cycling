function [vcerel] = Muscle_F_V_VU(fcerel,q,lcerel,parms,varargin)
%Muscle_F_V_VU Force velocity muscle contraction VU
%   input arguments
%       (1): normalized muscle force
%       (2): active state
%       (3): normalized fiber length
%       (4): parameters
%       (5): varargin:
%           (5.1): bool_ActivationDependent
%           (5.2): bool_LengthDependent

                


bool_ActivationDependent = true;
if ~isempty(varargin)
    bool_ActivationDependent = varargin{1};
end

bool_LengthDependent = true;
if length(varargin)>1
    bool_LengthDependent = varargin{2};
end

% function for smooth transition
sigma = @(x,x0,w) 1./(1+exp(-(x-x0)/w));

% get maximal force for given length
if bool_LengthDependent
    fisomrel = Muscle_F_L_VU(lcerel);
else
    fisomrel = ones(size(fcerel));
    lcerel = ones(size(fcerel));
end

% dependence on active state
q0_b = (log(1/parms.vfactmin-1)+parms.q0*22)/22;
if bool_ActivationDependent
    brel = parms.brel_c./(1+exp(-22*(q-q0_b)));
else    
    brel = parms.brel_c;
end

sig_arel = sigma(lcerel,1,0.01);
arel = (1-sig_arel).*parms.arel_c + sig_arel.*(parms.arel_c.*fisomrel);
arel = arel.*fisomrel;

dvdf_isom_con=brel./(q.*(fisomrel+arel)); % slope in the isometric point at wrt concentric part
dvdf_isom_ecc=dvdf_isom_con./parms.slopfac; % slope in the isometric point at wrt eccentric part
dFdvcon0=1./dvdf_isom_con;
s_as = 1./parms.sloplin;
p1 = -(fisomrel.*q.*(parms.fasymp - 1))./(s_as - dFdvcon0.*parms.slopfac); 
p3 =  -parms.fasymp.*fisomrel.*q;
p2 =  (fisomrel.^2.*q.^2.*(parms.fasymp - 1).^2)./(s_as - dFdvcon0.*parms.slopfac); 
p4 =  -s_as;

% Concentric part
sig_c1 = sigma(fcerel,fisomrel.*q,0.01);
lcereld_c = (1-sig_c1).*brel.*(fcerel-q.*fisomrel)./(fcerel+q.*arel);
sig_c2 = sigma(dvdf_isom_con,parms.sloplin,0.01);
lcereld_c = lcereld_c + sig_c2.*parms.sloplin.*(fcerel-q.*fisomrel);

% Eccentric part
sig_e1 = sigma(fcerel,fisomrel.*q,0.01);

lcereld_e = sig_e1.*(-(fcerel + p3 + p1.*p4 + (fcerel.^2 - 2*fcerel.*p1.*p4 +...
    2*fcerel.*p3 + p1.^2.*p4.^2 - 2*p1.*p3.*p4 + p3.^2 + ...
    4*p2.*p4).^(1/2))./(2.*p4));
sig_e2 = sigma(dvdf_isom_ecc,parms.sloplin./parms.slopfac,0.01);
lcereld_e = lcereld_e + sig_e2.*(parms.sloplin./parms.slopfac).*(fcerel-q.*fisomrel);

% Final relationship
vcerel = lcereld_c + lcereld_e;




end