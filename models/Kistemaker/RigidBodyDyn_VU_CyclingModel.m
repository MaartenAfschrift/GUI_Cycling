function [CasF,parms] = RigidBodyDyn_VU_CyclingModel(S)
%RigidBodyDyn_VU_CyclingModel Gets rigid body dynamics for VU cycling model

if ~isfield(S,'Cycling')
    S.Cycling.CrankLength = 0.17;
else
    if ~isfield(S.Cycling,'CrankLength')
        S.Cycling.CrankLength = 0.17;
    end
end
        


%% Skeletal dynamics parameters
l = [S.Cycling.CrankLength 0.165 0.458 0.4851];
d = [S.Cycling.CrankLength/2 0.120 0.260 0.275];
m = [0.200 1.234 3.540 8.470];
j = [0.001 0.010 0.068 0.209];
p = l-d;

parms.segparms.L = l;
parms.segparms.d = d;
parms.segparms.m = m;
parms.segparms.j = j;
parms.segparms.p = p;
parms.segparms.g = -9.81;

%% equations of motion

% Dinant derived these equations of motion (using the simbolic toolbox). It
% is essentially a 2D 4 link bar system with a 
nseg = length(l);
import casadi.*
fi     = SX.sym('fi', nseg,1);
fid    = SX.sym('fid', nseg,1);
fidd   = SX.sym('fidd',nseg,1);
Fr5x   = SX.sym('Fr5x',1,1);
Fr5y   = SX.sym('Fr5y',1,1);
M      = SX.sym('M',1,1);
Tmus   = SX.sym('Tmus',nseg-1,1);

tor = Tmus;

As = ...
[-j(1)-d(1)^(2)*m(1)*cos(fi(1))^(2)-d(1)^(2)*m(1)*sin(fi(1))^(2)-l(1)*m(2)*cos(fi(1))*(d(1)*cos(fi(1))-cos(fi(1))*(d(1)-l(1)))-l(1)*m(3)*cos(fi(1))*(d(1)*cos(fi(1))-cos(fi(1))*(d(1)-l(1)))-l(1)*m(4)*cos(fi(1))*(d(1)*cos(fi(1))-cos(fi(1))*(d(1)-l(1)))-l(1)*m(2)*sin(fi(1))*(d(1)*sin(fi(1))-sin(fi(1))*(d(1)-l(1)))-l(1)*m(3)*sin(fi(1))*(d(1)*sin(fi(1))-sin(fi(1))*(d(1)-l(1)))-l(1)*m(4)*sin(fi(1))*(d(1)*sin(fi(1))-sin(fi(1))*(d(1)-l(1))),-d(2)*m(2)*cos(fi(2))*(d(1)*cos(fi(1))-cos(fi(1))*(d(1)-l(1)))-l(2)*m(3)*cos(fi(2))*(d(1)*cos(fi(1))-cos(fi(1))*(d(1)-l(1)))-l(2)*m(4)*cos(fi(2))*(d(1)*cos(fi(1))-cos(fi(1))*(d(1)-l(1)))-d(2)*m(2)*sin(fi(2))*(d(1)*sin(fi(1))-sin(fi(1))*(d(1)-l(1)))-l(2)*m(3)*sin(fi(2))*(d(1)*sin(fi(1))-sin(fi(1))*(d(1)-l(1)))-l(2)*m(4)*sin(fi(2))*(d(1)*sin(fi(1))-sin(fi(1))*(d(1)-l(1))),-d(3)*m(3)*cos(fi(3))*(d(1)*cos(fi(1))-cos(fi(1))*(d(1)-l(1)))-l(3)*m(4)*cos(fi(3))*(d(1)*cos(fi(1))-cos(fi(1))*(d(1)-l(1)))-d(3)*m(3)*sin(fi(3))*(d(1)*sin(fi(1))-sin(fi(1))*(d(1)-l(1)))-l(3)*m(4)*sin(fi(3))*(d(1)*sin(fi(1))-sin(fi(1))*(d(1)-l(1))),-d(4)*m(4)*cos(fi(4))*(d(1)*cos(fi(1))-cos(fi(1))*(d(1)-l(1)))-d(4)*m(4)*sin(fi(4))*(d(1)*sin(fi(1))-sin(fi(1))*(d(1)-l(1)));
                                                                                             -l(1)*m(3)*cos(fi(1))*(d(2)*cos(fi(2))-cos(fi(2))*(d(2)-l(2)))-l(1)*m(4)*cos(fi(1))*(d(2)*cos(fi(2))-cos(fi(2))*(d(2)-l(2)))-l(1)*m(3)*sin(fi(1))*(d(2)*sin(fi(2))-sin(fi(2))*(d(2)-l(2)))-l(1)*m(4)*sin(fi(1))*(d(2)*sin(fi(2))-sin(fi(2))*(d(2)-l(2)))-d(2)*l(1)*m(2)*cos(fi(1))*cos(fi(2))-d(2)*l(1)*m(2)*sin(fi(1))*sin(fi(2)),                                                         -j(2)-d(2)^(2)*m(2)*cos(fi(2))^(2)-d(2)^(2)*m(2)*sin(fi(2))^(2)-l(2)*m(3)*cos(fi(2))*(d(2)*cos(fi(2))-cos(fi(2))*(d(2)-l(2)))-l(2)*m(4)*cos(fi(2))*(d(2)*cos(fi(2))-cos(fi(2))*(d(2)-l(2)))-l(2)*m(3)*sin(fi(2))*(d(2)*sin(fi(2))-sin(fi(2))*(d(2)-l(2)))-l(2)*m(4)*sin(fi(2))*(d(2)*sin(fi(2))-sin(fi(2))*(d(2)-l(2))),-d(3)*m(3)*cos(fi(3))*(d(2)*cos(fi(2))-cos(fi(2))*(d(2)-l(2)))-l(3)*m(4)*cos(fi(3))*(d(2)*cos(fi(2))-cos(fi(2))*(d(2)-l(2)))-d(3)*m(3)*sin(fi(3))*(d(2)*sin(fi(2))-sin(fi(2))*(d(2)-l(2)))-l(3)*m(4)*sin(fi(3))*(d(2)*sin(fi(2))-sin(fi(2))*(d(2)-l(2))),-d(4)*m(4)*cos(fi(4))*(d(2)*cos(fi(2))-cos(fi(2))*(d(2)-l(2)))-d(4)*m(4)*sin(fi(4))*(d(2)*sin(fi(2))-sin(fi(2))*(d(2)-l(2)));
                                                                                                                                                                                                     -l(1)*m(4)*cos(fi(1))*(d(3)*cos(fi(3))-cos(fi(3))*(d(3)-l(3)))-l(1)*m(4)*sin(fi(1))*(d(3)*sin(fi(3))-sin(fi(3))*(d(3)-l(3)))-d(3)*l(1)*m(3)*cos(fi(1))*cos(fi(3))-d(3)*l(1)*m(3)*sin(fi(1))*sin(fi(3)),                                                                                                                                                      -l(2)*m(4)*cos(fi(2))*(d(3)*cos(fi(3))-cos(fi(3))*(d(3)-l(3)))-l(2)*m(4)*sin(fi(2))*(d(3)*sin(fi(3))-sin(fi(3))*(d(3)-l(3)))-d(3)*l(2)*m(3)*cos(fi(2))*cos(fi(3))-d(3)*l(2)*m(3)*sin(fi(2))*sin(fi(3)),                                                         -j(3)-d(3)^(2)*m(3)*cos(fi(3))^(2)-d(3)^(2)*m(3)*sin(fi(3))^(2)-l(3)*m(4)*cos(fi(3))*(d(3)*cos(fi(3))-cos(fi(3))*(d(3)-l(3)))-l(3)*m(4)*sin(fi(3))*(d(3)*sin(fi(3))-sin(fi(3))*(d(3)-l(3))),-d(4)*m(4)*cos(fi(4))*(d(3)*cos(fi(3))-cos(fi(3))*(d(3)-l(3)))-d(4)*m(4)*sin(fi(4))*(d(3)*sin(fi(3))-sin(fi(3))*(d(3)-l(3)));
                                                                                                                                                                                                                                                                                                             -d(4)*l(1)*m(4)*cos(fi(1))*cos(fi(4))-d(4)*l(1)*m(4)*sin(fi(1))*sin(fi(4)),                                                                                                                                                                                                                                                              -d(4)*l(2)*m(4)*cos(fi(2))*cos(fi(4))-d(4)*l(2)*m(4)*sin(fi(2))*sin(fi(4)),                                                                                                                                                      -d(4)*l(3)*m(4)*cos(fi(3))*cos(fi(4))-d(4)*l(3)*m(4)*sin(fi(3))*sin(fi(4)),                                                         -m(4)*d(4)^(2)*cos(fi(4))^(2)-m(4)*d(4)^(2)*sin(fi(4))^(2)-j(4)];
 
 
bs = ...
[(d(1)*sin(fi(1))-sin(fi(1))*(d(1)-l(1)))*(l(1)*m(2)*cos(fi(1))*fid(1)^(2)+d(2)*m(2)*cos(fi(2))*fid(2)^(2))-M(1)-(d(1)*cos(fi(1))-cos(fi(1))*(d(1)-l(1)))*(l(1)*m(4)*sin(fi(1))*fid(1)^(2)+l(2)*m(4)*sin(fi(2))*fid(2)^(2)+l(3)*m(4)*sin(fi(3))*fid(3)^(2)+d(4)*m(4)*sin(fi(4))*fid(4)^(2)-Fr5y-(981*m(4))/100)+(d(1)*sin(fi(1))-sin(fi(1))*(d(1)-l(1)))*(l(1)*m(4)*cos(fi(1))*fid(1)^(2)+l(2)*m(4)*cos(fi(2))*fid(2)^(2)+l(3)*m(4)*cos(fi(3))*fid(3)^(2)+d(4)*m(4)*cos(fi(4))*fid(4)^(2)-Fr5x)-(d(1)*cos(fi(1))-cos(fi(1))*(d(1)-l(1)))*(l(1)*m(3)*sin(fi(1))*fid(1)^(2)+l(2)*m(3)*sin(fi(2))*fid(2)^(2)+d(3)*m(3)*sin(fi(3))*fid(3)^(2)-(981*m(3))/100)+(d(1)*sin(fi(1))-sin(fi(1))*(d(1)-l(1)))*(l(1)*m(3)*cos(fi(1))*fid(1)^(2)+l(2)*m(3)*cos(fi(2))*fid(2)^(2)+d(3)*m(3)*cos(fi(3))*fid(3)^(2))-(d(1)*cos(fi(1))-cos(fi(1))*(d(1)-l(1)))*(l(1)*m(2)*sin(fi(1))*fid(1)^(2)+d(2)*m(2)*sin(fi(2))*fid(2)^(2)-(981*m(2))/100)+d(1)*cos(fi(1))*(- d(1)*m(1)*sin(fi(1))*fid(1)^(2)+(981*m(1))/100)+d(1)^(2)*fid(1)^(2)*m(1)*cos(fi(1))*sin(fi(1));
                                                                                                                                       tor(1)-(d(2)*cos(fi(2))-cos(fi(2))*(d(2)-l(2)))*(l(1)*m(4)*sin(fi(1))*fid(1)^(2)+l(2)*m(4)*sin(fi(2))*fid(2)^(2)+l(3)*m(4)*sin(fi(3))*fid(3)^(2)+d(4)*m(4)*sin(fi(4))*fid(4)^(2)-Fr5y-(981*m(4))/100)+(d(2)*sin(fi(2))-sin(fi(2))*(d(2)-l(2)))*(l(1)*m(4)*cos(fi(1))*fid(1)^(2)+l(2)*m(4)*cos(fi(2))*fid(2)^(2)+l(3)*m(4)*cos(fi(3))*fid(3)^(2)+d(4)*m(4)*cos(fi(4))*fid(4)^(2)-Fr5x)-(d(2)*cos(fi(2))-cos(fi(2))*(d(2)-l(2)))*(l(1)*m(3)*sin(fi(1))*fid(1)^(2)+l(2)*m(3)*sin(fi(2))*fid(2)^(2)+d(3)*m(3)*sin(fi(3))*fid(3)^(2)-(981*m(3))/100)+(d(2)*sin(fi(2))-sin(fi(2))*(d(2)-l(2)))*(l(1)*m(3)*cos(fi(1))*fid(1)^(2)+l(2)*m(3)*cos(fi(2))*fid(2)^(2)+d(3)*m(3)*cos(fi(3))*fid(3)^(2))-d(2)*cos(fi(2))*(l(1)*m(2)*sin(fi(1))*fid(1)^(2)+d(2)*m(2)*sin(fi(2))*fid(2)^(2)-(981*m(2))/100)+d(2)*sin(fi(2))*(l(1)*m(2)*cos(fi(1))*fid(1)^(2)+d(2)*m(2)*cos(fi(2))*fid(2)^(2));
                                                                                                                                                                                                                                                                                                                         tor(2)-tor(1)-(d(3)*cos(fi(3))-cos(fi(3))*(d(3)-l(3)))*(l(1)*m(4)*sin(fi(1))*fid(1)^(2)+l(2)*m(4)*sin(fi(2))*fid(2)^(2)+l(3)*m(4)*sin(fi(3))*fid(3)^(2)+d(4)*m(4)*sin(fi(4))*fid(4)^(2)-Fr5y-(981*m(4))/100)+(d(3)*sin(fi(3))-sin(fi(3))*(d(3)-l(3)))*(l(1)*m(4)*cos(fi(1))*fid(1)^(2)+l(2)*m(4)*cos(fi(2))*fid(2)^(2)+l(3)*m(4)*cos(fi(3))*fid(3)^(2)+d(4)*m(4)*cos(fi(4))*fid(4)^(2)-Fr5x)-d(3)*cos(fi(3))*(l(1)*m(3)*sin(fi(1))*fid(1)^(2)+l(2)*m(3)*sin(fi(2))*fid(2)^(2)+d(3)*m(3)*sin(fi(3))*fid(3)^(2)-(981*m(3))/100)+d(3)*sin(fi(3))*(l(1)*m(3)*cos(fi(1))*fid(1)^(2)+l(2)*m(3)*cos(fi(2))*fid(2)^(2)+d(3)*m(3)*cos(fi(3))*fid(3)^(2));
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              tor(3)-tor(2)+Fr5x*sin(fi(4))*(d(4)-l(4))-d(4)*cos(fi(4))*(l(1)*m(4)*sin(fi(1))*fid(1)^(2)+l(2)*m(4)*sin(fi(2))*fid(2)^(2)+l(3)*m(4)*sin(fi(3))*fid(3)^(2)+d(4)*m(4)*sin(fi(4))*fid(4)^(2)-Fr5y-(981*m(4))/100)+d(4)*sin(fi(4))*(l(1)*m(4)*cos(fi(1))*fid(1)^(2)+l(2)*m(4)*cos(fi(2))*fid(2)^(2)+l(3)*m(4)*cos(fi(3))*fid(3)^(2)+d(4)*m(4)*cos(fi(4))*fid(4)^(2)-Fr5x)-Fr5y*cos(fi(4))*(d(4)-l(4))];
 
CasF.f_skeldyn = Function('f_skeldyn',{fi,fid,fidd,Fr5x,Fr5y,M,...
    Tmus},{As*fidd-bs},{'fi','fid','fidd','Fr5x','Fr5y','M','Tmus'},...
    {'Terr'});


% Hip position and velocity
rh = 0;
vh = 0;
for k=1:4
    rh = rh+l(k)*[cos(fi(k));sin(fi(k))];
    vh = vh+l(k)*fid(k)*[-sin(fi(k));cos(fi(k))];
end
CasF.f_Hip_kin = Function('f_Hip_kin', {fi,fid}, {rh, vh},...
    {'fi','fid'}, {'rh', 'vh'});


end