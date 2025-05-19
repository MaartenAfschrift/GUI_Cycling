function [rJoint] = getJointPositions(parms,fi)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


[nseg, nfr] = size(fi);
if nseg>nfr
    fi = fi';
    [nseg, nfr] = size(fi);
end

% loop over all frames
rJoint = zeros(2,nseg+1,nfr);
for i = 1:nfr
    for k=1:nseg
        r_prox = squeeze(rJoint(:,k,i));
        rJoint(:,k+1,i) = r_prox+parms.segparms.L(k)*[cos(fi(k,i));sin(fi(k,i))];
    end
end


       
end