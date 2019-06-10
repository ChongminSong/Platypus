function [ F ] = AddNodalForces(BC_Frc, F)
%Assembly of prescribed nodal forces to load vector

%Inputs:
%  BC_Frc(i,:)   - one force component per row [Node Dir F]
%  F             - nodal force vector
%
%Output:
%  F             - nodal force vector

ndn = 2; %ltx 2 DOFs per node
if ~isempty(BC_Frc)
    fdof = (BC_Frc(:,1)-1)*ndn + BC_Frc(:,2); %ltx DOFs 
    F(fdof) = F(fdof) + BC_Frc(:,3); %ltx accumulate forces
end

end