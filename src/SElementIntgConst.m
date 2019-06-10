function [ sdIntgConst ] = SElementIntgConst( d, sdSln )
%Integration constants of S-elements
%
% Inputs:
%   d           - nodal displacements
%   sdSln       - solutions for S-elements 
% Outputs:
%   sdIntgConst - vector of intergration constants

Nsd = length(sdSln); %ltx total number of S-elements
sdIntgConst = cell(Nsd,1); %ltx initialisation of output argument
for isd = 1:Nsd %ltx loop over S-elements
    
    %ltx {\bf Integration constants}
    %ltx global DOFs of nodes in an S-element
    dof = reshape([2*sdSln{isd}.node-1;2*sdSln{isd}.node], [], 1);
    dsp = d(dof); %ltx nodal displacements at boundary of S-element 
    sdIntgConst{isd} = (sdSln{isd}.v)\dsp; %ltx integration constants, see Eq.~\eqref{eq:EigenMethod-internalSln-c}
   
end