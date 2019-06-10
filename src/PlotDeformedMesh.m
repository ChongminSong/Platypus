function varargout = PlotDeformedMesh(d, coord, sdConn, opt)
%Plot deformed mesh
%
%Inputs:
%  coord(i,:)   - coordinates of node i
%  sdConn{isd,:}(ie,:)  - S-element conncetivity. The nodes of 
%                         line element ie in S-element isd. 
%  d            - nodal displacements
%
%  options:
%  opt=struct('MagnFct', 0.1, 'Undeformed',':b');
%  where the options are:
%    opt.MagnFct    : magnification factor of deformed mesh
%    opt.Undeformed : style of undeformed mesh

if nargin > 3 && isfield(opt, 'MagnFct')
   magnFct = opt.MagnFct; 
else
   magnFct = 0.1; %ltx default magnification factor
end

Umax = max(abs(d)); %ltx maximum displacement
Lmax = max(max(coord)-min(coord)); %ltx maximum dimension of domain
fct = magnFct*Lmax/Umax; %ltx factor to magnify the displacement 
%ltx augment nodal coordinates
deformed = coord + fct*(reshape(d,2,[]))'; 
hold on
%ltx plot undeformed mesh
if nargin > 3 && isfield(opt, 'Undeformed') && ...
        ~isempty(opt.Undeformed)
    %ltx plotting option of undeformed mesh
    undeformedopt = struct('LineSpec',opt.Undeformed); 
    PlotSBFEMesh(coord, sdConn, undeformedopt);
end
title('DEFORMED MESH')
%ltx plot deformed mesh
deformedopt = struct('LineSpec','-k'); %ltx plotting option
PlotSBFEMesh(deformed, sdConn, deformedopt);
varargout={deformed, fct};%ltx outputs