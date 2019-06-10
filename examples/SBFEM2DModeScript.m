%ltx {\bf Solution of subdomains and assemblage of global stiffness}
%ltx {\bf and mass matrices}
[sdSln, K, M] = SBFEMAssemblage(coord, sdConn, sdSC, mat);

%ltx {\bf Compute natrural frequencies and modal shapes}

if ~exist('NMode','var')
    NMode = 1; %ltx number of modes to extract
end
if ~exist('PltMode','var')
    PltMode = 1; %ltx number of modes to plot
end

[freqAngular, modeShape] = SolverMode(NMode, K, M, BC_Disp);

%ltx  output natrual frequency (radians/second) and period (seconds)
[ (1:NMode)' freqAngular 2*pi./freqAngular]

%ltx plot modal shape
for imode = PltMode; %ltx mode number
    figure
    U = modeShape(:,imode); %ltx set the modal shape as displacement
    Umax = max(abs(U)); %ltx maximum displacement
    Lmax = max(max(coord)-min(coord)); %ltx maximum dimension of domain
    fct = 0.1*Lmax/Umax; %ltx factor to magnify the displacement 
    %ltx augment nodal coordinates
    deformed = coord + fct*(reshape(U,2,[]))';
    hold on
    opt = struct('LineSpec',':b'); %ltx plotting option
    PlotSBFEMesh(coord, sdConn, opt); %ltx undeformed mesh
    opt = struct('LineSpec','-k'); %ltx plotting option
    PlotSBFEMesh(deformed, sdConn, opt); %ltx modal shape
    title(['MODE SHAPE: ',int2str(imode)]);
end