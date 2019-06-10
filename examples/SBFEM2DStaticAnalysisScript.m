%ltx {\bf Solution of S-elements and assemblage of global stiffness}
%ltx {\bf and mass matrices}
[sdSln, K, M] = SBFEMAssembly(coord, sdConn, sdSC, mat);

%ltx {\bf Static solution of nodal displacements and forces}
[d, F] = SolverStatics(K, BC_Disp, F);

%ltx {\bf Plot deformed mesh}
figure
opt = struct('MagnFct', 0.1, 'Undeformed','--k');
PlotDeformedMesh(d, coord, sdConn, opt)