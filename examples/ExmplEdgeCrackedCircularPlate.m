%cbgnltx
%\leftadjust  {\bf Matlab script: Edge-cracked circular body under uniform radial tension}  
%cendltx
clearvars; close all; dbstop if error

%ltx {\bf Input} 
R = 1; %ltx radius of circular body
p = 1000; %ltx radial surface traction (KPa)
E = 10E6; %ltx  \texttt{E} in KPa
nu = 0.25; %ltx Poisson's ratio
den = 2; %ltx mass density in $\unit{Mg/m^{3}}$

a = 0.75*R; %ltx crack length\label{Exmpl_EdgeCrackedCircularPlate_crackLength}
nq = 4; %ltx number of elements on one quadrant

%ltx {\bf Mesh}
%ltx nodal coordinates. One node per row \texttt{[x y]}
n = 4*nq; %ltx number of element on boundary
dangle = 2*pi/n; %ltx angular increment of an element $\Delta_\theta$
angle = -pi:dangle:pi+dangle/5; %ltx angular coordinates of nodes
%ltx Note: there are 2 nodes at crack mouth with the same coordinates
%ltx nodal coordinates $ x=R\cos(\theta)$, $ y=R\sin(\theta)$ 
coord = R*[cos(angle);sin(angle)]';
%cbgnltx
% Input S-element connectivity as a cell array (One S-element per cell). 
% In a cell, the connectivity of line elements is given by
% one element per row \texttt{[Node-1 Node-2]}.
%cendltx
sdConn = { [1:n; 2:n+1]' }; %ltx Note: elements form an open loop
%ltx select scaling centre at crack tip
sdSC = [a-R 0];

%ltx {\bf Materials}
%ltx \texttt{E}: Young's modulus; \texttt{p}: Poisson's ratio
ElasMtrx = @(E, p) E/(1-p^2)*[1 p 0;p 1 0;0 0 (1-p)/2];
mat.D = ElasMtrx(E, nu); 
mat.den = den; 

%ltx {\bf Boundary conditions} 
%ltx displacement constrains. One constrain per row: \texttt{[Node Dir Disp]}
BC_Disp = [  nq+1 1  0; 2*nq+1 2  0; ...
            3*nq+1 1  0]; %ltx constrain rigid-body motion
eleF = R*dangle*p; %ltx resultant radial force on an element $p R\Delta_\theta$
%ltx assemble radial nodal forces $\{F_r\}$
nodalF = [eleF/2, eleF*ones(1,n-1), eleF/2]; 
%ltx nodal forces. One node per row: \texttt{[Node Dir F]} 
       %ltx $F_x = F_r \cos(\theta)$, $F_y = F_r \sin(\theta)$
BC_Frc = [1:n+1 1:n+1; ones(1,n+1) 2*ones(1,n+1); ...
            nodalF.*cos(angle) nodalF.*sin(angle)]'; 

%ltx {\bf Plot mesh}
h1 =  figure; 
opt=struct('LineSpec','-k', 'sdSC',sdSC, ...
    'PlotNode',1, 'LabelNode', 1,...
     'BC_Disp',BC_Disp); %ltx plotting options
PlotSBFEMesh(coord, sdConn, opt);
title('MESH');

%ltx solution of S-elements and global stiffness and mass matrices
[sdSln, K, M] = SBFEMAssembly(coord, sdConn, sdSC, mat);

%ltx{\bf Assemblage external forces} 
ndn = 2; %ltx 2 DOFs per node
NDof = ndn*size(coord,1); %ltx number of DOFs
F = zeros(NDof,1); %ltx initialising right-hand side of equation $[K]{u} = {F}$
F = AddNodalForces(BC_Frc, F); %ltx add prescribed nodal forces 

%ltx {\bf Static solution}
[U, F] = SolverStatics(K, BC_Disp, F);

CODnorm = (U(end)-U(2))*E/(p*R);
disp(['Normalised crack openning displacement = ', num2str(CODnorm)])    

%ltx plot deformed shape
Umax = max(abs(U)); %ltx maximum displacement
fct = 0.2/Umax; %ltx factor to magnify the displacement to 0.2 m
%ltx argument nodal coordinates
deformed = coord + fct*(reshape(U,2,[]))'; 
hold on
%ltx plotting options
opt = struct('LineSpec','-ro', 'LabelNode', 1); 
PlotSBFEMesh(deformed, sdConn, opt);
title('DEFORMED MESH');

%ltx {\bf Internal displacements and stresses}\label{Exmple_EdgeCrackedCircularPlate_internal}

%ltx strain modes of S-elements
sdStrnMode =  SElementStrainMode2NodeEle( sdSln );
%ltx integration constants of S-element 
sdIntgConst = SElementIntgConst( U, sdSln );

isd = 1; %ltx S-element number
xi = (1:-0.01:0).^2; %ltx radial coordinates \label{Exmple_EdgeCrackedCircularPlate_xi}
%ltx initialisation of variables for plotting
X = zeros(length(xi), length(sdSln{isd}.node));
Y = X; C = X;
%ltx displacements and strains at the specified raidal coordinate
for ii= 1:length(xi)
    [nodexy, dsp, strnNode, GPxy, strnEle] = ...
        SElementInDispStrain(xi(ii), sdSln{isd},  ...
                              sdStrnMode{isd}, sdIntgConst{isd});
    deformed = nodexy + fct*(reshape(dsp,2,[]))';
    %ltx coordinates of grid points
    X(ii,:) =  deformed(:,1)';
    Y(ii,:) =  deformed(:,2)';
    strsNode = mat.D*strnNode; %ltx nodal stresses
    C(ii,:) =  strsNode(2,:); %ltx store $\sigma_{yy}$ for plotting
end

%ltx plot stress contour
h2 = figure('Color','white');
contourf(X,Y,C, (-1000:1000:10000), 'LineStyle','none');
hold on
axis off; axis equal;
xlabel('x'); ylabel('x');
% % %ltx gray scale plot for inclusion in book. Remove for colour plotting
% % colormap( 0.9*(1-gray)); 
colormap(jet);colorbar;