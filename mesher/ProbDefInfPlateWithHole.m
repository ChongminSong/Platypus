function [x] = ProbDefInfPlateWithHole(Demand,Arg)
%Problem definition function of a circular hole in
% an infinite plane under remote uniaxial tension
%
%inputs
%   Demand    - option to access data in local functions
%   Arg       - argument to pass to local functions
%output
%   x         - cell array of output.

%ltx \texttt{L}: Length; \texttt{D}: Height; \texttt{R}: Radius of hole;
%ltx \texttt{E}: Young's modulus; \texttt{pos}: Poisson's ratio; \texttt{den}: mass density
para = struct('L',2, 'a', 0.4, ...
    'E', 1E3, 'pos',0.25, 'den', 2, 'trac', 1);
BdBox = [-para.L/2 para.L/2 -para.L/2 para.L/2]; %ltx $[x_1,x_2,y_1,y_2]$
switch(Demand)
    case('para');    x = para; %ltx parameters of problem definitions
    %ltx direct input triangular mesh
    case('TriMesh'); x = [];
       
    %ltx for use with polygon mesh generator PolyMesher 
    case('Dist');  x = DistFnc(Arg,para,BdBox);
    case('BC');    x = {[],[]};
        
    %ltx for use with triangular mesh generator DistMesh 
    case('Dist_DistMesh');   x = DistFnc(Arg,para,BdBox);
        x = x(:,end); %ltx only use the last column
    case('fh');
        x = huniform(Arg); %ltx element size \label{Code_ProbDefInfPlatWithHole_elemSize}
        % %         x = 1+1*dcircle(Arg,0,0,para.a);
    case('pfix');
        x = [BdBox(1) BdBox(3); BdBox(1) BdBox(4); ...
            BdBox(2) BdBox(3); BdBox(2) BdBox(4)];
        
    %ltx for use with PolyMesher and DistMesh 
    case('BdBox');    x = BdBox;
        
    case('MAT'); 
        x = struct('D',IsoElasMtrx(para.E, para.pos), ...
            'den', para.den);  %ltx material  constants
        
    case('BCond');  x = BoundryConditions(Arg, para, BdBox);
        
    case('Output'); x = OutputRequests(Arg, BdBox);
        
    case('EXACT');  x =  ExactSln(Arg, para);
end
end

%ltx \textbf{Signed distance function}\label{Code_ProbDefInfPlatWithHole_DistFnc}
function Dist = DistFnc(P,plate,BdBox)
d1 = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
d2 = dCircle(P, 0, 0, plate.a); %ltx Centre and radius of circle $x_c$, $y_c$, $r$
Dist = dDiff(d1,d2);
end

%ltx \textbf{Boundary conditions}  \label{Code_ProbDefInfPlatWithHole_BC}
function [x] = BoundryConditions(Arg, para, BdBox)
coord = Arg{1}; %ltx nodal coordinates
eps = 0.01*sqrt((BdBox(2)-BdBox(1))*...
    (BdBox(4)-BdBox(3))/size(coord,1)); %ltx tolerance
%ltx displacement constrains (prescribed acceleration in a response
%ltx history analysis). One constrain per row: \texttt{[Node Dir]}
d = dRectangle(coord,BdBox(1),BdBox(2), ...
    BdBox(3),BdBox(4)); %ltx distance from nodes to the square boundary
bNodes = find(abs(d(:,end))<eps); %ltx nodes on the square boundary
% % wNodes = find(abs(coord(:,1)-BdBox(1))<eps);
% % eNodes = find(abs(coord(:,1)-BdBox(2))<eps);
% % sNodes = find(abs(coord(:,2)-BdBox(3))<eps);
% % nNodes = find(abs(coord(:,2)-BdBox(4))<eps);
% % bNodes = unique([wNodes; eNodes; sNodes; nNodes]);
ex = ExactSln(coord(bNodes,:), para); %ltx exact solution
n1 = length(bNodes);
BC_Disp = [ [bNodes; bNodes], [ones(n1,1); 2*ones(n1,1)], ...
    [ex(:,1); ex(:,2)]]; %ltx displacement boundary condition

%ltx nodal forces
F = zeros(2*size(coord,1),1); %ltx initialising force vector

%ltx output
x = {BC_Disp, F};
end

%ltx {\bf Output Requests}
function [outputs] = OutputRequests(Arg, BdBox)
outputs = {};
end

%ltx \textbf{Exact solution} 
function x = ExactSln(xy, para)
p = para.trac; %ltx remote tension
a = para.a; %ltx radius of hole
pos = para.pos; G = para.E/(2*(1+pos)); %ltx material constants
x = xy(:,1);  y = xy(:,2);
ar = a./sqrt(x.^2 + y.^2); c = atan2(y,x); %ltx polar coordinates
%ltx displacements (see Eq.~\eqref{eq:ProbDefInfPlateWithHoleExactDisp})
ka = (3-pos)/(1+pos); %ltx Kolosov constant for plane stress
ux = p*a/8/G *( 1./ar*(1+ka).*cos(c) + ...
    2*ar.*((1+ka)*cos(c)+cos(3*c)) - 2*ar.^3.*cos(3*c) );
uy = p*a/8/G *( 1./ar*(ka-3).*sin(c) + ...
    2*ar.*((1-ka)*sin(c)+sin(3*c)) - 2*ar.^3.*sin(3*c) );
%ltx stresses (see Eq.~\eqref{eq:ProbDefInfPlateWithHoleExactStress})
sx  = p/2*( 2 - ar.^2.*(3*cos(2*c) + (2-3*ar.^2).*cos(4*c)) );
sy  = -p*ar.^2/2.*(cos(2*c) - (2-3*ar.^2).*cos(4*c));
sxy = -p*ar.^2/2.*(sin(2*c) + (2-3*ar.^2.).*sin(4*c));
x = [ux, uy, sx, sy, sxy]; %ltx output
end
