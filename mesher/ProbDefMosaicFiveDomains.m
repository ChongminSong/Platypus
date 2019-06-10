function [x] = ProbDefMosaicFiveDomains(Demand,Arg)
%Problem definition function of a mosaic composed of five domains
%
%inputs
%   Demand    - option to access data in local functions
%   Arg       - argument to pass to local functions
%output
%   x         - cell array of output.

switch(Demand) %ltx access local functions according to \texttt{Demand}
    case('Domains'); x = definitionOfDomains();
    case('MAT'); x = ProbDefInfPlateWithHole('MAT'); %ltx material
    case('BCond');  x = BoundryConditions(Arg);
% %     case('BCondModal');  x = BoundryConditions(Arg);
    case('Output'); x = [];
    case('EXACT');  x = [];
end
end

%ltx {\bf Definition of composition of domains}\label{ProbDefMosaicFiveDomainsDomain}
function [x] = definitionOfDomains ()

prodef ={@ProbDefLShapedPlate, @ProbDefInfPlateWithHole, ...
    @ProbDefInfPlateWithHole,  @ProbDefLShapedPlate, ...
    @ProbDefPureBendingBeam};
%ltx amount of shifting of x- and y- coordinates
cshft = [-1 -1; 2 0;  0.5 2; 1.5 1; 0.5 3.5]; 

%ltx Definitions of the interfaces between the domains
onLine = {@(xy, eps) find(abs(xy(:,1)-1) < eps), ... %ltx line $x=1$
    @(xy, eps) find(abs(xy(:,2)-1) < eps), ... %ltx line $y=1$
    @(xy, eps) find( (abs(xy(:,2)-1) < eps) | ...
    (abs(xy(:,1)-1.5) < eps)), ... %ltx lines $y=1$ and $x=1.5$
    @(xy, eps) find(abs(xy(:,2)-3) < eps)}; ... %ltx line $y=3$

x = struct('prodef',{prodef}, 'cshft',cshft, ...
    'onLine',{onLine}); %ltx output
end

%ltx {\bf Prescribe boundary conditions}\label{ProbDefMosaicFiveDomainsBC}
function [x] = BoundryConditions(Arg)
coord = Arg{1}; %ltx nodal coordinates
sdConn = Arg{2}; %ltx element connectivity

%ltx displacement constrains (prescribed acceleration in a response
%ltx  history analysis). One constrain per row: \texttt{[Node Dir]}
eps = 0.001*(max(coord(:,1))-min(coord(:,1)));
bNodes = find(abs(coord(:,2)-min(coord(:,2)))<eps);
BC_Disp = [ [bNodes; bNodes], [ones(length(bNodes),1);...
    2*ones(length(bNodes),1)], 0*[bNodes; bNodes] ];

%ltx surface traction
%ltx edges of polygon S-elements
[ MeshEdges ] = meshConnectivity( sdConn );
%ltx find edges at the top
ymax = max(coord(:,2));
tEdges = MeshEdges((abs(coord(MeshEdges(:,1),2)-ymax)+ ...
    abs(coord(MeshEdges(:,2),2)-ymax))<eps,:);
%ltx nodal forces
F = zeros(2*size(coord,1),1); %ltx initialising force vector
%ltx add nodal forces equivalent to surface traction
F = addSurfTraction(coord, tEdges, [0; 1; 0; 1], F);

%ltx output
x = {BC_Disp, F};
end