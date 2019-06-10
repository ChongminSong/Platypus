function [x] = ProbDefCombinedPlateWithHole(Demand,Arg)
%Plate with a hole composed of two domains
%
%inputs
%   Demand    - option to access data in local functions
%   Arg       - argument to pass to local functions
%output
%   x         - cell array of output. 

switch(Demand) %ltx access local functions according to \texttt{Demand}
    case('Domains'); x = definitionOfDomains();
    case('MAT');     x = ProbDefInfPlateWithHole('MAT'); 
    case('BCond');   x = BoundryConditions(Arg);
    case('BCondModal');   x = BoundryConditionsModal(Arg);
    case('Output');  x = [];
    case('EXACT');   x =  ExactSln(Arg);
end
end

%ltx {\bf  Definition of composition of domains}\label{ProbDefCombinedPlateWithHoleDomain}
function [x] = definitionOfDomains ()

%ltx handle of problem definition functions of  constituent domains
prodef ={@ProbDefPureBendingBeam, @ProbDefInfPlateWithHole}; 
cshft = [ -1 -1.5;  0 0 ]; %ltx amount of shifting of x- and y- coordinates

%ltx Definition of the interface between the two domains
onLine = {@(xy, eps) find(abs(xy(:,2)+1) < eps)}; %ltx lines $y=-1$

x = struct('prodef',{prodef}, 'cshft',cshft, 'onLine',{onLine}); %ltx output
end

%ltx {\bf Prescribe boundary conditions}\label{ProbDefCombinedPlateWithHoleBC}
function [x] = BoundryConditions(Arg)
coord = Arg{1}; %ltx nodal coordinates

%ltx displacement constrains (prescribed acceleration in a response
%ltx  history analysis). One constrain per row: \texttt{[Node Dir]}
eps = 1d-4*(max(coord(:,1))-min(coord(:,1)));
bNodes = find(abs(coord(:,1)-min(coord(:,1)))<eps | ...
    abs(coord(:,1)-max(coord(:,1)))<eps | ...
    abs(coord(:,2)-min(coord(:,2)))<eps | ...
    abs(coord(:,2)-max(coord(:,2)))<eps );
%ltx get exact solution at specified nodes
ex = ProbDefInfPlateWithHole('EXACT', coord(bNodes,:)); 
BC_Disp = [ [bNodes; bNodes ], [ones(length(bNodes),1);...
    2*ones(length(bNodes),1)], [ex(:,1); ex(:,2)] ]; 

%ltx nodal forces
F = zeros(2*size(coord,1),1); %ltx initialising force vector

%ltx output
x = {BC_Disp, F};
end

%ltx {\bf Prescribe boundary conditions}\label{ProbDefCombinedPlateWithHoleBCModal}
function [x] = BoundryConditionsModal(Arg)
coord = Arg{1}; %ltx nodal coordinates

%ltx displacement constrains (prescribed acceleration in a response
%ltx  history analysis). One constrain per row: \texttt{[Node Dir]}
eps = 1d-4*(max(coord(:,1))-min(coord(:,1)));
bNodes = find(abs(coord(:,2)-min(coord(:,2)))<eps);
BC_Disp = [ [bNodes; bNodes ], [ones(length(bNodes),1);...
    2*ones(length(bNodes),1)], zeros(2*length(bNodes),1) ]; 

%ltx output
x = {BC_Disp};
end

function x = ExactSln(xy) %ltx exact solution \label{ProbDefCombinedPlateWithHolesln}
x = ProbDefInfPlateWithHole('EXACT', xy);
end