function [ F ] = addSurfTraction(coord, edge, trac, F)
%Assembly of surface tractions as equivalent nodal forces
% to load vector
%
%Inputs:
%  coord(i,:)   - coordinates of node i
%  edge(i,1:2)  - the 2 nodes of edge i
%  trac         - surface traction
%     when trac has only one column
%       trac(1:4) - surface traction at the 2 nodes of all edges
%     when trac has more than one column
%       trac(1:4,i) - surface traction at the 2 nodes of edge i
%
%Outputs:
%  F            - global load vector 

if size(trac,2) == 1 %ltx expand uniform surface traction to all edges
    trac = trac(:, ones(1,length(edge)));
end

%ltx equivalent nodal forces
fmtx = [ 2 0 1 0; 0 2 0 1; 1 0 2 0; 0 1 0 2]; %ltx see Eq.~\eqref{eq:2D-2NodeElem-trac-NodalForce}
edgeLen = sqrt(sum( ( coord(edge(:,2),:) - ...
                     coord(edge(:,1),:) ).^2, 2 )); %ltx edge length
nodalF = 1/6*fmtx*trac.*edgeLen(:,ones(1,4))'; %ltx Eq.~\eqref{eq:2D-2NodeElem-trac-NodalForce}

%ltx assembly of nodal forces
for ii = 1:size(edge,1)
    dofs = [2*edge(ii,1)-1 2*edge(ii,1) ...
            2*edge(ii,2)-1 2*edge(ii,2)];
    F(dofs) = F(dofs) + nodalF(:,ii);
end

end