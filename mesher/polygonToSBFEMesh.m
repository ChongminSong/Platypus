function [sdConn, sdSC] = polygonToSBFEMesh(coord, polygon)
%Convert a polygon mesh to an S-element mesh 
%
%Inputs:
%  coord(i,:)   - coordinates of node i
%  polygon{i}   - array of vertices of polygon i.
%
%Outputs:
%  sdConn{isd,:}(ie,:)  - S-element conncetivity. The nodes of 
%                         line element ie in S-element isd. 
%  sdSC(isd,:)  - coordinates of scaling centre of S-element isd

%cbgnltx
% \texttt{sdConn} : S-element connectivity stored as a cell array (One S-element per cell). 
% In a cell, the connectivity of line elements is given by
% one element per row \texttt{[Node-1 Node-2]}.
%cendltx
nsd = length(polygon); %ltx number of S-elements
sdConn = cell(nsd,1);  %ltx initialising connectivity
sdSC = zeros(nsd,2);   %ltx scaling centre
for isub =1:nsd
    %ltx build connectivity 
    sdConn{isub}=[polygon{isub}; ...
                  polygon{isub}(2:end) polygon{isub}(1)]';
    %ltx scaling centre at centroid of polygon
    sdSC(isub,:) = polygonCentroid(coord(polygon{isub},:));
end

end