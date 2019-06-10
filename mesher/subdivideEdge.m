function [coord, sdConn ] = subdivideEdge(edgeXi, ...
                            coord,  meshEdge, sdEdge )
%Subdivide selected edges of S-elements into multiple elements
%
%Inputs
%   edgeXi(i)        - When it is an array: number of subdivision
%                           of edge i
%                    - When it is a cell array: column vector of   
%                      local parametric coordinate (between 0  
%                      and 1) of subdivision of edge i
%   coord(i,:)       - coordinates of node i
%   meshEdge(i,1:2)  - the 2 nodes on line i in a mesh.  
%                      The node number of the first node is 
%                      larger than that of the 2nd one
%   sdEdge{i}        - lines forming S-element i, 
%                      >0 when a line follows anti-clockwise  
%                         direction around the scaling centre.   
%                      <0 clockwise direction 
%
%Outpus
%   coord  - nodal coordinates.  New nodes are padded to the end
%   sdConn{isd,:}(ie,:)  - S-element conncetivity. The nodes of
%                         line element ie in S-element isd.

nEdge0 = length(meshEdge); %ltx number of edges in original mesh
nExtraPnt = 0; %ltx number of extra points after subdivision
%ltx local parametric coordinate $\xi$ (between 0 and 1) of subdivision
if ~iscell(edgeXi) 
    edgeDiv = edgeXi; %ltx save array of number of subdivision of edges 
    edgeXi = cell(nEdge0,1);
    for ii = 1:nEdge0
        if edgeDiv(ii) > 1
            nDiv = edgeDiv(ii); %ltx number of subdivision
            edgeXi{ii} = (1:1:nDiv-1)'/nDiv; %ltx uniform subdivision
            nExtraPnt = nExtraPnt + nDiv-1;
        end
    end
else
    nExtraPnt = sum( cellfun(@length, edgeXi) );
end

%ltx create new nodes and mesh edges
newEdges = cell(nEdge0,1); %ltx to store edges after division
ib = length(coord) + 1; %ltx index for  new nodes
coord = [coord; zeros(nExtraPnt, 2)]; %ltx including new nodes
for ii = 1:nEdge0 %ltx dividing edges in orginial mesh
    xi = edgeXi{ii}; %ltx local parametric coordinate $\xi$ (0, 1)
    if ~isempty(xi)
        ie = ib + length(xi) - 1; %ltx index for the last new node
        xyb = coord(meshEdge(ii,1),:); %ltx coordinates of start point
        xye = coord(meshEdge(ii,2),:); %ltx coordinates of end point
        coord(ib:ie,:) = bsxfun(@times,xyb,1-xi) + ...
            bsxfun(@times,xye,xi); %ltx new nodes by interpolation
        newEdges{ii} = [meshEdge(ii,1), ib:ie; ...
            ib:ie meshEdge(ii,2)]'; %ltx new edges on an old one
        ib = ie + 1;
    else
        newEdges{ii} = meshEdge(ii,:); %ltx no division
    end
end

%ltx construct S-element connectivity
nsd = length(sdEdge); %ltx number of S-elements
sdConn = cell(nsd,1);  %ltx initialising connectivity
for isd = 1:nsd
    sdNewEdges = newEdges(abs(sdEdge{isd})); %ltx update connectivity
    %ltx find edges that follow clockwise direction
    idx = find(sdEdge{isd}(:) < 0);
    for ii = idx'
        %ltx reverse the element orientation
        sdNewEdges{ii} = [sdNewEdges{ii}(end:-1:1,2), ...
            sdNewEdges{ii}(end:-1:1,1)];
    end
    %ltx new element connectivity of an S-element
    sdConn{isd} = vertcat(sdNewEdges{:}); 
end

end