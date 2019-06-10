%Pure bending of beam with subdivision of edges
close all; clearvars; dbstop error;

probdef = @ProbDefPureBendingBeam; % function handler

s = cell(2,1);

mesher = 'RectMesh'; % keyword to access square grid of subdomains
para.h0 = 0.5; % size of square subdomains
para.nDiv = 1; % number of elements per edge
[ s{1}.coord, s{1}.sdConn, s{1}.sdSC ] = createSBFEMesh(probdef, mesher, para);
figure;
opt=struct('sdSC', s{1}.sdSC, 'LabelSC',1,   ...
    'PlotNode',1, 'MarkerSize', 6, 'LabelNode',11);
PlotSBFEMesh(s{1}.coord, s{1}.sdConn, opt);
title('MESH');
hold on


[ s{1}.meshEdge, s{1}.sdEdge, edge2sd, node2Edge, node2sd] = ...
    meshConnectivity( s{1}.sdConn );

hold on
PlotEdges( s{1}.coord, s{1}.meshEdge )


para.h0 = 0.2; % size of square subdomains
para.nDiv = 1; % number of elements per edge
[ s{2}.coord, s{2}.sdConn, s{2}.sdSC ] = createSBFEMesh(probdef, mesher, para);

cshft = [ 1.008 1]; % shift coordinates
nshft = length(s{1}.coord);
s{2}.coord = bsxfun(@plus,s{2}.coord, cshft);
s{2}.sdSC = bsxfun(@plus,s{2}.sdSC, cshft);% shift scaling centre

%scaling coordinates about a point
cscl = [ 0.5 0 0.5]; % [scalingFactor x- and y- coordiantes of orgin]
s{2}.coord = cscl(1)*bsxfun(@minus,s{2}.coord, cscl(2:3)); %scaling about a point
s{2}.coord = bsxfun(@plus,s{2}.coord, cscl(2:3)); % recover original position after scaling
s{2}.sdSC = cscl(1)*bsxfun(@minus,s{2}.sdSC, cscl(2:3));
s{2}.sdSC = bsxfun(@plus,s{2}.sdSC, cscl(2:3));

opt=struct('sdSC',s{2}.sdSC, 'LabelSC',1, ...
    'PlotNode',1, 'MarkerSize',6, 'LabelNode',11);
PlotSBFEMesh(s{2}.coord, s{2}.sdConn, opt);

[ s{2}.meshEdge, s{2}.sdEdge, edge2sd, node2Edge, node2sd] = ...
    meshConnectivity( s{2}.sdConn );

hold on
PlotEdges( s{2}.coord, s{2}.meshEdge )

eps =1d-3;
ctLine = cell(2,1);
for ii = 1:2
    edgeCentre = 0.5*(s{ii}.coord(s{ii}.meshEdge(:,1),:)+s{ii}.coord(s{ii}.meshEdge(:,2),:));
    ctLine{ii}.domain = ii;
    ctLine{ii}.edgeID = find(abs(edgeCentre(:,2)-0.5)< eps);
    ctLine{ii}.edge = s{ii}.meshEdge(ctLine{ii}.edgeID,:);
end

for ii = 1:2
    ctLine{ii}.node = unique(ctLine{ii}.edge(:));
    ctLine{ii}.ib = setdiff(ctLine{ii}.edge(:,1),ctLine{ii}.edge(:,2));
    ctLine{ii}.ie = setdiff(ctLine{ii}.edge(:,2),ctLine{ii}.edge(:,1));
    mdnode = setdiff(ctLine{ii}.node, [ctLine{ii}.ib ctLine{ii}.ie]);
    id = ctLine{ii}.domain;
    s{id}.coord(mdnode,:) = round( s{id}.coord(mdnode,:)/eps ) * eps;
    ctLine{ii}.xy = s{id}.coord(ctLine{ii}.node,:);
end


edgeXi{1} = InsertPointsOnEdge( ctLine{1}.edge, s{1}.coord, ctLine{2}.xy, eps );
edgeXi{2} = InsertPointsOnEdge( ctLine{2}.edge, s{2}.coord, ctLine{1}.xy, eps);

for ii = 1:2
    edgeDiv = cell(length(s{ii}.meshEdge),1);
    edgeDiv(ctLine{ii}.edgeID) = edgeXi{ii};
    n1 = length(s{ii}.coord);
    [s{ii}.coord, s{ii}.sdConn ] = subdivideEdge(edgeDiv, s{ii}.coord, s{ii}.meshEdge, s{ii}.sdEdge);
    ctLine{ii}.node = [ctLine{ii}.node; (n1+1:length(s{ii}.coord))'];
end

% assemble sub-structures
n1 = length(s{1}.coord);
s{1}.node = 1:n1;
n2 = n1 + length(s{2}.coord);
s{2}.node = n1+1:n2;
coord = [s{1}.coord; s{2}.coord];
n1 = length(s{1}.sdConn);
s{1}.sd = 1:n1;
n2 = n1 + length(s{2}.sdConn);
s{2}.sd = n1+1:n2;
sdSC = [s{1}.sdSC; s{2}.sdSC];
[s{2}.sdConn] = shift_sdConn(length(s{1}.coord), s{2}.sdConn);
sdConn = [ s{1}.sdConn; s{2}.sdConn];

%contact pair
ctLine{1}.node = s{1}.node(ctLine{1}.node);
ctLine{2}.node = s{2}.node(ctLine{2}.node);
for ii = 1:length(ctLine{1}.node)
    
end

figure;
opt=struct('sdSC', sdSC, 'LabelSC',12,   ...
    'fill', [0.9 0.9 0.9], 'PlotNode',1, 'LabelNode',12);
PlotSBFEMesh(coord, sdConn, opt);
title('MESH');


% %  analysis
% % [U, sdSln] = SBFEPoly2NSolver(probdef, coord, sdConn, sdSC);



