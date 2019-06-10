%Combine 2 rectangular domains 
close all; clearvars; dbstop error;

probdef = @ProbDefTwoRectDomains; % function handler
mesh = probdef('MESH', []); %ltx get mesh
coord = mesh{1}; sdConn = mesh{2}; sdSC = mesh{3};

% d = cell(2,1);
% mesher = 1;
% para = struct('h0',0.25, 'nPolygon',32, 'nDiv',2);
% [ d{1}.coord, d{1}.sdConn, d{1}.sdSC ] = ...
%     createSBFEMesh(@ProbDefPureBendingBeam, mesher, para);
% [ d{1}.meshEdge, d{1}.sdEdge] = meshConnectivity( d{1}.sdConn );
% 
% mesher = 1;
% para = struct('h0',0.2, 'nPolygon',64, 'nDiv',2);
% [ d{2}.coord, d{2}.sdConn, d{2}.sdSC ] = ...
%     createSBFEMesh(@ProbDefPureBendingBeam, mesher, para);
% 
% %scaling coordinates about a point
% cscl = [ 0.4 0 -0.5]; % [scalingFactor x- and y- coordiantes of orgin]
% d{2}.coord = cscl(1)*bsxfun(@minus,d{2}.coord, cscl(2:3)); %scaling about a point
% d{2}.coord = bsxfun(@plus,d{2}.coord, cscl(2:3)); % recover original position after scaling
% d{2}.sdSC = cscl(1)*bsxfun(@minus,d{2}.sdSC, cscl(2:3));
% d{2}.sdSC = bsxfun(@plus,d{2}.sdSC, cscl(2:3));
% 
% cshft = [ 0.5 1]; % shift coordinates
% d{2}.coord = bsxfun(@plus,d{2}.coord, cshft);
% d{2}.sdSC = bsxfun(@plus,d{2}.sdSC, cshft);% shift scaling centre
% 
% [ d{2}.meshEdge, d{2}.sdEdge] = meshConnectivity( d{2}.sdConn );
% 
% eps = 1d-4;
% onLine  = @(x, eps) find(abs(x(:,2)-0.5) < eps);
% [ coord, sdConn, sdSC ] = combineDomains( d, onLine, eps );

figure; % plot mesh
opt=struct('sdSC', sdSC, 'PlotNode',1,'MarkerSize', 3);
PlotSBFEMesh(coord, sdConn, opt); 
title('MESH');

% analysis
[U, sdSln] = SBFEPoly2NSolver(probdef, coord, sdConn, sdSC);