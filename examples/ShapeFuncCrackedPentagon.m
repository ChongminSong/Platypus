%% Solution of pentagon subdomain

addpath('..\'); clearvars; close all;
dbstop if error

%ltx {\bf Mesh}
%ltx nodal coordinates. One node per row \texttt{[x y]}
coord = [cosd(-126:72:180);sind(-126:72:180)]';
%coord = [cosd(18:72:360);sind(18:72:360)]';
conn = [1:5 ; 2:5 1]';

%ltx insert crack
crkedNode = 3;
dspNode = 5;
coord = [coord(1:crkedNode,:); ...
    coord(crkedNode,:); coord(crkedNode+1:end,:)]; 
conn(conn>crkedNode) = conn(conn>crkedNode)+1;
conn(crkedNode,1) = crkedNode + 1;
conn = [conn(crkedNode:end, :); conn(1:crkedNode-1, :)];
nodes = [conn(:,1) ; conn(end,2)];
sdConn = { conn };

sdSC = [ 0 0]; %ltx scaling centre

%ltx {\bf Materials}
%ltx \texttt{E}: Young's modulus; \texttt{p}: Poisson's ratio
ElasMtrx = @(E, p) E/(1-p^2)*[1 p 0;p 1 0;0 0 (1-p)/2];
D = ElasMtrx(10E6, 0.25); %ltx  \texttt{E} in KPa
den = 2000; %ltx mass density


%ltx {\bf Global stiffness matrix}
[sdSln, K, M] = SBFEM_Assemblage(coord, sdConn, sdSC, D, den);

nd = size(K,1);
U = zeros(nd,1);
U(2*dspNode) = 1;
%ltx integration constants and strain modes
sdPstP = SubdomainPostP2NodeEle( U, sdSln );

isd = 1; %ltx subdomain number
nNode = size(coord,1);
xi = (1:-0.05:0).^2; %ltx radial coordinate
X = zeros(length(xi), nNode); Y = X; Z = X; 
%ltx displacements and strains at the specified raidal coordinate
for ii= 1:length(xi)
    [nodexy, dsp, strnNode, GPxy, strnEle] = ...
        SubdomainInDispStrain(xi(ii), sdSln{isd}, sdPstP{isd});
    X(ii,:) =  (nodexy(nodes,1)');
    Y(ii,:) =  (nodexy(nodes,2)');
    Z(ii,:) =  (dsp(2*nodes)');
end

figure('Color','white');
surf(X,Y,Z,'FaceColor','interp',... 
           'EdgeColor','none',  ...
           'FaceLighting','phong');
daspect([1 1 2]);
view(40, 30);
hold on
text(1.1*(coord(nodes,1)-0.02), 1.1*coord(nodes,2), Z(1,1:end)'+0.05,int2str(nodes));
plot3(X(1,:), Y(1,:), Z(1,:), '-k'); %ltx edges
axis off; xlabel('x'); ylabel('x'); zlabel('N');
%ltx gray scale shading
invGray = 0.95*(1-gray); 
colormap(invGray)

% figure
% contourf(X,Y,Z, 20, 'LineStyle','none');
% hold on
% text(1.05*(coord(nodes,1)-0.02), 1.05*coord(nodes,2),int2str(nodes));
% plot(X(1,:), Y(1,:), '-ko');
% axis off; xlabel('x'); ylabel('x');
% colormap(invGray)
% 
% % figure(1)
% % surf(X,Y,Z,'FaceColor','interp',...
% %    'EdgeColor','none',  ...
% %    'FaceLighting','phong');
% % axis equal; xlabel('x'); ylabel('y'); zlabel('N');
% % colorbar;
% % 
% % figure(2)
% % contourf(X,Y,Z, 20, 'LineStyle','none');
% % axis equal; colorbar;
% % xlabel('x'); ylabel('x');
