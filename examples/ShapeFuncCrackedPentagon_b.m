%% Solution of pentagon subdomain

addpath('..\'); clearvars; close all;

%ltx {\bf Mesh}
%ltx nodal coordinates. One node per row \texttt{[x y]}
%coord = [cosd(-126:72:180+76);sind(-126:72:180+76)]';
coord = [cosd(-126:72:240);sind(-126:72:240)]';
%cbgnltx
% Input subdomain connectivity as a cell array (One subdomain per cell).
% In a cell, the connectivity of line elements is given by
% one element per row \texttt{[Node-1 Node-2]}.
%cendltx
sdConn = { [1:5 ; 2:6]' }; %ltx subdomain 1
%ltx coordinates of scaling centres of subdomains.
sdSC = [ 0 0]; %ltx one subdomain per row

%ltx {\bf Materials}
%ltx \texttt{E}: Young's modulus; \texttt{p}: Poisson's ratio
ElasMtrx = @(E, p) E/(1-p^2)*[1 p 0;p 1 0;0 0 (1-p)/2];
D = ElasMtrx(10E6, 0.25); %ltx  \texttt{E} in KPa
den = 2000; %ltx mass density


%ltx {\bf Global stiffness matrix}
[sdSln, K, M] = SBFEM_Assemblage(coord, sdConn, sdSC, D, den);

nd = size(K,1);
U = zeros(nd,1);
U(8) = 1;
% integration constants and strain modes
sdPstP = SubdomainPostP2NodeEle( U, sdSln );

isd = 1; % subdomain number
nNode = size(coord,1);
xi = [1:-0.02:0].^2; %radial coordinate
X = zeros(length(xi), nNode); Y = X; Z = X; 
% displacements and strains at the specified raidal coordinate
for ii= 1:length(xi)
    [nodexy, dsp, strnNode, GPxy, strnEle] = ...
        SubdomainInDispStrain(xi(ii), sdSln{isd}, sdPstP{isd});
    X(ii,:) =  [nodexy(:,1)' ];
    Y(ii,:) =  [nodexy(:,2)' ];
    Z(ii,:) =  [dsp(2:2:end)'];
end

figure(1)
surf(X,Y,Z,'FaceColor','interp',...
   'EdgeColor','none',  ...
   'FaceLighting','phong');
axis equal; xlabel('x'); ylabel('y'); zlabel('N');
colorbar;

figure(2)
contourf(X,Y,Z, 20, 'LineStyle','none');
axis equal; colorbar;
xlabel('x'); ylabel('x');
