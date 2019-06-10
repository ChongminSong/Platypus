%cbgnltx
%\leftadjust  {\bf Matlab script: A shape function of a cracked pentagon}  
%cendltx
clearvars; close all;
dbstop if error

%ltx {\bf Mesh}
%ltx nodal coordinates. One node per row \texttt{[x y]}
coord = [ -0.5878   -0.8090;    0.5878   -0.8090;
           0.9511    0.3090;    0.9511    0.3090;
                0    1.0000;   -0.9511    0.3090];
%ltx connectivity
sdConn = { [  4     5;    5     6;
              6     1;    1     2;
              2     3 ] };
sdSC = [ 0.2 0]; %ltx scaling centre

%ltx {\bf Materials}
%ltx \texttt{E}: Young's modulus; \texttt{p}: Poisson's ratio
ElasMtrx = @(E, p) E/(1-p^2)*[1 p 0;p 1 0;0 0 (1-p)/2];
D = ElasMtrx(10E6, 0.25); %ltx  \texttt{E} in KPa
den = 2000; %ltx mass density

%ltx {\bf Global stiffness matrix}
[sdSln, K, M] = SBFEMAssemblage(coord, sdConn, sdSC, D, den);

%ltx unit y-displacement at Node 5
nd = size(K,1);
U = zeros(nd,1); U(10) = 1;  
%ltx integration constants for specified displacement
sdPstP = SubdomainPostP2NodeEle( U, sdSln );

isd = 1; %ltx subdomain number
nNode = size(coord,1);
nodes = [4 5 6 1 2 3 ]'; %ltx nodes in sequence to form a open loop
xi = (1:-0.02:0).^2; %ltx radial coordinate
X = zeros(length(xi), nNode); Y = X; Z = X; %ltx initialising variables 
%ltx displacements and strains at the specified raidal coordinate
for ii= 1:length(xi)
    %ltx displacements at specified $\xi$
    [nodexy, dsp, strnNode, GPxy, strnEle] = ...
        SubdomainInDispStrain(xi(ii), sdSln{isd}, sdPstP{isd});
    X(ii,:) =  (nodexy(nodes,1)'); %ltx $x$-coordinates
    Y(ii,:) =  (nodexy(nodes,2)'); %ltx $y$-coordinates
    Z(ii,:) =  (dsp(2*nodes)');  %ltx $y$-displacement
end

%ltx
figure('Color','white');
surf(X,Y,Z,'FaceColor','interp', 'EdgeColor','none',  ...
           'FaceLighting','phong');
daspect([1 1 2]);
view(40, 30);
hold on
text(1.1*(coord(nodes,1)-0.02), 1.1*coord(nodes,2), Z(1,1:end)'+0.05,int2str(nodes));
plot3(X(1,:), Y(1,:), Z(1,:), '-k'); %ltx edges
axis equal; axis on;
xlabel('x'); ylabel('y'); zlabel('N');
title('SHAPE FUNCTION')
%ltx gray scale shading
invGray = 0.9*(1-gray); 
colormap(invGray)