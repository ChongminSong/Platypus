%ltx {\bf Shape function of a pentagon S-Element}
clearvars; close all;

%ltx {\bf Mesh}
%ltx nodal coordinates
coord = [cosd(-126:72:180);sind(-126:72:180)]';
%ltx connectivity
sdConn = { [1:5 ; 2:5 1]' }; 
sdSC = [ 0 0]; %ltx one S-Element per row

%ltx {\bf Materials}
%ltx \texttt{E}: Young's modulus; \texttt{p}: Poisson's ratio
ElasMtrx = @(E, p) E/(1-p^2)*[1 p 0;p 1 0;0 0 (1-p)/2];
mat.D = ElasMtrx(10E6, 0.25); %ltx  \texttt{E} in KPa
mat.den = 2; %ltx mass density in $\unit{Mg/m^3}$

%ltx {\bf S-Element solution}
[sdSln, K, M] = SBFEMAssembly(coord, sdConn, sdSC, mat);

%ltx {\bf Shape function}
%ltx prescribed unit displacement
U = zeros(size(K,1),1);
U(2*4) = 1;

%ltx strain modes of S-Element
sdStrnMode =  SElementStrainMode2NodeEle( sdSln );

%ltx integration constants 
sdIntgConst = SElementIntgConst( U, sdSln );

isd = 1; %ltx S-Element number
xi = 1:-0.01:0; %ltx radial coordinates 
%ltx initialisation of variables for plotting
X = zeros(length(xi), length(sdSln{isd}.node)+1);
Y = X; Z = X;
%ltx displacements and strains at the specified raidal coordinate
for ii= 1:length(xi)
    [nodexy, dsp, strnNode, GPxy, strnEle] = ...
        SElementInDispStrain(xi(ii), sdSln{isd},  ...
                              sdStrnMode{isd}, sdIntgConst{isd});
    %ltx coordinates of grid points forming a close loop
    X(ii,:) =  [nodexy(:,1)' nodexy(1,1)];
    Y(ii,:) =  [nodexy(:,2)' nodexy(1,2)];
    Z(ii,:) =  [dsp(2:2:end)' dsp(2)]; %ltx store $u_{y}$ for plotting
end

%ltx plot the shape function as a surface
figure('Color','white')
surf(X,Y,Z,'FaceColor','interp', 'EdgeColor','none',  ...
           'FaceLighting','phong');
view(-110, 15); %ltx set direction of viewing
hold on
text(1.1*(coord(:,1)-0.02), 1.1*coord(:,2), ...
     Z(1,1:end-1)'+0.05,int2str((1:5)')); %ltx label the nodes
axis equal, axis off; 
xlabel('x'); ylabel('y'); zlabel('N'); %ltx label the axes
colormap(jet)
% % invGray = 0.9*(1-gray); 
% % %ltx gray scale plot for inclusion in book. Remove for colour plotting
% % colormap(invGray)
plot3(X(1,:), Y(1,:), Z(1,:), '-b'); %ltx plot edges 

%ltx contour of the shape function
h = figure('Color','white');
contourf(X,Y,Z,10, 'LineStyle','none'); %ltx 10 contour lines
hold on
text(1.05*(coord(:,1)-0.02), 1.05*coord(:,2),int2str((1:5)'));
axis equal; axis off;
%ltx show a colorbar indicating the value of the shape function
colormap(jet)
caxis([0 1]); colorbar;
% % %ltx gray scale plot for inclusion in book. Remove for colour plotting
% % colormap(invGray); 
plot(X(1,:), Y(1,:), '-b'); %ltx plot edges
