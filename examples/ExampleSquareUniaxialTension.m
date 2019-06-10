%% Square plate (The horizontal displacement of left edge and the vertical
%  displacement of the middle of the left edge are fixed. 
%  A uniform tension p=1 is applied on the right edge)

%%
%clear all
clearvars;
close all
dbstop if error

%% x coordinates and y coordinates of the node are defined in appropriate
% order. Example :- Node 1 is located at (0,0), node 2 at (0.5,0) and so on
coord =[   0         0;
    0.5000         0;
    1.0000         0;
    0    0.5000;
    0.5000    0.5000;
    1.0000    0.5000;
    0    1.0000;
    0.5000    1.0000;
    1.0000    1.0000;
    0.2500         0;
    0    0.2500;
    0.2500    0.2500;
    0.5000    0.2500;
    0.2500    0.5000;
    0    0.7500;
    0.2500    0.7500;
    0.5000    0.7500;
    0.2500    1.0000;
    0.1250         0;
    0    0.1250;
    0.1250    0.1250;
    0.2500    0.1250;
    0.1250    0.2500;
    0.1250    0.7500;
    0    0.8750;
    0.1250    0.8750;
    0.2500    0.8750;
    0.1250    1.0000];

% Nodes that constitute an element are defined here following anti-clockwise
% direction (including the hanging nodes)

ele = { [1    19    21    20];
    [ 2     3     6     5    13];
    [ 4    14    16    24    15];
    [ 5     6     9     8    17];
    [  10     2    13    12    22];
    [  11    23    12    14     4];
    [  12    13     5    14];
    [   14     5    17    16];
    [   15    24    26    25];
    [   16    17     8    18    27];
    [   19    10    22    21];
    [   20    21    23    11];
    [   21    22    12    23];
    [   24    16    27    26];
    [   25    26    28     7];
    [   26    27    18    28]};

% The support conditions are defined here. The first column indicates the
% node number. The second number represents the constraint. '1' is used to
% represent u_x = 0 and '2' represents u_y = 0
BC_Disp = [7 1 0;
    25 1 0;
    15 1 0;
    4 1 0;
    4 2 0;
    11 1 0;
    20 1 0;
    1  1 0];
%The applied loads are defined here. The first column represents the node
%number. The second represents the load in x-direction and the thirs
%defines load in y-direction.
BC_Frc = [9 1 0.25;
    6 1 0.5;
    3 1 0.25];

%The material constants for the isotropic material are defined here and the
%material matrix is constructed.
E = 10;
pos = 0;
D=E/(1-pos^2)*[1 pos 0;pos 1 0;0 0 (1-pos)/2];
den = 0;


nele = length(ele);
sdConn = cell(nele,1);
sdSC = zeros(nele,2);
for ie =1:length(ele)
    sdConn{ie}=[ele{ie}; ele{ie}(2:end) ele{ie}(1)]';
    sdSC(ie,:) = mean(coord(ele{ie},:));
end

%ltx solution of subdomains and global stiffness and mass matrices
[sdSln, K, M] = SBFEMAssemblage(coord, sdConn, sdSC, D, den);

%ltx {\bf Static solution}
[U, F] = SolverStatics(K, BC_Disp, BC_Frc);

% ltx plot deformed shape
figure
Umax = max(abs(U)); %ltx maximum displacement
fct = 0.2/Umax; %ltx factor to magnify the displacement to 0.2 m
%ltx augment nodal coordinates
deformed = coord + fct*(reshape(U,2,[]))'; 
hold on
%ltx plot undeformed mesh
opt = struct('LineSpec','-b'); %ltx plotting option
PlotSBFEMesh(coord, sdConn, opt);
%ltx plot deformed mesh
opt = struct('LineSpec','-k'); %ltx plotting option
PlotSBFEMesh(deformed, sdConn, opt);
title('DEFORMED MESH');

%% ------------------------------------------------------------------------

