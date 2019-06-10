close all; clearvars; dbstop error;

p = [    -0.0000      0.4854;    -0.0000      1.5174;
   -0.0000      1.0308;     0.0000      0.0000;
    0.0000      2.0000;     0.4828      2.0000;
    0.4960      0.7059;     0.4980      0.0000;
    0.5994      1.3949;     0.9679      2.0000;
    1.0000      0.4734;     1.0000      0.0000;
    1.0000      1.0000;     1.2864      1.5042;
    1.5142      2.0000;     1.5267      1.0000;
    2.0000      1.0000;     2.0000      2.0000;
    2.0000      1.5021 ];
%ltx triangles
t = [ 19      14     16 ;    16      14     13 ;
    13       7     11 ;    19      18     15 ;
    15      14     19 ;    13      14      9 ;
     2       3      9 ;     9       7     13 ;
     9       3      7 ;     7       3      1 ;
    19      16     17 ;     8      11      7 ;
     8       1      4 ;     7       1      8 ;
    14      15     10 ;    10       9     14 ;
    11       8     12 ;     9      10      6 ;
     6       5      2 ;     2       9      6 ];

 % Plot mesh of triangular elements
opt =struct('LabelEle', 11, 'LabelNode', 10, 'LineWidth',1); % plotting options
figure ('Position', [100 100 340 340])
PlotTriFEMesh( p, t, opt)
axis off

 [ coord,  sdConn, sdSC ] = triToSBFEMeshWithPlot( p, t );
 
figure ('Position', [100 100 340 340])
  opt=struct('LineSpec','-k', 'sdSC', sdSC, 'LabelSC', 10, ...
      'fill', [0.8 0.8 0.8], 'PlotNode', 1, 'LabelNode', 9, ...
      'MarkerSize', 6);
PlotSBFEMesh(coord, sdConn, opt); 
axis off;
