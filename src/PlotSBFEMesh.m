function PlotSBFEMesh(coord, sdConn, opt)
%Plot polygon S-element mesh
%
%Inputs:
%  coord(i,:)           - coordinates of node i
%  sdConn{isd}(ie,:)    - S-element conncetivity. The nodes of
%                         line element ie in S-element isd.
%  plot options:
%   opt=struct('LineSpec','-b', 'sdSC', sdSC, 'LabelSC', 14, ...
%       'fill', [0.9 0.9 0.9], 'PlotNode', 1, 'LabelNode', 14, ...
%       'MarkerSize', 6, 'BC_Disp',BC_Disp, 'BC_Frc',BC_Frc);
%  where the options are:
%    opt.sdSC     : scaling centers of S-elements. 
%    opt.LabelSC  : If specified, plot a marker at
%                     the scaling centre
%                   = 0, do not label S-element
%                   > 0, show S-element number 
%                     If > 2, it also specifies the font size.  
%    opt.fill     : =[r g b]. Fill an S-element with color.
%                   opt.sdSC has also to be given.
%    opt.LineSpec : LineSpec of 'plot' function in Matlab
%    opt.LineWidth: LineWidth of 'plot' function in Matlab
%    opt.PlotNode : = 0, do not plot node symbol; otherwise, plot
%    opt.LabelNode: = 0, do not label nodes;
%                   > 0, show nodal number. If opt.LabelNode > 5,
%                        it specifies the font size.
%                   < 0, draw a marker only
%    opt.MarkerSize: marker size of nodes and scaling centres
%    opt.BC_Disp   : if specified, plot a marker at a node
%                    with prescribed displacement(s)
%    opt.BC_Frc    : if specified, plot a marker at a node
%                    with applied force(s)

LineWidth = 1;
LineSPec = '-k';
% use specified LineSpec if present
if nargin > 2 
    if isfield(opt, 'LineSpec') && ~isempty(opt.LineSpec)
        LineSPec = opt.LineSpec;
    end
    if isfield(opt, 'LineWidth') && ~isempty(opt.LineWidth)
        LineWidth = opt.LineWidth;
    end
end

nsd = length(sdConn); %ltx number of S-elements

%ltx fill S-elements by treating scaling centre and an edge as a triangle
meshEdge = vertcat(sdConn{:});
if nargin > 2
    if isfield(opt, 'sdSC') && ~isempty(opt.sdSC)
        if isfield(opt, 'fill') && ~isempty(opt.fill)
            p = [coord; opt.sdSC]; %ltx points
            nNode = length(coord); %ltx number of nodes
            %ltx appending scaling centres
            nEdge = cellfun(@length, sdConn);
            %ltx initilisation of array of scaling centre
            cnt = zeros(sum(nEdge),1); 
            ib = 1; %ltx starting index
            for ii = 1:nsd;
                ie = ib-1 + nEdge(ii); %ltx ending index
                cnt(ib:ie) = nNode + ii; %ltx scaling centre
                ib = ie + 1;
            end
            t = [meshEdge cnt]; %ltx triangles
            patch('Faces',t,'Vertices',p, ...
                'FaceColor',opt.fill,'LineStyle','none');
        end
    end
end

%ltx plot mesh 
meshEdge = vertcat(sdConn{:});
meshEdge = unique(sort(meshEdge,2),'rows');
hold on
X = [coord(meshEdge(:,1),1)'; coord(meshEdge(:,2),1)'];
Y = [coord(meshEdge(:,1),2)'; coord(meshEdge(:,2),2)'];
plot(X,Y,LineSPec,'LineWidth',LineWidth);
xlabel('x'); ylabel('y'); %ltx label axes

%ltx apply plot options
if nargin > 2
    if isfield(opt, 'MarkerSize') 
        markersize = opt.MarkerSize;
    else 
        markersize = 5;
    end
    if isfield(opt, 'sdSC') && ~isempty(opt.sdSC)
        if isfield(opt, 'LabelSC') && ~isempty(opt.LabelSC) 
            %ltx plot scaling centre
            plot(opt.sdSC(:,1), opt.sdSC(:,2), 'r+', ...
             'MarkerSize', markersize, 'LineWidth', LineWidth);
            plot(opt.sdSC(:,1), opt.sdSC(:,2), 'ro', ...
             'MarkerSize', markersize, 'LineWidth', LineWidth);
            if opt.LabelSC > 1
                if opt.LabelSC > 5
                    fontsize = opt.LabelSC;
                else 
                    fontsize = 12;
                end
                text(opt.sdSC(:,1), opt.sdSC(:,2), ...
                    [blanks(nsd)' int2str((1:nsd)')], ...
                    'Color','r', 'FontSize',fontsize);
            end
        end
    end
    if isfield(opt, 'PlotNode') && opt.PlotNode
        %ltx showing nodes by plotting a circle
        plot(coord(:,1), coord(:,2), 'ko', ...
            'MarkerSize', markersize, 'LineWidth', LineWidth);
    end
    if isfield(opt, 'LabelNode') && ~isempty(opt.LabelNode) ...
            && opt.LabelNode
        nNode = length(coord);
        if opt.LabelNode > 1
            fontsize = opt.LabelNode;
        else
            fontsize = 12;
        end
        %ltx showing nodes by plotting a circle
        plot(coord(:,1), coord(:,2), 'ko', ...
            'MarkerSize', markersize, 'LineWidth', LineWidth);
        %ltx label nodes with nodal number
        text(coord(:,1),coord(:,2), ...
            [blanks(nNode)' int2str((1:nNode)')], ...
            'FontSize',fontsize);
    end
    if isfield(opt, 'BC_Disp') && ~isempty(opt.BC_Disp)
        %ltx show fixed DOFs by a marker at the nodes
        Node = opt.BC_Disp(:,1);
        plot(coord(Node,1),coord(Node,2),'b>', ...
            'MarkerSize',8, 'LineWidth', LineWidth);
    end
    if isfield(opt, 'BC_Frc') && ~isempty(opt.BC_Frc)
        %ltx show DOFs carrying external forces by a marker at the nodes
        Node = opt.BC_Frc(:,1);
        plot(coord(Node,1),coord(Node,2),'m^', ...
            'MarkerSize',8, 'LineWidth', LineWidth);
    end
end
axis equal;
end