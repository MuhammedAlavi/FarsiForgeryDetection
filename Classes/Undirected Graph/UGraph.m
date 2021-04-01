classdef UGraph < handle
    %
    % UGraph : Undirected graph for an image skeleton
    %    
    properties
        G  % [n x 2] node coordinates
        E  % [n x n boolean] adjacency matrix
        EI % [n x n cell], each cell is that edge's index into trajectory list
           %  S. It should be a cell array, because there could be two paths
        S  % [k x 1 cell] edge paths in the image
        I  % [nimg x nimg boolean] original image
    end
    properties (Dependent = true)
        % number of nodes in a graph
        NNode 
        % number of edges in a graph
        NEdge
        % image size
        imsize 
        % [k x 2] for each edge (row), the source and dest. node,
        % given the direction it is listed in S
        link_ei_to_ni    
        Degree% degree of nodes
    end
    properties (Constant = true)
       % circular masks for merging critical points
       ListMask = loadListMask();
    end
    methods
        function this = UGraph(G,E,EI,S,I)
            % constructor for making an undirected graph
           if ~exist('G','var') || isempty(G)
               this.G = [];
           else
               this.G = G;
           end
           if ~exist('E','var') || isempty(E)
               this.E = [];
           else
               this.E = E;
           end
           if ~exist('EI','var') || isempty(EI)
               this.EI = [];
           else
               this.EI = EI;
           end
           if ~exist('S','var') || isempty(S)
               this.S = [];
           else
               this.S = S;
           end
           if ~exist('I','var') || isempty(I)
               this.I = [];
           else
               this.I = I;
           end
        end
        function y = get.NNode(this)
           y = size(this.G,1);            
        end
        function e = get.NEdge(this)
            e = length(this.S);
        end
        function d = get.Degree(this)
            d = sum(this.E,2);
        end
        function sz = get.imsize(this)
            sz = size(this.I);            
        end
        function Y = get.link_ei_to_ni(this)
            % get what nodes the edges are connceted to
           k = numel(this.S);
           Y = zeros(k,2);
           for i=1:k
               source_pt = this.S{i}(1,:);
               dest_pt = this.S{i}(end,:);
               source_ni = this.map_pts_to_ni(source_pt);
               dest_ni = this.map_pts_to_ni(dest_pt);
               Y(i,1) = source_ni;
               Y(i,2) = dest_ni;
           end
        end
        function plot_skel(this)
            % Plot the skeleton on-top of the image
           this.plot_on_image(~this.I);
        end
        function plot_circle(this)
            % Plot the maximum circle analysis, with the
            % original image and the skeleton
            I_cluster = this.max_circle_mask;
            I_mask = double(this.I);
            I_mask(I_cluster) = 0.5;
            this.plot_on_image(1-I_mask);
        end
        function clean_skeleton(this)
            % Merge groups of graph nodes connected
            % by the maximum circle criterion
           I_mask = this.max_circle_mask();
           
           % cluster nodes
           C = this.cluster_G(I_mask);
           nc = numel(unique(C));
           
           % until we no longer have things to remove
           while has_rmv(C)
              
              % for each cluster index, look for 
              % things to merge
              for i=1:nc
                  sel = C == i;
                  if sum(sel)>1
                     this.merge(sel,I_mask);
                     
                     % cluster nodes
                     C = this.cluster_G(I_mask);
                     [~,nc] = has_rmv(C);
                     break
                  end                  
              end
               
           end
           
           % do we still have groups to remove?
           %
           % max_c: maximum cluster index
           function [yyyy,max_c] = has_rmv(myC)
               myC(myC==0) = [];
               yyyy = numel(myC) ~= numel(unique(myC));
               max_c = max(myC);
           end           
           
        end
        function I_mask = max_circle_mask(this)
            % Image based on maximum circles surrounding
            % each critical point
            I_mask = false(this.imsize);
            for indx=1:this.NNode
                I_mask = this.max_circle(indx,I_mask);
            end
        end
        function merge(this,sel,I_mask)
            % merge a collection of nodes
            %
            % Input
            %   sel: [n x 1 logical] or [k x 1 index]
            %
            % process input
            if islogical(sel)
                assert(numel(sel)==this.NNode);
            else
                choose = false(this.NNode,1);
                choose(sel) = true;
                sel = choose;
            end
            
            % partition
            E_sel_to_not = this.E(sel,~sel);
            EI_sel_to_sel = this.EI(sel,sel);
            EI_sel_to_not = this.EI(sel,~sel);
            
            % create the new row in adjacency matrix, which 
            % includes all the nodes connected to the selected set/pair
            new_row_E = any(E_sel_to_not,1);
            new_row_EI = cell(size(new_row_E));
            for i=1:size(EI_sel_to_not,2);
                new_row_EI{i} = tovec(EI_sel_to_not(:,i));
            end
            list_ei_sel_to_sel = tovec(EI_sel_to_sel);
            list_ei_sel_to_sel = unique(list_ei_sel_to_sel);
            
            % see which of the edges, from sel to sel, should be removed.
            % we don't want to remove edges that make a long arc and then return
            nss = length(list_ei_sel_to_sel);
            rmv_edge = false(nss,1);
            for i=1:nss
               indx = list_ei_sel_to_sel(i);       
               
               % remove edge if all the pixels are in the circle criterion mask
               pix = this.S{indx};
               sz = size(I_mask);
               ind_px = sub2ind(sz,pix(:,1),pix(:,2));
               if all(I_mask(ind_px));
                    rmv_edge(i) = true; 
               end
               
            end
            keep_ei_sel_to_sel =  list_ei_sel_to_sel(~rmv_edge);
            rmv_ei_sel_to_sel = list_ei_sel_to_sel(rmv_edge);

            % remove node positions
            G_sel = this.G(sel,:);
            new_node = mean(G_sel,1);
            Grmv = this.G(sel,:);
            this.G(sel,:) = [];
            this.G = [this.G; new_node];

            % modify adjacency matrix with the new row
            this.E(sel,:) = [];
            this.E(:,sel) = [];
            this.E(end+1,:) = new_row_E;
            this.E(:,end+1) = [new_row_E,false];
            if ~isempty(keep_ei_sel_to_sel)
                this.E(end,end) = true;
            end
            
            % modify the list of edges with each cell
            this.EI(sel,:) = [];
            this.EI(:,sel) = [];
            this.EI(end+1,:) = new_row_EI;
            this.EI(:,end+1) = [new_row_EI, {keep_ei_sel_to_sel} ];  

            % modify the paths
            this.S(rmv_ei_sel_to_sel) = [];
            for i=1:numel(this.EI)
               mycell = this.EI{i};

               % update counts in adjacency matrix 
               for j=1:length(mycell)
                  el = mycell(j);
                  this.EI{i}(j) = el - sum(rmv_ei_sel_to_sel < el);
               end
            end

            % update the paths, stored in this.S, such that all 
            % nodes we have now replaced are updated with their
            % new coordinates
            for i=1:length(this.S)
                for j=1:size(Grmv,1);
                    node_rmv = Grmv(j,:);
                    d = pdist2(node_rmv,this.S{i});
                    swap = d < .001;
                    nswap = sum(swap);
                    this.S{i}(swap,:) = repmat(new_node,[nswap 1]);
                end
            end
            
            % Make sure all the trajectories are roughly uniform distance
% %             for i=1:length(this.S)
% %                 this.S{i} = expand_unif_interp(this.S{i},1);
% %             end
            
            this.assert_valid_graph();
        end
        function assert_valid_graph(this)
            % Check to make sure this is a valid graph structure
            % Throw assert if it is not
        
            
            n = this.NNode;
    
            % check that we have the correct sizes 
            % of adjacency matrices
            assert(size(this.G,1)==n);
            assert(size(this.E,1)==n);
            assert(size(this.EI,1)==n);

            % make sure all of the paths in S
            % have a corresponding edge
            v = tovec(this.EI);
            uq = unique(v);
            assert(isequal(uq(:),(1:length(uq))'));    

            % check on the paths throught the graph
            for i=1:size(this.EI,1) % for each edge in the graph
                for j=1:size(this.EI,2)

                    v = this.EI{i,j}; % get the associated path
                    if isempty(v)
                       assert(this.E(i,j)==false);
                    else
                       assert(this.E(i,j)==true);
                    end

                    for k=1:length(v)
                       indx = v(k);

                       % check that the start and end points
                       % align with the coordinates stored in G
                       mystart = this.S{indx}(1,:);
                       myend = this.S{indx}(end,:);               
                       assert(aeq(mystart,this.G(i,:)) || aeq(mystart,this.G(j,:)));
                       assert(aeq(myend,this.G(i,:)) || aeq(myend,this.G(j,:)));
                    end

                end
            end

            % check that the matrix E and EI are in correspondence
            for i=1:size(this.EI,1) % for each edge in the graph
                for j=1:size(this.EI,2)
                    if this.E(i,j)
                       assert(~isempty(this.EI(i,j))); 
                    end
                end
            end
                   
        end
        function vni = map_pts_to_ni(this,vpts)
            %
            % Map list of nodes (as points) to their indices
            %
            % Input
            %  vpts: [n x 2] points
            %
            % Output
            %  vni: [n x 1] indices
            
            [n,dim] = size(vpts);
            assert(dim==2);
            vni = zeros(n,1);
            for i=1:n
                pt = vpts(i,:);
                dist = pdist2(pt,this.G);
                assert(isvector(dist));
                dist = dist(:);
                [minval,minindx] = min(dist);
                assert(aeq(minval,0));
                vni(i) = minindx;
            end   
        end
        function [vei,vei_flip] = get_branches(this,pt)
            % Given the current point, where can we go?
            %
            % Input: pt: [1x2] current point we want to analyze
            %
            % Output
            %  vei: [m x 1] list of edges we can go to
            %  vei_flip: [m x 1 boolean] do we need to flip that edge
            %    to align it to the current point
            assert(numel(pt)==2);
            ni = this.map_pts_to_ni(pt);
            list_ei = this.EI(ni,:);
            vei = concat_cell(list_ei);
            vei = vei(:);
            nn = numel(vei);
            list_traj = this.S(vei);
            vei_flip = false(size(vei));
            for i=1:nn
               traj = list_traj{i};
               if isequal(pt,traj(1,:))
                   vei_flip(i) = false;
               elseif isequal(pt,traj(end,:))
                   vei_flip(i) = true;
               else
                   error('cannot map point to its paths'); 
               end
            end            
        end              
        
        % TRACE FUNCTIONS ------------------------
        function [Trace,StkType] = makeTrace(this,N)
            % TRACE ALGORITHM:
            % if one_node || more_than_N_nodes
            % Trace = -1;
            % StkType = 'Unknown';
            % elseif one_edge
            % this.simpleTrace();
            % elseif less_than_N_edge && more_than_1_edge
            % complexTrace()
            % end
            if this.NNode == 1 || this.NEdge > N
                Trace = -1;
                StkType = StrokeType.Unknown;
            elseif this.NEdge == 1
                Trace = this.simpleTrace;
                StkType = StrokeType.Simple;
            elseif this.NEdge > 1 && this.NEdge <= N
                Trace = this.complexTrace;
                if Trace < 0
                    StkType = StrokeType.Unknown;
                else
                    StkType = StrokeType.Complex;
                end
            end
        end
    end
    methods (Access = private)
        function Trace = simpleTrace(this)
            % if there is just one edge in stroke, trace is just traverse
            % of single edge from Right to Left
            RM = this.rightMostPoint;
            if all(this.S{1}(1,:) == RM)
                Trace = this.S{1};
            else
                Trace = this.S{1}(end:-1:1,:);
            end
        end
        function Trc = complexTrace(this)
            % if thre is more than one edge in out graph, trace will be
            % according to this algorithm:
            % #1: find right-most point with Degree 1 and traverse edge
            % #2: go to next single edge with shortest length
            RM = this.rightMostPoint;
            Edges = 1:this.NEdge;
            Trc = [];
            l = this.link_ei_to_ni;
            [Status,T,NextNode,Edges] = this.Trace(RM,Edges,l);
            if Status < 0
                Trc = -1;
                return
            end
            Trc = cat(1,Trc,T);
            c = 1;% a counter to avoid infinite loop
            while ~isempty(Edges) && c < 10000
                [Status,T,NextNode,Edges] = this.Trace(NextNode,Edges,l);
                if Status < 0
                    Trc = -1;
                    return
                end
                Trc = cat(1,Trc,T);
                c = c+1;
            end
            
            % there is something with algorithm
            if c >= 10000
                Trc = -1;
                return
            end
        end
        function RightMost = rightMostPoint(this)
            %RIGHMOSTPOINT for starting Trace.
            [~,Idx] = sort(this.G(:,2),'descend');
            for i=1:length(Idx)
                if this.Degree(Idx(i)) == 1
                    RightMost = this.G(Idx(i),:);
                    return
                end
            end
        end
        function [Status,Trace,NextNode,Edges] = Trace(this,Node,Edges,L)
            MEdge = this.middleEdges(Node,L);
            MEdge = intersect(MEdge,Edges);
            if numel(MEdge) > 1
                Status = -1;
                Trace = -1;
                NextNode = -1;
                Edges = -1;
                return
            end
            SEdge = this.singleEdges(Node,L);
            SEdge = intersect(SEdge,Edges,'stable');
            SEdge = SEdge(:)';% row vector
            Trace = [];
            
            % determine last edge to go to next node
            if numel(MEdge) > 0
                LastEdge = MEdge;
            else
                LastEdge = SEdge(end);
                SEdge = SEdge(1:end-1);
            end
            for i=SEdge
                if isequal(this.S{i}(1,:),Node)
                    t = [this.S{i};this.S{i}(end:-1:1,:)];
                    Trace = cat(1,Trace,t);
                elseif isequal(this.S{i}(end,:),Node)
                    t = [this.S{i}(end:-1:1,:);this.S{i}];
                    Trace = cat(1,Trace,t);
                end
            end
            
            % concatenate middle edge to trace
            if isequal(this.S{LastEdge}(1,:),Node)
                t = this.S{LastEdge};
                Trace = cat(1,Trace,t);
            elseif isequal(this.S{LastEdge}(end,:),Node)
                t = this.S{LastEdge}(end:-1:1,:);
                Trace = cat(1,Trace,t);
            end
            
            % remove traced edges
            VisEdges = union(SEdge,LastEdge);
            Edges = setdiff(Edges,VisEdges);
            % next Node is last point of trace
            NextNode = Trace(end,:);
            Status = 1; % every thing is OK!
        end
        function Edges = singleEdges(this,N,L)
            % return index of all single edge(s) connected to node N
            % Input:
            %       N: cordination of node
            %       L: link of edges by nodes(this.link_ei_to_ni)
            
            % find node index
            N = and(N(1) == this.G(:,1),N(2) == this.G(:,2));
            N = find(N);
            
            % there is just one path from node N
            member = ismember(L,N);
            if sum(member(:)) == 1
                member = sum(member,2);
                Edges = find(member);
                return
            end
            Edges = zeros(this.NEdge,1);
            c = 1;
            for i=1:length(L)
                if L(i,1) == N || L(i,2) == N
                    if L(i,1) == N
                        otherPt = L(i,2);
                        member = ismember(L,otherPt);
                        if sum(member(:)) == 1
                            Edges(c) = i;
                            c = c+1;
                        end
                    elseif L(i,2) == N
                        otherPt = L(i,1);
                        member = ismember(L,otherPt);
                        if sum(member(:)) == 1
                            Edges(c) = i;
                            c = c+1;
                        end
                    end
                end
            end
            if c > 1
                Edges = Edges(1:c-1);
            else
                Edges = [];
            end
            
            % sort edges
            len = zeros(size(Edges));
            for i=1:length(Edges)
                len(i) = length(this.S{i});
            end
            [~,Index] = sort(len,'ascend');
            Edges = Edges(Index);
            
        end
        function Edges = middleEdges(this,N,L)
            % return index of all edges connected to N 
            % Input:
            %       N: cordination of node
            %       L: link of edges by nodes(this.link_ei_to_ni)
            
            % find node index
            N = and(N(1) == this.G(:,1),N(2) == this.G(:,2));
            N = find(N);
            
            member = ismember(L,N);
            if sum(member(:)) == 1
                Edges = [];
                return
            end
            
            Edges = zeros(this.NEdge,1);
            c = 1;
            for i=1:length(L)
                if L(i,1) == N || L(i,2) == N
                    if L(i,1) == N
                        otherPt = L(i,2);
                        member = ismember(L,otherPt);
                        if sum(member(:)) > 1
                            Edges(c) = i;
                            c = c+1;
                        end
                    elseif L(i,2) == N
                        otherPt = L(i,1);
                        member = ismember(L,otherPt);
                        if sum(member(:)) > 1
                            Edges(c) = i;
                            c = c+1;
                        end
                    end
                end
            end
            if c >= 1
                Edges = Edges(1:c-1);
            else
                Edges = [];
            end
            
            % sort edges
            len = zeros(size(Edges));
            for i=1:length(Edges)
                len(i) = length(this.S{i});
            end
            [~,Index] = sort(len);
            Edges = Edges(Index);
        end
        % --------------------------------------
        function I_cluster = max_circle(this,indx,I_cluster)
            % find the maximum circle in the image surrounding a point
            pt = this.G(indx,:);
            pt = round(pt);
            nmask = numel(this.ListMask);
            bool_fit = false(nmask,1);
            for i=1:nmask
               x_mask = this.ListMask{i};               
               x_mask(:,1) = x_mask(:,1) + pt(1);
               x_mask(:,2) = x_mask(:,2) + pt(2);               
               x_mask = this.cut_image_plane(x_mask);
               xlind = sub2ind(this.imsize,x_mask(:,1),x_mask(:,2));      
               bool_fit(i) = all(this.I(xlind));
               if bool_fit(i)
                   I_cluster(xlind) = true; 
               else
                   return
               end
            end
        end 
        function x = cut_image_plane(this,x)
            % remove pixels that are out of the image plane
            sz = this.imsize;
            xx = x(:,1);
            xy = x(:,2);
            rmvx = xx<=0 | xx>sz(1);
            rmvy = xy<=0 | xy>sz(2);
            rmv = rmvx | rmvy;
            x(rmv,:) = [];
        end
        function plot_on_image(this,I)
            %  Plot the graph skeleton ontop of an image
            sz = size(I);
            if size(I,3) == 1
                I = repmat(I,[1 1 3]);
            end            
            hold on
            ns = length(this.S);
            image([1 sz(1)],[1 sz(2)],I);
            for i=1:ns
                stk = this.S{i};
                %color = rand(3,1);
                plot_traj(stk,'g');
            end    
            for i=1:this.NNode
               plot(this.G(i,2),this.G(i,1),'r.','MarkerSize',7); 
            end
            set(gca,'YDir','reverse','XTick',[],'YTick',[]);
            xlim([1 sz(1)]);
            ylim([1 sz(2)]);
            
        end
        function C = cluster_G(this,I_mask)
            % Cluster the vertices in G
            % into clusters based on the maximum circle criterion
            L = bwlabel(I_mask,4);
            GG = round(this.G);
            lind_G = sub2ind(this.imsize,GG(:,1),GG(:,2));
            C = L(lind_G);            
        end
    end
end
function v = concat_cell(vcell)
% flatten a cell array into a regular array
% by vert cat each cell

    assert(iscell(vcell));
    v = [];
    for i=1:numel(vcell)
        v = [v; vcell{i}];
    end
end
function plot_traj(stk,color)
% plot a stroke trajectory in image space,
% where the x and y dimensions are reversed
    ystk = stk(:,2);
    stk(:,2) = stk(:,1);
    stk(:,1) = ystk;       
    plot(stk(:,1),stk(:,2),'Color',color,'LineWidth',1);
end