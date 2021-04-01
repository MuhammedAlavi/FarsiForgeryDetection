classdef TargetHandWrit < HandWrit
    %HANDWRIT class for each handwritten image
    properties
        % BASE CLASS PROPERTIES
        %RGB     % [Z x P x 3 (double)] RGB image
        %Gray    % [Z x P (int)] gray level image
        %BW      % [Z x P (bool)] black and white image
        
        %DERIVED PROPERTIES
        %CC      % [N x 1 (ConnComp)] connected components
        %Skew    % [1 x 1 (double)] angle of skewness of image
        %DistanceMat %[NCC x NCC (int)] distance between CCs
        ClusterID   % [N x 1 (int)] number of cluster for each CC
        ClusterInfo %[N x 1 (struct)] clustering more informations        
    end
    properties(Dependent = true)
        % BASE DEPENDENT PARAMETERS
        %nCC     % [1 x 1 (int)] number of connected components
        %nTrace  % [1 x 1 (int)] number of connected components with Trace
        %ImgSize % [1 x 2 (int)] size of image
        %K       % [scalar (int)] count of clusters
        ClustOrd % [ vector (int)] order of cluster by it's area size
    end
    %INITIZLIZTION
    methods
        function this = TargetHandWrit(Img,SqrSz,MaxJunction,MaxEdge)
            %TARGETHANDWRIT constructor for a target hand-written
            this@HandWrit(Img,SqrSz,MaxJunction,MaxEdge);
        end
    end
    %GET FUNCTIONS
    methods
        function s = get.ClustOrd(this)
            C = this.ClusterInfo.C;
            C = C(:,1) .* C(:,2);
            [~,s] = sort(C,'ascend');
        end
        function CCID = getCCID(this,ClNum)
            %GETCCID return CC identifiers in cluster ClNum
            % Input:
            %       this: current object
            %       clNum: [scalar (int)] number of cluster
            % Output:
            %       CCID: [this.nCC x 1 (logical)] identifier of clusters
            
            CCID = this.ClusterID == ClNum;
        end
        function CCID = getAllID(this)
            %GETALLID return all CC IDs that is available in our aim CC set
            % Input:
            %       this: current object
            % Output:
            %       CCID: [this.nCC x 1 (logical)] id of all CC in aim CC
            %       set
            CCID = false(this.nCC,1);
            for i=1:this.nCC
               CCID(i) = this.CC(i).ST ~= StrokeType.Unknown; 
            end
        end
    end
    %COMPARISON
    methods
        function cluster(this,K,Verbose)
            %CLUSTER, clustering hand-written by it's stroke size and
            % speciry CCs that are point (both dot point and line points)
            
            if ~exist('Verbose','var')
                Verbose = false;
            end
            
            [G,IDOfCCs] = this.getThumbSize;
            
            % center of gravity of CC sizes
            M = mean(G);
            
            % divide all element sizes by their mean
            G(:,1) = G(:,1) ./ M(1);
            G(:,2) = G(:,2) ./ M(2);
            
            % clustering
            [Idx,MoreInf] = clusterHelper(G,K);
            
            % show clusters
            if Verbose
                figure,hold on;
                for i=1:numel(unique(Idx))
                    clusterIdx = Idx == i;
                    scatter(G(clusterIdx,1),G(clusterIdx,2),'fill');
                end
            end
            
            this.ClusterID = zeros(this.nCC,1);
            this.ClusterID(IDOfCCs) = Idx;
            this.ClusterInfo = MoreInf;
            
             % set value of clusters by their size
             S = this.ClustOrd;
             NewClust = zeros(size(this.ClusterID));
             for i=1:numel(S)
                 I = this.ClusterID == S(i);
                 NewClust(I) = i;
             end
             this.ClusterID = NewClust;
        end
        function CmpArray = search(this,AnotherHW,I,Z,W,F)
            %CMPARRAY compare each CC in cluster I and return number of Z 
            % Cmp object with minimum distance.
            % Input:
            %       this: current object
            %       AnotherHW: the other hand-written object
            %       I: [scalar (int)] number of cluster
            %       Z: [scalar (int)] number of most resemble CCs
            %       W: [vector (int)] list of window sizes
            %       F: [vector (double)] fraction of Target CC size with
            %       Test CC size
            % Output:
            %       CmpArray: [ M x 1 (cmp object)] comparison object
            
            thisID = this.getCCID(I);
            anotherID = AnotherHW.getCCID(I);
            
            % make all logical indexing to linear indexing
            thisID = find(thisID);
            anotherID = find(anotherID);
            
            N = numel(thisID);
            M = numel(anotherID);
            
            % there is and exception: if Train or Test have a cluster
            % without member, then comparison failed.
            if N == 0
                CmpArray = [];
                return
            elseif M == 0
                CmpArray = [];
                return
            end
            CmpArray = cell(N,1);
            
            % initialize CC array list 
            CCList = ConnComp;
            CCList(M) = ConnComp;
            
            for j=1:M
                CCList(j) = AnotherHW.CC(anotherID(j));
            end
            
            c = 1;
            for i=1:N
                % compare each CC in Target with a CC in CCList
                CmpResult = connCompSearch(this.CC(thisID(i)),CCList,Z,W,F);
                CmpArray{c} = CmpResult;
                c = c + 1;
            end
           CmpArray = cell2mat(CmpArray); 
        end
        function CmpArray = fullSearch(this,AnotherHW,Z,W)
            %CMPARRAY compare each CC in cluster I and return number of Z 
            % Cmp object with minimum distance.
            % Input:
            %       this: current object
            %       AnotherHW: the other hand-written object
            %       Z: [scalar (int)] number of most resemble CCs
            %       W: [vector (int)] list of window sizes
            %       F: [vector (double)] fraction of Target CC size with
            %       Test CC size
            % Output:
            %       CmpArray: [ M x 1 (cmp object)] comparison object
            
            thisID = this.getAllID();
            anotherID = AnotherHW.getAllID();
            
            % make all logical indexing to linear indexing
            thisID = find(thisID);
            anotherID = find(anotherID);
            
            N = numel(thisID);
            M = numel(anotherID);
            
            % there is and exception: if Train or Test have a cluster
            % without member, then comparison failed.
            if N == 0
                CmpArray = [];
                return
            elseif M == 0
                CmpArray = [];
                return
            end
            CmpArray = cell(N,1);
            
            % initialize CC array list 
            CCList = ConnComp;
            CCList(M) = ConnComp;
            
            for j=1:M
                CCList(j) = AnotherHW.CC(anotherID(j));
            end
            
            c = 1;
            for i=1:N
                % compare each CC in Target with a CC in CCList
                CmpResult = connCompSearch(this.CC(thisID(i)),CCList,Z,W);
                CmpArray{c} = CmpResult;
                c = c + 1;
            end
           CmpArray = cell2mat(CmpArray); 
        end
    end
    %DISPLAY
    methods
        function displayClusters(this,SpecificCluster)
            %DISPLAYCLUSTERS, show all cluster CCs
            % Input:
            %       SpecificCluster:[1 x 1 (int)] number of specific cluster
            if ~exist('SpecificCluster','var')
                UniqueClusterID = unique(this.ClusterID);
                UniqueClusterID = UniqueClusterID(:)';
                for i=UniqueClusterID
                    Idx = this.ClusterID == i;
                    figure
                    this.displayImg(Idx);
                   
                    % announce number of strokes
                    L = length(Idx);
                    ILen = sum(Idx);
                    A = @(String) num2str(String);
                    TitleText = strcat('Cluster Number:',A(i),...
                        ' Number of StrokesL',A(ILen),...
                        '(',num2str(round((ILen/L)*100)),'%) of ',A(L));
                    title(TitleText);
                end
            else
                ClusterNumbers = unique(this.ClusterID);
                if ~ismember(SpecificCluster,ClusterNumbers) && SpecificCluster >= 0
                    error('Not such cluster ID!')
                end
                Idx = this.ClusterID == SpecificCluster;
                figure
                this.displayImg(Idx);
                
                % announce information
                L = length(Idx);
                ILen = sum(Idx);
                s = SpecificCluster;
                A = @(String) num2str(String);
                TitleText = strcat('Cluster Number:',A(s),...
                    ' Number of StrokesL',A(ILen),...
                    '(',num2str(round((ILen/L)*100)),'%) of ',A(L));
                title(TitleText);
            end
        end
    end 
    
end

function [ClustIdx,MoreInfo] = clusterHelper(Data,K)
%CLUSTERHELPER is a separate function for clustering
% Input:
%       Data: [N x M ] Data points
%       K : [scalar (int)] number of clusterss
% Output:
%       ClustIdx: [N x 1 (int)] cluster index
%       Moreinfo: [N (structure)] more information about clustering

[ClustIdx,C,sumd,D] = kmedoids(Data,K,'Distance','mahalanobis','Replicates',10);

% more information about clustering
MoreInfo.C = C;
MoreInfo.sumd = sumd;
MoreInfo.D = D;
end
