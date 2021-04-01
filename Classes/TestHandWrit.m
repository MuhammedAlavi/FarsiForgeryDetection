classdef TestHandWrit < HandWrit
    %HANDWRIT class for each handwritten image
    properties
        %BASE PARAMETERS
        %RGB     % [Z x P x 3 (double)] RGB image
        %Gray    % [Z x P (int)] gray level image
        %BW      % [Z x P (bool)] black and white image
        
        %DERIVED PROPERTIES
        %CC      % [N x 1 (ConnComp)] connected components
        %Skew    % [1 x 1 (double)] angle of skewness of image
        %DistanceMat %[NCC x NCC (int)] distance between CCs
        ClusterID   % [N x 1 (cell)] number of cluster for each CC
        ClusterInfo %[N x 1 (struct)] clustering more informations
    end
    properties(Dependent = true)
        % BASE DEPENDENT PARAMETERS
        %nCC     % [1 x 1 (int)] number of connected components
        %nTrace  % [1 x 1 (int)] number of connected components with Trace
        %ImgSize % [1 x 2 (int)] size of image
        NClust   % [scalar (int)] last cluster (Excess CCs)
    end
    %INITIZLIZTION
    methods
        function this = TestHandWrit(Img,SqrSz,MaxJunction,MaxEdge)
            %TESTHANDWRIT constructor for test hand-writtens
            this@HandWrit(Img,SqrSz,MaxJunction,MaxEdge);
        end
    end
    %GET FUNCTION
    methods
        function k = get.NClust(this)
            k = max(this.ClusterID(:));
        end
        function CCID = getCCID(this,ClNum)
            %GETCCID return CC identifiers in cluster ClNum
            % Input:
            %       this: current object
            %       clNum: [scalar (int)] number of cluster
            % Output:
            %       CCID: [this.nCC x 1 (logical)] identifier of clusters
            
            CCID = this.ClusterID(:,ClNum);
            CCID = CCID > 0;
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
        function clusterByCircle(this,M,O,R,Epsilon,Verbose)
            %CLUSTERBYCIRCLE cluster data with circles.
            % Input:
            %       this: current object
            %       M: [scalar (double)] center of gravity in Target CC
            %       sizes
            %       O:[N x 2 (double)] center of circles
            %       R:[vector (double)] radius of circles
            % NOTE: length of O is number of clusters(suppose L). 
            % and if any of data didn't match inside of any circle,
            % it will go to (L + 1)Th cluster. 
            
            % YOU CAN ADD A PLAN TO CALCULATE EPSILON 
            if ~exist('Epsilon','var')
                Epsilon = 0;
            end
            
            if ~exist('Verbose','var')
                Verbose = false;
            end
            
            NCl = size(O,1); % number of clusters
            
            [G,IDOfCCs] = this.getThumbSize();
            
            % center of gravity in current object CC size 
            ThisM = mean(G);
            
            % move center of gravity of this CC size to Target CC size
            V = M - ThisM;
            G(:,1) = G(:,1) + V(1);
            G(:,2) = G(:,2) + V(2);
            G = round(G);
            
            Viz = isnan(IDOfCCs);% Trick:logical matrix with same size
            
            % visualize the result
            if Verbose
                figure,hold on;
                scatter(G(:,1),G(:,2),'fill');
%                 scatter(O(:,1),O(:,2),'k','d',12);
                for i=1:size(O,1)
                    center = O(i,:);
                   viscircles(center,R(i));
                   scatter(O(i,1),O(i,2),'fill','black');
                end
            end
            
            % initialize ID Of Clusters
            this.ClusterID = zeros(this.nCC,NCl);
            
            for i=1:NCl
                Idx = circularSampling(O(i,:),R(i),G,Epsilon(i));
                
                this.ClusterID(IDOfCCs(Idx),i) = i;
                Viz(Idx) = true;
            end
            
            % if there is a CC without cluster
            if ~all(Viz)
                this.ClusterID(IDOfCCs(~Viz),NCl+1) = NCl + 1;
            end
            
        end
    end
    %DISPLAY
    methods
        function displayClusters(this,SpecificCluster)
            %DISPLAYCLUSTERS, show all cluster CCs
            % Input:
            %       SpecificCluster:[1 x 1 (int)] number of specific cluster
            if ~exist('SpecificCluster','var')
                for i=1:this.NClust
                    Idx = this.ClusterID(:,i);
                    Idx = Idx > 0;
                    
                    % if there is not such cluster in cluster ID
                    if all(~Idx)
                        warning('Not such cluster ID!')
                        return
                    end
                    
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
                    warning('Not such cluster ID!')
                    return
                end
                Idx = this.ClusterID(:,SpecificCluster);
                Idx = Idx>0;
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
    %PRIVATE
    methods(Access = private)
    end
end
