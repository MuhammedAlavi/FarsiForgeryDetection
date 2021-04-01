classdef HandWrit < handle
    %HANDWRIT class for each handwritten image
    properties
        %BASE PROPERTIES
        RGB     % [Z x P x 3 (double)] RGB image
        Gray    % [Z x P (int)] gray level image
        BW      % [Z x P (bool)] black and white image
        K       % [scalar (int)] cluster number
        
        %DERIVED PROPERTIES
        CC      % [N x 1 (ConnComp)] connected components
        Skew    % [1 x 1 (double)] angle of skewness of image
    end
    properties(Dependent = true)
        nCC     % [1 x 1 (int)] number of connected components
        nTrace  % [1 x 1 (int)] number of connected components with Trace
        ImgSize % [1 x 2 (int)] size of image
    end
    %INITIZLIZTION
    methods
        function this = HandWrit(Img,SqrSz,MaxJunction,MaxEdge)
            %HANDWRIT is object of hand-written image
            % Input:
            %       RGB: RGB image of handwritten
            %       SqrSz: size of indicator square
            
            [RGB,Skew,Frame] = extractHandWrit(Img,SqrSz);
            Gray = rgb2gray(RGB);
            % enhance gray level image for better BW image extraction
            EnhancedGray = grayEnhance(Gray);
            Thresh = graythresh(EnhancedGray);
            BW = im2bw(EnhancedGray,Thresh);
            
            % complement both gray and black and white image
            BW = imcomplement(BW);
            Gray = imcomplement(Gray);
            
            this.Skew  = Skew;
            this.RGB = RGB;
            this.BW = BW;
            this.Gray = Gray;
            BWCC = bwconncomp(BW);
            %  REMOVE CCs ON BORDER !!!!!
            BWCC = removeOnBorder(BW,BWCC,Frame);
            
            this.CC = ConnComp;
            this.CC(BWCC.NumObjects) = ConnComp;
            for i=1:BWCC.NumObjects
                P = BWCC.PixelIdxList{i};
                Sz = size(Gray);
                this.CC(i) = ConnComp(P,Sz,i);
                this.CC(i).initConnComp(Skew,MaxJunction,MaxEdge);
                this.CC(i).initInkModel(this.Gray);
            end
        end
    end
    %GET FUNCTION
    methods
        
        function n = get.nCC(this)
            n = length(this.CC);
        end
        function n = get.nTrace(this)
            n = 0;
            for i=1:this.nCC
                if this.CC(i).ST ~= StrokeType.Unknown
                    n = n + 1;
                end
            end
        end
        function s = get.ImgSize(this)
            s = size(this.BW);
        end
        
        function [G,IDOfCCs] = getThumbSize(this)
            % return matrix of Thumb size (image of CC) and their
            % corresponding index.
            
            G = zeros(this.nCC,2);
            IDOfCCs = zeros(this.nCC,1);
            
            c = 1;
            for i=1:this.nCC
                
                % if stroke has no without trace, continue
                if this.CC(i).ST == StrokeType.Unknown
                    continue
                end
                
                %extrace size of CCs
                G(c,1:2) = (this.CC(i).ThumbSize(:))';
                IDOfCCs(c) = i;
                c = c + 1;
            end
            
            IDOfCCs = IDOfCCs(1:c-1,:);
            G = G(1:c-1,:);
            
        end
        function [S,I] = sortCCBySize(this,type)
            
            % validate type of sorting
            type = lower(type);
            type = validatestring(type,{'ascend','descend'});
            
            % capture size of CCs
            ThumbSize = this.getData('ThumbSize');
            
            % capture area of CC rectangle
            CCArea = ThumbSize(:,1).*ThumbSize(:,2);
            
            % sorted Array and their Index
            [S,I] = sort(CCArea,type);
            
            % remove zero elements in S and I (Those are [0,0], are
            % unknown CCs)
            if strcmp(type,'ascend')
                FirstIndex = find(S,1,'first');
                S = S(FirstIndex:end);
                I = I(FirstIndex:end);
            else
                LastIndex = find(S,1,'last');
                S = S(1:LastIndex);
                I = I(1:LastIndex);
            end
        end
        function Index = getCCWithType(this,StrkType)
            %GETCCWITHTYPE return index of all CCs with specifi Stroke
            %types
            % Input:
            %       this: current object
            %       StrkType:[enum of class StrokeType]
            % Output:
            %       Index: [this.nCC x 1 (logical)] index of CCs with
            %       specific stroke type
            
            Index = false(this.nCC,1);
            for i=1:this.nCC
                for j=1:length(StrkType)
                    if this.CC(i).ST == StrkType(j)
                        Index(i) = true;
                    end
                end
            end
        end
        function Index = getCCWithClusterNum(this,CNum)
            %GETCCWITHCLUSTERNUM return CCs with a spacific cluster number.
            % Input:
            %       this: current object
            %       CNum: [scalar (int)] Number of cluster
            ClusterElements = this.ClusterID == CNum;
            Index = find(ClusterElements);
        end
    end
    %DISPLAY
    methods
        function displayImg(this,Index)
            if ~exist('Index','var')
                Index = 1:this.nCC;
            end
            if islogical(Index)
                Index = find(Index);
            end
            
            % be sure that Index is a row vector
            Index = Index(:);
            Index = Index';
            
            % create a blank page
            Img = false(this.ImgSize);
            
            for i=Index
                Img(this.CC(i).PixelIdx) = true;
            end
            imshow(Img);
            
        end
        function displayGraph(this,Index)
            if ~exist('Index','var')
                Index = 1:this.nCC;
            end
            if islogical(Index)
                Index = find(CCIndex);
            end
            
            % be sure that Index is a row vector
            Index = Index(:);
            Index = Index';
            
            Img = zeros(this.ImgSize);
            for c = Index
                Img(this.CC(c).PixelIdx) = 1;
            end
        end
    end
    %PRIVATE
    methods(Access = private)
        function [NewClusters,K1] = DataWithLinearityFeature...
                (this,Idx,IDOfCCs,K1,K2,CurrentCluster)
            %DATAWITHLINEARITYFEATURE ignore dot cluster and add a new
            % feature to Data.
            % Input:
            %       this: current object
            %       K: [scalar (int)] number of cluster
            %       Idx: [N x 1 (int)] index of clusters
            %       ID: [N x 1 (int)] corresponding CC Id of cluster elements
            %       C: [K x 2 (double)] center point of each cluster
            % Output:
            %       G: [ N x 3 (int)] size and linearity feature
            %       ID:[ N x 1 (int)] id of CCs
            %       K :[scalar (int)] cluster numbers
            %
            % Note: linearity feature is Square error in Least Square Error
            % algorithm.
            
            % CHECK EXCEPTIONS:
            % 1. if number of elements in a cluster is less than sub
            % cluster number
            if sum(Idx) < K2
                NewClusters = Idx * CurrentCluster;
                return
            end
            
            
            % gather all the data
            ClusterIndex = find(Idx);
            
            % initialize data vector
            G = zeros(length(ClusterIndex),2);
            
            c = 1;
            ClusterIndex = ClusterIndex(:)';
            for i=ClusterIndex
                % linearity feature
                % fit a line to trace data
                X = this.CC(IDOfCCs(i)).Trace(:,2);
                Y = this.CC(IDOfCCs(i)).Trace(:,1);
                
                % fit a line into data points
                p = polyfit(X,Y,1);
                Yk = polyval(p,X);
                % because it is image and discrete so we make it round
                Yk = round(Yk);
                % Square Error
                LSQ =   sum(  ( Y - Yk ).^2  )  ;
                G(c,1) = LSQ/length(X);
                G(c,2) = p(1);
                c = c + 1;
            end
            % %             boxplot(G(:,1));figure
            % %             boxplot(G(:,2));figure
            % %             this.displayImg(IDOfCCs(ClusterIndex));
            % %             close all;
            % now divide each cluster into K sub cluster
            SubIdx = kmeans(G,K2);
            NewClusters = Idx * 0;
            ClusterIndex = ClusterIndex(:);
            
            Temp = SubIdx == 1;
            Temp = ClusterIndex(Temp);
            NewClusters(Temp) = CurrentCluster;
            
            % new cluster called with a new floating point number
            % for example 1
            for i=2:K2
                Temp = SubIdx == i;
                Temp = ClusterIndex(Temp);
                NewClusters(Temp) = K1+1;
                K1 = K1 + 1;
            end
        end
        function subClusterDistance(this,ExceptionID)
            %SUBCLUSTERDISTANCE is distance between each CC inside of each
            % cluster.
            % Input:
            %        ExceptionID: [ vector (int)] list of point cluster that
            %        should not take part in comparison
            
            for c=1:this.K
                
                % if exception happened
                if ismember(c,ExceptionID)
                    continue
                end
                
                % for each cluster
                ClElements = (this.ClusterID == c);
                NElements = sum(ClElements);
                ClElements = find(ClElements);
                
                % if there is just one CC in a cluster
                if NElements == 1
                    % there is no distance
                    return
                end
                
                % for each cluster, measure distance of CCs one by one
                for i=1:NElements
                    Ang1 = this.CC(ClElements(i)).Ang;
                    for j=i+1:NElements
                        Ang2 = this.CC(ClElements(j)).Ang;
                        this.DistanceMat(ClElements(i),ClElements(j)) = dtw(Ang1,Ang2);
                        this.DistanceMat(ClElements(j),ClElements(i)) = ...
                            this.DistanceMat(ClElements(i),ClElements(j));
                    end
                end
                
            end
        end
    end
end
