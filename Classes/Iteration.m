classdef Iteration < handle
    %VERIFYITERATION verification iteratoin.
    properties
        Z  % count of most resemble stroke
        K
        W
        MaxEdge
        MaxJunction
        SqrSz
    end
    properties(Constant = true)
        G = loadGuideTable();
    end
    methods
        function this = Iteration(W,SqrSz,Z,MaxEdge,MaxJunction)
                this.Z = Z;
                this.W = W;
                this.MaxEdge = MaxEdge;
                this.MaxJunction = MaxJunction;
                this.SqrSz = SqrSz;
        end
        %COMPARISON
        function Cmp = compareWithID(this,TargetID,TestID,TestType)
            TargetHW = this.loadHandWrit(TargetID,'T');
            TestHW = this.loadHandWrit(TestID,TestType);
            Cmp = this.compare(TargetHW,TestHW);
        end
        function Clust = clusterWithID(this,ID,nCluster)
            % load specific hand-writen
            Target = this.loadHandWrit(ID,'T');
            % cluster target and test handwritten with K 
            Target.cluster(nCluster);
            Clust = Target.ClusterID(Target.ClusterID > 0);
        end
        function Cmp = fullCompareWithID(this,TargetID,TestID,TestType)
            TargetHW = this.loadHandWrit(TargetID,'T');
            TestHW = this.loadHandWrit(TestID,TestType);
            Cmp = this.fullCompare(TargetHW,TestHW);
        end
        function Cmp = compare(this,Target,Test)
            %SEARCH search through test hand-writtens and find strokes with
            % most resemblance to strokes in train.
            % Input:
            %       TargetImg:[M x N x 3] RGB Target Image
            %       TestImg: [M x N x 3] Test Image
            % Output:
            %       Cmp:[N x M (double)] comparison matrix
            
            Cmp = cell(this.K,1);
            % extract hand written object
%             Target = TargetHandWrit(TargetImg,this.SqrSz,this.MaxJunction,this.MaxEdge);
%             Test = TestHandWrit(TestImg,this.SqrSz,this.MaxJunction,this.MaxEdge);
            
            % cluster target and test handwritten with K 
            [O,R,Eps,M] = Target.cluster(this.K);
            Test.clusterByCircle(M,O,R,Eps);
            
            % for each cluster
            for i=1:this.K
                % for each hand-written, return a Cmp object array
                MainCmpResult = Target.search(Test,i,this.Z,this.W,F);
                S = find(Target.ClustOrd == i);
                % add cluster id at the end of matrix
                ClusterID = ones(size(MainCmpResult,1),1) * S;
                MainCmpResult = cat(2,MainCmpResult,ClusterID);
                Cmp{i,1} = MainCmpResult;
            end
            
            % remove empty cells
            Empty = cellfun(@isempty,Cmp);
            Cmp = Cmp(~Empty,:);
            
            Cmp = cell2mat(Cmp);
            
        end        
        function Cmp = fullCompare(this,Target,Test)
            %FULLSEARCH search all of Train CCs with all Test CCs
            % Input:
            %       TargetImg:[M x N x 3] RGB Target Image
            %       TestImg: [M x N x 3] Test Image
            % Output:
            %       Cmp:[N x M (double)] comparison matrix
            
            
            
            % for each hand-written, return a Cmp object array
            Cmp = Target.fullSearch(Test,this.Z,this.W);
        end
        function showClusters(this,HWName,ClustNum)
            HW = this.returnByName(HWName);
            if ~exist('ClustNum','var')
                HW.displayClusters();
            else
                HW.displayClusters(ClustNum);
            end
        end
        function showCorrespondingClusters(this,i)
            %SHOWCORRESPONDINGCLUSTERS show corresponding clusters
            this.TrainHandWrit.displayClusters(i);
            this.MainHandWrit.displayClusters(i);
            this.Forge1HandWrit.displayClusters(i);
            this.Forge2HandWrit.displayClusters(i);
            this.Forge3HandWrit.displayClusters(i);
        end
        function hw = loadHandWrit(this,ID,Type)
            for i=1:length(this.G)
                if this.G(i).Writer == ID && strcmp(this.G(i).Type,Type)
                    hw = load(this.G(i).Address);
                    hw = hw.HW;
                    return
                end
            end
        end
    end
end