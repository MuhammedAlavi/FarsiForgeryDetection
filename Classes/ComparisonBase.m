classdef ComparisonBase < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant = true, Hidden = true)
        % CONSTANT SCALAR
        Z = 1;
        SqrSz = 20;
        MaxJunction = 5;
        MaxEdge = 5;
        nEpoch = 1;
        
        % CONSTANT VECTOR
        WinSz = [2,4,8,16,-1];
        K = 10;
        SampleIDValues = 1:62;
        
        % HEURISTIC ALGORITHM CONSTRAINTS
        MinWeight = 0;
        MaxWeight = 1;
        FunctionRepeat = 5;
        Classifier = 'STREE';
        
        % DATASET PARAMTERS
        nForgeryHandWrit = 3;    % number of forgery hand-written
    end
    properties(SetAccess = private)
        CurSampleID     %[ scalar (int)] index of current Target in lookup table
        t               %[ 1 x 1 (Iteration)]
        Clusters        %[ 1 x this.I (cell)] cluster ID for each hand writ ID
        BestFeatures    %[vector (logical)] feature vector
        BestFeaturesACC %[scalar (double)] best accuracy for best feature
    end
    properties(SetAccess = protected)
        Cmp             %[ N x nEpoch (cell)] result of comparison for each epoch
        LookupTable     %[ Z x 2 (cell)] table of target and test IDs
        
        CurWin          %[ scalar (int)] current window size
        CurWinNum       %[ scalar (int)] current window number
        
        CurEpoch        %[ scalar (int)] current epoch number
        
        CmpResult       %[struct] comparison result
        ClassificationResult % [struct] result of classification after optimization
        
        OptimFlag = true   %[(bool)] flag to set optimization algorithm
    end
    properties(Dependent = true)
        
        % lookup table counters
        I       % number of Target identifiers
        J       % number of Test identifiers
        
        nWin
        nK
        nSampleIDValues
        
        % is it end?
        EndEpoch
        EndWindow
        EndSampleID
    end
    methods
        % CONSTRUCTOR
        function this = ComparisonBase()
            % current number of window size and sample hand-written number
            this.CurWinNum = 1;
            this.CurWin = this.WinSz(this.CurWinNum);
            this.CurSampleID = 1;
            
            % current solution values
            this.CurEpoch = 1;
            
            % initialize iteration object
            ME = this.MaxEdge;
            MJ = this.MaxJunction;
            this.t = Iteration(this.CurWin,this.SqrSz,this.Z,ME,MJ);
        end
        % MAIN FUNCTION FOR COMPARISON
        function D = compare(this)
            if isempty(this.LookupTable)
                error('Define a senario to start Comparison!');
            end
            
            % set iteration object
            this.setIteraion();
            
            % comparison preallocation
            D = cell(this.I,2);
            while ~this.EndSampleID
                
                Temp = cell(this.J,2);
                TargetID = this.LookupTable{this.CurSampleID,1};
                for j = 1:this.J
                    TestID = this.LookupTable{this.CurSampleID,2}{j,1};
                    TestType = this.LookupTable{this.CurSampleID,2}{j,2};
                    % comparison of current smaple ID with another Test ID
                    Temp{j,2} = this.t.fullCompareWithID(TargetID,TestID,TestType);
                    Temp{j,1} = this.tag(TargetID,TestType,TestID);
                end
                D{this.CurSampleID,1} = TargetID;
                D{this.CurSampleID,2} = Temp;
                this.nextSampleID
                % show where is algorithm?
                this.showTheSteps();
            end
        end
        % INITIALIZATION AND NEXT FUNCTIONS
        function initIteration(this)
            this.resetEpoch();
            this.resetWindow();
            this.resetSampleID();
        end
        function setIteraion(this)
            ME = this.MaxEdge;
            MJ = this.MaxJunction;
            this.t = Iteration(this.CurWin,this.SqrSz,this.Z,ME,MJ);
        end
        function nextWindow(this)
            
            % change window size to next
            if this.CurWinNum >= this.nWin
                this.CurWinNum = this.CurWinNum + 1;
            else
                this.CurWinNum = this.CurWinNum + 1;
                this.CurWin = this.WinSz(this.CurWinNum);
            end
        end
        function nextEpoch(this)
            this.CurEpoch = this.CurEpoch + 1;
        end
        function nextSampleID(this)
            this.CurSampleID = this.CurSampleID + 1;
        end
        function resetWindow(this)
            this.CurWinNum = 1;
            this.CurWin = this.WinSz(this.CurWinNum);
        end
        function resetEpoch(this)
            this.CurEpoch = 1;
        end
        function resetSampleID(this)
            this.CurSampleID = 1;
        end
        % GET FUNCTIONS
        function n = get.nWin(this)
            n = length(this.WinSz);
        end
        function n = get.nK(this)
            n = numel(this.K);
        end
        function n = get.nSampleIDValues(this)
            n = length(this.SampleIDValues);
        end
        function YN = get.EndEpoch(this)
            if this.CurEpoch > this.nEpoch
                YN = true;
            else
                YN = false;
            end
        end
        function YN = get.EndWindow(this)
            if this.CurWinNum > this.nWin
                YN = true;
            else
                YN = false;
            end
        end
        function YN = get.EndSampleID(this)
            if this.CurSampleID > this.I
                YN = true;
            else
                YN = false;
            end
        end
        function i = get.I(this)
            i = size(this.LookupTable,1);
        end
        function j = get.J(this)
            Temp = this.LookupTable{this.CurSampleID,2};
            j = size(Temp, 1);
        end
        % DISPLAY
        function showTheSteps(this)
            % verbosely Display Steps
            clc;
            disp(strcat('current epoch: ',num2str(this.CurEpoch),'/',num2str(this.nEpoch)));
            disp(strcat('current Window:',num2str(this.CurWinNum),'/',num2str(this.nWin)));
            disp(strcat('current ID:', num2str(this.CurSampleID),'/',num2str(this.I)));
        end
    end
    methods (Access = protected)
        function Data = createMatrix(this,W,F)
            % CREATEMATRIX, create the data with weight W and 
            % features of 'F' //
            D = this.createDataMat(W);
            Data = D(:,[find(F),size(D,2)]);
        end
        function [A,FVal,GAMat,PSOMat] = Classify(this)
            
            % cluster each hand written in K clusters
            this.createCluster();
            
            % HEURISTIC FOR BEST FEATURE
            W = ones(1,this.K);
            Mat = this.createDataMat(W);
            
            % if not optimize
            if ~this.OptimFlag
                A = -1;
                FVal = -1;
                return
            end
            
            % find best features
            disp('genetic algorithm for best features!');
            options = optimoptions('ga','display','iter','MaxTime',3600);
            Obj = BestFeat(Mat(:,1:end-1),Mat(:,end),this.FunctionRepeat,options);
            [this.BestFeatures,this.BestFeaturesACC] = Obj.run(this.Classifier);
            GAMat = Mat(:,[find(this.BestFeatures),size(Mat,2)]);
            
            % objective function
            f = @(Weight)this.weightedClassification(Weight);
            

            % HEURISTIC PROPERTIES
            N = this.K;
            LB = ones(N,1) .* this.MinWeight;
            UB = ones(N,1) .* this.MaxWeight;
%             IntCon = 1:N;
            options = optimoptions('particleswarm','display','iter','MaxTime',7200);
            disp('PC algorithm for best weight of clusters');
            [A,FVal] = particleswarm(f,N,LB,UB,options);
%             [A,FVal] = ga(f,N,[],[],[],[],LB,UB,[],IntCon,options);
            PSOMat = this.weightedClassificationData(A);
        end
        function A = weightedClassification(this,W)
            % create data matrix
            Mat = this.createDataMat(W);
            % find best matched features
            f = this.BestFeatures;
            f = logical(f);
            
            Temp = zeros(this.FunctionRepeat,1);
            for i=1:this.FunctionRepeat
                [~,Temp(i)] = trainClassifierSTREE(Mat,f);
            end
            A = mean(Temp);
            A = -A;
        end
        function M = weightedClassificationData(this,W)
            % create data matrix
            Mat = this.createDataMat(W);
            % find best matched features
            f = this.BestFeatures;
            M = Mat(:,[find(f),size(Mat,2)]);
        end
        function t = tag(this,TargetID,TestType,TestID)
            if isa(this,'Ssenario')
                switch TestType
                    case 'T2'
                        t = -1;
                    case 'M'
                        t = 1;
                    case {'F1','F2','F3'}
                        t = 0;
                end
            elseif isa(this,'USkillSenario') || isa(this,'SandUsenarios')
                if TestID == TargetID
                    switch TestType
                        case 'T2'
                            t = -1;
                            return
                        case 'M'
                            t = 1;
                            return
                        case {'F1','F2','F3'}
                            t = 0;
                            return
                    end
                else
                    t = 0;
                end
            end
            
        end
        function D = createDataMat(this,W)
           
            % create feature matrix for each iteration
            D = cell(this.I,1);
            for i=1:this.I
                ClusterID = this.Clusters{i};
                D{i} = createFeatureMat(this.CmpResult(i).Cmp,ClusterID,W);
            end
            D = cell2mat(D);
        end
        function createCluster(this)
            % cluster each hand written ID in cmp matrix
            this.Clusters = cell(this.I,1);
            for i=1:this.I
                this.Clusters{i} = this.t.clusterWithID(this.CmpResult(i).ID,this.K);
            end
        end
    end
end

function Mat = createFeatureMat(DataPiece,ClusterID,W)

N = size(DataPiece,1);
Mat = cell(N,1);
for i=1:N
    Tag = DataPiece{i,1};
    AllVar = DataPiece{i,2};
    Mat{i} = featureExtraction(Tag,AllVar,ClusterID,W);
end
Mat = cell2mat(Mat);

% find base vector
I = Mat(:,end) == -1;
V = Mat(I,:);
for i=1:size(Mat,1)
    if i == find(I)
        continue
    else
        Mat(i,1:end-1) = Mat(i,1:end-1) - V(1:end-1);
    end
end
Mat = Mat(Mat(:,end)>=0,:);
end
