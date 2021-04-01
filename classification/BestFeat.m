classdef BestFeat < handle
    %BESTFEAT is a class for all classification process
    % at the end of work we expect two value to get out of our 
    % classification process. 1. BEST FEATURE SET  2. BEST ACCURACY
    
    properties(SetAccess = private)
        % DATA VALUES
        Data
        Label
        
        % CLASSIFIER VALUES
        ClassifierRepeat = 5;      %[scalar (int)] number of repeating classifier funtions
        options = optimoptions('ga','display','iter','FunctionTolerance',1e-2);
        
        % HEURISTICS PARAMTERS
        BestFeatures 
        BestFeatureACC
    end
    properties(Dependent = true)
        nVector
        nFeature
        
    end
    methods
        % CONSTRUCTOR
        function this = BestFeat(Data,Label,R,options)
            this.Data = Data;
            this.Label = Label;
            if exist('R','var')
                if ~isempty(R)
                    this.ClassifierRepeat = R;
                end
            end
            if exist('option','var')
                if ~isempty(options)
                    this.options = options;
                end
            end
        end
        
        % GET FUNCTINOS
        function f = get.nFeature(this)
            f = size(this.Data,2);
        end
        function v = get.nVector(this)
            v = size(this.Data,1);
        end
        function [bFeat,bAccuracy] = run(this,FName)
            this.bestFeatures(FName)
            bFeat = this.BestFeatures;
            bAccuracy = this.BestFeatureACC;
        end
    end
    methods(Access = private)
        
        % CLASSIFIERS
        function A = simpleTree(this,Weight)
            Weight = logical(Weight);
            Temp = zeros(this.ClassifierRepeat,1);
            D = [this.Data,this.Label];
            for i=1:this.ClassifierRepeat
                [~,Temp(i)] = trainClassifierSTREE(D,Weight);
            end
            A = mean(Temp(:));
            A = -A;
        end
        function A = medianGSVM(this,Weight)
            Weight = logical(Weight);
            Temp = zeros(this.ClassifierRepeat,1);
            D = [this.Data,this.Label];
            for i=1:this.ClassifierRepeat
                [~,Temp(i)] = trainClassifierSTREE(D,Weight);
            end
            A = mean(Temp(:));
            A = -A;
        end
        
        % HEURISTICS
        function  bestFeatures(this,Func)
            %BESTFEATURES use Genetic algorithm to find best possible
            % feature
            % Input:
            %       this:[current object]
            %       Fun: [(string)] classifier function
            % Output:
            %       X:[vector (logical)] best features
            %       FVal:[scalar(double)] best accuracy
            
            switch Func
                case 'STREE'
                    f = @(W) this.simpleTree(W);
                case 'MGSVM'
                    f = @(W) this.medianGSVM(W);
                otherwise
                    error('there is no such name for classificaton function!');
            end
                    
            
            N = this.nFeature;
            LB = zeros(N,1);
            UB = ones(N,1);
            [X,FVal] = ga(f,N,[],[],[],[],LB,UB,[],1:N,this.options);
            this.BestFeatures = X;
            this.BestFeatureACC = FVal;
        end
        
        % DISPLAY
        function histObservation(this,ObsNum)
            %HISTOBSERVATION 
            I = this.Label;
            I = logical(I);
            TrueLabelData = this.Data(I,ObsNum);
            FalseLabelData = this.Data(~I,ObsNum);
            
            % plot
            subplot(2,2,1);
            hist(TrueLabelData);
            title('True Label Data');
            subplot(2,2,3);
            boxplot(TrueLabelData,'Orientation','horizontal');
            
            subplot(2,2,2);
            hist(FalseLabelData);
            title('False Label Data');
            subplot(2,2,4);
            boxplot(FalseLabelData,'Orientation','horizontal');
        end
        
    end
    
end


