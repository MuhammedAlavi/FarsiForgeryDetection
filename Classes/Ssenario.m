classdef Ssenario < ComparisonBase
    %SKILLEDFORGECLASSIFICATION skilled forgery senario.
    
    methods
        %CONSTRUCTOR
        function this = Ssenario()
            this@ComparisonBase();
            
            % CREATE A SENARIO FOR SKILLED FORGERY
            this.LookupTable = cell(this.nSampleIDValues,2);
            for i=1:this.nSampleIDValues
                this.LookupTable{i,1} = i;
                Temp = {i,'T2';i,'M';i,'F1';i,'F2';i,'F3'};
                this.LookupTable{i,2} = Temp;
            end
            
            this.CmpResult = struct('ID',[],'Cmp',[]);
            this.CmpResult(this.I) = this.CmpResult;
            
            this.ClassificationResult = struct('Epoch',[],'WinSize',[],'Weight',[],'ACC',[]);
            this.ClassificationResult(this.nEpoch * this.nWin) = this.ClassificationResult;
        end   
        function run(this)
            % initialize
            this.initIteration();
            
            c = 1;
            % run comparison
            while ~this.EndEpoch
                this.resetWindow();
                while ~this.EndWindow
                    
                    % compare samples according to lookup table
                    this.comparison();
                    
                    % find best accuracy
                    [Weight,ACC,GAData,PSOData]= this.Classify();
                    
                    % strore result
                    this.ClassificationResult(c).Epoch = this.CurEpoch;
                    this.ClassificationResult(c).WinSize = this.CurWin;
                    this.ClassificationResult(c).Weight = Weight;
                    this.ClassificationResult(c).ACC = ACC;
                    this.ClassificationResult(c).GAData = GAData;
                    this.ClassificationResult(c).PSOData = PSOData;
                    
                    % best features to select
                    this.ClassificationResult(c).BestFeature = this.BestFeatures;
                    this.ClassificationResult(c).BestFeatureACC = this.BestFeaturesACC;
                    c = c + 1;
                    % next window size
                    this.nextWindow();
                    this.resetSampleID();
                end
                this.nextEpoch();
            end
        end
        function comparison(this)
            % show steps
            this.showTheSteps();
           % call compare function from super class ComparisonBase
            D = this.compare();
            for i=1:size(D,1)
                this.CmpResult(i).ID = D{i,1};
                this.CmpResult(i).Cmp = D{i,2};
            end
        end  
        function generate(this)
            % generate the data for a weight vector 'W' and 
            % features 'F'//
            % initialize
            
            this.initIteration();
            
            c = 1;
            % run comparison
            while ~this.EndEpoch
                this.resetWindow();
                while ~this.EndWindow
                    
                    % compare samples according to lookup table
                    this.comparison();
                    
                    % acquire the weights and features
                    W = this.ClassificationResult(c).Weight;
                    F = this.ClassificationResult(c).BestFeature;
                    
                    % find best accuracy
                    Data= this.createMatrix(W,F);
                    
                    % strore result
                    this.ClassificationResult(c).Data = Data;
                    c = c + 1;
                    % next window size
                    this.nextWindow();
                    this.resetSampleID();
                end
                this.nextEpoch();
            end
        end
    end
end
