classdef USkillSenario < ComparisonBase
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N = 10; % [scalar (int)] number of unskilled forgery samples   
    end
    
    methods
        function this = USkillSenario()
            this@ComparisonBase();
            
            % CREATE A SENARIO FOR SKILLED FORGERY
            this.LookupTable = cell(this.nSampleIDValues,2);
            for i=1:this.nSampleIDValues
                this.LookupTable{i,1} = i;
                Temp = cell(this.N+2,2);
                R = randomValue(1:this.nSampleIDValues,this.N,i);
                
                % preallocate T2 and M hand-written
                Temp(1,:) = {i,'T2'};
                Temp(2,:) = {i,'M'};
                for j=1:this.N
                    Temp{j+2,1} = R(j);
                    Temp{j+2,2} = 'M';
                end
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
                    [Weight,ACC]= this.Classify();
                    
                    % strore result
                    this.ClassificationResult(c).Epoch = this.CurEpoch;
                    this.ClassificationResult(c).WinSize = this.CurWin;
                    this.ClassificationResult(c).Weight = Weight;
                    this.ClassificationResult(c).ACC = ACC;
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
                    
                    % best features to select
                    this.ClassificationResult(c).BestFeature = this.BestFeatures;
                    this.ClassificationResult(c).BestFeatureACC = this.BestFeaturesACC;
                    this.ClassificationResult(c).GAData = GAData;
                    this.ClassificationResult(c).PSOData = PSOData;
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

function r = randomValue(Range,Num,Exception)
% random number inside of Range except Exception values
% Input:
%       Range:[vector (double)] values
%       Num:[scalar (int)] number of random numbers
%       Exception:[vector (double)] values that should remove in Range
% Output:
%       r:[vector (int)] random values
if ~exist('Exception','var')
    Exception = [];
end

Range = setdiff(Range,Exception);
L = length(Range);
I = randperm(L,Num);
r = Range(I);
end