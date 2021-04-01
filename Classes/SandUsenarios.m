classdef SandUsenarios < ComparisonBase
    %SKILLEDANDUNSKILLEDCLASSIFICATION combine both skilled and unskilled
    %senarios
    properties(Constant = true)
        nSkilled = 2;    %[scalar (int)] number of skilled hand-written
        nUnSkilled = 2;  %[scalar (int)] number of unskilled hand-written
    end
    
    methods
        function this = SandUsenarios()
            this@ComparisonBase();
            
            % CREATE A SENARIO FOR SKILLED FORGERY
            this.LookupTable = cell(this.nSampleIDValues,2);
            for i=1:this.nSampleIDValues
                MainTags = {i,'T2';i,'M'};
                this.LookupTable{i,1} = i;
                R1 = randomValue(1:this.nSampleIDValues,this.nUnSkilled,i);
                R2 = randperm(this.nForgeryHandWrit,this.nUnSkilled);
                
                % convert numbers to Tag
                SkilledTags = cell(length(R2),2);
                for l=1:length(R2)
                    SkilledTags{l,1} = i;
                    if R2(l) == 1
                        SkilledTags{l,2} = 'F1';
                    elseif R2(l) == 2
                        SkilledTags{l,2} = 'F2';
                    elseif R2(l) == 3
                        SkilledTags{l,2} = 'F3';
                    end
                end
                
                % preallocate T2 and M hand-written
                USkilledTags = cell(length(R1),2);
                for j=1:length(R1)
                    USkilledTags{j,1} = R1(j);
                    USkilledTags{j,2} = 'M';
                end
                this.LookupTable{i,2} = cat(1,MainTags,SkilledTags,USkilledTags);
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
                    [Weight,ACC,GAMat,PSOMat]= this.Classify();
                    
                    % strore result
                    this.ClassificationResult(c).Epoch = this.CurEpoch;
                    this.ClassificationResult(c).WinSize = this.CurWin;
                    this.ClassificationResult(c).Weight = Weight;
                    this.ClassificationResult(c).ACC = ACC;
                    % best features to select
                    this.ClassificationResult(c).BestFeature = this.BestFeatures;
                    this.ClassificationResult(c).BestFeatureACC = this.BestFeaturesACC;
                    % save the data
                    this.ClassificationResult(c).GAData = GAMat;
                    this.ClassificationResult(c).PSOData = PSOMat;
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
