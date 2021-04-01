classdef IdentifyIteration < handle
    %SAMPLE a class for each iteration of algorithm
    properties
        
        % raw images
        TrainImg       % train image
        Forge1Img       % forgery images
        Forge2Img
        Forge3Img
        MainImg        % main image( owner handwriting)
        
        % hand-written objects
        TrainHandWrit
        Forge1HandWrit
        Forge2HandWrit
        Forge3HandWrit
        MainHandWrit
    end
    properties(Constant = true)
%         DS = DataSet();
        SqrSz = 20;
        MaxJunction = 5;
        MaxEdge = 5;
        K = 10; 
    end
    methods
        %INITIALIZATION
        function this = VerifyIteration(Address)
            if ~exist('Address','var')
                error('address of images are essential!');
            end
            [this.TrainImg,...
                this.MainImg,...
                this.Forge1Img,...
                this.Forge2Img,...
                this.Forge3Img] = readFile(Address);
        end
        function convert2HandWrit(this)
            %CONVERT2HANDWRIT convert images into handwritten objects. 
            % Input:
            %       this: object
            
            SS = this.SqrSz;
            MJ = this.MaxJunction;
            ME = this.MaxEdge;
            
            % train handwritten
            this.TrainHandWrit = TargetHandWrit(this.TrainImg,SS,MJ,ME);
            this.TrainImg = [];
            
            % capture handwritten object of forgery images
            this.Forge1HandWrit = TestHandWrit(this.Forge1Img,SS,MJ,ME);
            this.Forge1Img = [];
            this.Forge2HandWrit = TestHandWrit(this.Forge2Img,SS,MJ,ME);
            this.Forge2Img = [];
            this.Forge3HandWrit = TestHandWrit(this.Forge3Img,SS,MJ,ME);
            this.Forge3Img = [];
            
            % main image(owner handwritten)
            this.MainHandWrit = TestHandWrit(this.MainImg,SS,MJ,ME);
            this.MainImg = [];
        end
        function inkModelExtraction(this,StrkType)
            %INKMODELEXTRACTION trigger handwrittens inkModelExtraction
            %methods
            I = {'Train','Main','Forge1','Forge2','Forge3'};
            for i = I
                HW = this.returnByName(i{1});
                
                % return index of these stroke types
                Index  = HW.getCCWithType(StrkType);
                % extract their ink model
                HW.inkModelExtraction(Index)
            end
        end
        function cluster(this,PrimeK)
            % cluster all of stroke of handwritten. clusterHelper adds a
            % new property(.clust) to data features structure.
            
            if ~exist('PrimeK','var')
                PrimeK = this.K;
            end
            [O,R,Epsilon] = this.TrainHandWrit.cluster(PrimeK);
            
            this.MainHandWrit.clusterByCircle(O,R,Epsilon);
            this.Forge1HandWrit.clusterByCircle(O,R,Epsilon);
            this.Forge2HandWrit.clusterByCircle(O,R,Epsilon);
            this.Forge3HandWrit.clusterByCircle(O,R,Epsilon);
        end
        function [D] = searchForMatch(this,i,CmpThresh)
            
            
            HWNames = {'Main','Forge1','Forge2','Forge3'};
            TargetHW = this.returnByName('Train');
            
            % 
            D = cell(TargetHW.nCC,2);
            
            for h=HWNames% for all hand-writtens
                HW = this.returnByName(HWNames{i});
                for i=1:TargetHW.CC
                    % check comparison condition
                    if 1==2
                    end
                    for j=1:HW.nCC
                        
                        %check comparison condition
                        if 2==2
                        end
                        
                    end
                end
            end
        end
        %GET FUNCTIONS
        function [HandWritObj] = returnByName(this,Name)
            %RETURNBYNAME, return handwritten object and features by name
            % of object.
            % Input:
            %       Name: [1 x N (char)] name of handwritten.
            %       Name can be:
            %           {'Train','Main','Forge1','Forge2','Forge3'}
            % Output:
            %       HandWritObj: [HandWrit object] 
            %       Features: [ feature data]
            switch Name
                case 'Train'
                    HandWritObj = this.TrainHandWrit;
                case 'Main'
                    HandWritObj = this.MainHandWrit;
                case 'Forge1'
                    HandWritObj = this.Forge1HandWrit;
                case 'Forge2'
                    HandWritObj = this.Forge2HandWrit;
                case 'Forge3'
                    HandWritObj = this.Forge3HandWrit;
                otherwise
                    error('bad input name!');
            end
        end
        %DISPLAY
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
    end
end

