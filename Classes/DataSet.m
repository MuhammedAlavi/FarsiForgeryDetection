classdef DataSet
    %DATASET for our goal strokes
    properties(Constant = true)
        Y = LoadStrokeData('Y');
        N = LoadStrokeData('N');
    end
    
    methods
        function this = DataSet()
        end
        function s = sampleData(this,Stroke,nSample,ExceptionID)
            %TAKESAMPLE, choose a sample with nSample members.
            % Input:
            %       Stroke: [1 x 1 (char)] type of stroke.
            %       nSample: [1 x 1 (int)] number of samples
            %       ExceptionID: [vector (int)] seperate sample with
            %       ExceptionID in it
            % Output:
            %       s:[1 x nSample (struct)] sample with nSample member
            
            switch Stroke
                case 'Y'
                    RandIndex = sampleDataHelper(this.Y,nSample,ExceptionID);
                    s = this.Y(RandIndex);
                case 'N'
                    RandIndex = sampleDataHelper(this.N,nSample,ExceptionID);
                    s = this.N(RandIndex);
                otherwise
                    error('there is no such data-set!');
            end
        end
    end    
end

% return random index 
function RandIndex = sampleDataHelper(Stroke,nSample,ExceptionID)
% take sample of data-set
Len = length(Stroke);

IndexList = zeros(Len,1);
c = 0;
% find index that is not in exception set
for i=1:Len
    if Stroke(i).ID == ExceptionID
        continue
    end
    c = c+1;
    IndexList(c) = i;
end

% return random selection of stroke set
nSample = min([nSample,Len]);
RandIndex = randperm(Len,nSample);
end
