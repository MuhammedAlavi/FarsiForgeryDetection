classdef Collection
    %COLLECTION is series of iteration.
    
    properties
        SeriesOfItr         % [N x 1 (Iteration)] all iterations
        N                   % [1 x 1 (int)] number of iteration
        TimeInfo            % (struct) time information
    end
    properties
        SqrSz = 20;
    end
    methods
        function this = Collection(n)
            this.N = n;
            this.SeriesOfItr = Iteration();
            this.SeriesOfItr(n) = Iteration();
            this.TimeInfo = struct('ReadFiles',[],'Convert2HW',[]);
        end
        function convert2HandWrit(this)
            for i=1:this.N
                tic;
                this.SeriesOfItr(i).convert2HandWrit(this.SqrSz);
                t = toc;
                this.TimeInfo(i).Convert2HW = t;
                disp(strcat('Convert to HandWrit: ',num2str(floor(i/this.N)),'%'));
            end
        end
    end
    
end

