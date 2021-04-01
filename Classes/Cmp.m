classdef Cmp
    %Cmp is all of information for comparison of two connected component(CC)
    
    properties
        
        % IS IT POSSIBLE TO COMPARE BOTH CC
        Comparable
        
        % BEZ COMPARISON PARAMETERS
        BezCx
        BezCy
        BezD
        
        % SIZE COMPARISON PARAMETERS
        Fraction
        S
        
        % INK MODEL PARAMETERS
        B
        %...
        
    end
    properties(Constant = true)
        % comparison threshold. 
        CmpThresh = [.5,1.5];
        StkType = {StrokeType.Simple,StrokeType.Complex};
        
        % bezier default parameters for dtw algorithm
        dtwMetric = 'euclidean';
        dtwMaxSamp = intmax('uint32');
        
    end
    
    methods
        function this = Cmp(TCC,ACCList,Z)
            %CMP constructor
            % Input:
            %       TCC[(ConnComp) object] target Connected Component
            %       ACCList[(ConnComp) objects] another Connected
            %       Z[scalar (int)] number of best match
            
            %       components
            
            % default constructor
            if ~exist('TCC','var') || ~exist('ACC','var')
                return
            end
            
            % now compare two CCs
            TAng = TCC.Ang;
            AAng = ACCList.Ang;
            
            % bezier angle dtw algorithm
            [this.BezD,~,~,this.BezCx,this.BezCy] = DTW(TAng,AAng);
            
            % Thumb size difference
            TCCRow = TCC.ThumbSize(1);
            ACCRow = ACCList.ThumbSize(1);
            TCCCol = TCC.ThumbSize(2);
            ACCCol = ACCList.ThumbSize(2);
            this.S = [abs(TCCRow-ACCRow),abs(TCCCol-ACCCol)];
            
            
        end
    end
    methods(Static = true)
        function YN = sizeMatch(TargetSz,AnotherSz,CmpThresh)
            
            MinSz = CmpThresh(1) .* TargetSz;
            MaxSz = CmpThresh(2) .* TargetSz;
            
            Cond1 = all(MinSz < AnotherSz);
            Cond2 = all(MaxSz > AnotherSz);
            
            YN = ~(Cond1 && Cond2);
        end
        function YN = stkTypeMatch(TargetST,AnotherST,StkTypes)
            N = length(this.StkType);
            for i=1:N
                if TargetST ~= this.StkType(i)
                    if TargetST ~= this.StkType
                        YN = false;
                    end
                end
            end
            if ismember(TargetST,StkTypes)&&ismember(AnotherST,StkTypes)
                YN = false;
            else
                YN = true;
            end
        end
    end
end



