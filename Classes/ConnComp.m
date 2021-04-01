classdef ConnComp < handle
    %CONNCOMP connected component
    properties(SetAccess = private)
        % BASE PROPERTIES
        ID          % [scalar (int)] identifier of ConnComp
        ImgSize     % [1 x 2 (int)] size of original image (main image that
        % CCs has cropped
        PixelIdx    % [Z x 1 (int)] index of pixle in Image
        
        % DERIVED PROPERTIES
        ST          % [1 x 1 (Stroketype)] stroke type
        Trace       % [1 x 1 (int)] trace of CC
        Bez         % [P x 2 (int) ] bezier curve
        Ang         % [P x 1 (int) ] angle of direction
        IM          % [K x 1 (int) ] ink model
        ThumbSize   % [1 x 2 (int) ] size of CC image
    end
    properties(Dependent = true)
        Img         % [N x M (bool)] image of connected component
        COM         % [1 x 2 (int) ] cneter of mass
        COMProp     % [1 x 2 (double)] center of mass proportion
        Slant       % [scalar (double)] slant of CC image
        BProjHist   % [vector (int)] buttom projection histogram
        SProjHist   % [vector (int)] side projection histogram
        % %         Thin        % [N x M (bool)] thinned image
        % %         Endpoint    % [N x M (bool)] endpoint in morphological operation
        % %         Branchpoint % [N x M (bool)] branch point in morpholocial operation
    end
    properties(Constant = true)
        ListMask = loadListMask();
        BezHelper = loadBezierHelperMatrix();
    end
    methods
        
        %INITIZLIZE 
        function this = ConnComp(Index,ImgSize,ID)
            if ~exist('Index','var') || ~exist('ImgSize','var')
                return
            end
            this.PixelIdx = Index;
            this.ImgSize = ImgSize;
            this.ID = ID;
        end
        function initConnComp(this,Skew,MaxJunction,MaxEdge)
            %CREATECONNCOMP create other propeties of ConnComp class
            % Input:
            %       Index: index of image
            %       GrayImg: gray scale image
            %       MaxEdge: maximum edge in the graph of image
            % Output:
            %       UG: undirected graph
            %       Bez: bezier curve
            %       IM: ink model
            %       ST: Stroke type
            
            % extract both BW and Gray image of CC
            CCImg = ind2img(this.PixelIdx,false(this.ImgSize));
            this.ThumbSize = size(CCImg);
              
            % check CCImg to not have a hole
            if hasLoop(CCImg)
                this.unknownObject();
                return
            end
            try
                % thinned image, junction points to create an undirected graph
                T = bwmorph(CCImg,'thin',inf);
%                 J1 = bwmorph(T,'endpoints');
%                 J2 = bwmorph(T,'branchpoints');
%                 J = or(J1,J2);
                J = extractJunctions(T);
                % an upper bound for count of junction points
                JIndex = sum(J(:));
                if JIndex > MaxJunction || JIndex == 1 % at least one junction point
                    % return unknown object
                    this.unknownObject();
                    return
                end
                
                % create graph
                UG = trace_graph(T,J,CCImg);
                
                % clean graph by circular critrion
                UG.clean_skeleton;
                
                % make a trace on undirected graph
                [Trc,StrokeType] = UG.makeTrace(MaxEdge);
                
                % after cleaning skeleton some of points are not integer
                Trc = floor(Trc);
                this.Trace = Trc;
                this.ST = StrokeType;
                
                % if trace exist(-1 means there is something with Trace)
                if Trc == -1
                    % if there is no trace for graph
                    this.unknownObject();
                    return
                end
                
                Y = Trc(:,1);
                X = Trc(:,2);
                
                % image has skew, to eleminate this skew, we rotate each CC around
                % center of gravity.
                c = regionprops(CCImg,'centroid');
                c = c.Centroid;
                B = bezier(Y,X,[],this.BezHelper);
                B = slantCorrection(B(:,1),B(:,2),c,Skew);
                this.Bez = B;
                
                Vec = diff(B);
                Angle = atan2d(Vec(:,2),Vec(:,1));
                % normalize angle
                sign = Angle < 0;
                Angle = abs(Angle);
                GrThan90 = Angle > 90;
                Angle(GrThan90) = 180 - Angle(GrThan90);
                Angle(sign) = -Angle(sign);
                this.Ang = Angle;
            catch Exc
                warning(Exc.message);
                this.unknownObject();
                return
            end
        end
        function initInkModel(this,GrayImg)
            % for the purpose of performance, ink model is called after
            % resemble strokes, extracted.
            if this.ST == StrokeType.Unknown
                this.IM = -1;
            else
                GrayCCImg = ind2img(this.PixelIdx,GrayImg);
                this.IM = inkModel(this.Trace,this.Img,GrayCCImg,this.ListMask);
            end
        end
        function unknownObject(this)
            %UNKOWNOBJECT return an empty object
            this.Trace = -1;
            this.ST = StrokeType.Unknown;
            this.Bez = -1;
            this.Ang = -1;
        end
        
        %GET FUNCTIONS
        function I = get.Img(this)
            I = ind2img(this.PixelIdx,false(this.ImgSize));
        end
        function c = get.COM(this)
            % center of gravity
            I = this.Img;
            c = regionprops(I,'centroid');
            c = c.Centroid;
        end
        function f = get.COMProp(this)
            f = this.ThumbSize./this.COM;
        end
        function t = get.Slant(this)
            %SLANT of CC image
            
            t = slantAng(this.Img);
        end
        function h = get.BProjHist(this)
            h = sum(this.Img,2);
        end
        function h = get.SProjHist(this)
            h = sum(this.Img,1);
            h = h(:);
        end
        %DISPLAY
        function displayBez(this)
            plot(-this.Trace(:,1),this.Trace(:,2),'k');hold on;
            scatter(-this.Trace(:,1),this.Trace(:,2),'k','fill');
             scatter(-this.Bez(1,1),this.Bez(1,2),'gr','fill');
            plot(-this.Bez(:,1),this.Bez(:,2),'blue');hold on;
            scatter(-this.Bez(:,1),this.Bez(:,2),'blue','fill');
%             scatter(this.Bez(1:10:end,1),this.Bez(1:10:end,2),'red','fill');
            scatter(-this.Bez(1,1),this.Bez(1,2),'gr','fill');
%             legend('Plot','Points','MileStone','StartPoint');
        end
    end
end

function YN = hasLoop(Img)
%HASLOOP, check Img and return true if it has a hole.
% Input:
%       Img: [N x M (logical)] black and white image(black == true)
% Output:
%       YN: [1 x 1 (logical)] true if loop exist else false

% covert border of Img to not false connected component
NewImg = padarray(Img,[5,5]);
cc = bwconncomp(~NewImg);
if cc.NumObjects > 1
    YN = true;
else
    YN = false;
end
end