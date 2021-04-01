function Model = inkModel(Trace,CCImg,GrayImg,ListMask)
assert(isequal(size(CCImg),size(GrayImg)));
nPt = size(Trace,1);
Model = zeros(nPt,1);

% video writer for output result as a avi file
% v = VideoWriter(strcat('inkModel_',num2str(rand()),'.mkv'));
% open(v);

for i=1:nPt
    Pt = Trace(i,:);
    ClustImg = zeros(size(GrayImg),class(GrayImg));  
    Cluster = ClusterFilter(Pt,CCImg,ClustImg,ListMask);
%     imshow(cat(3,CCImg.*255,Cluster.*255,zeros(size(CCImg))));%show
%     frame = cat(3,CCImg.*255,Cluster.*255,zeros(size(CCImg)));%frame to write
%     writeVideo(v,frame);% write frame    
    Cut = Cluster .* GrayImg;
    Cut = Cut(Cut>0);
    Model(i) = mean(Cut(:));
end

% close(v); % close avi file

end
function ClustImg = ClusterFilter(Pt,CCImg,ClustImg,ListMask)
    % find maximum filter inside of connected component
    
    [Row,Col] = size(CCImg);
    CCSize = [Row,Col];
    
    Cut = zeros(CCSize,class(ClustImg));
    
    % filter counter
    f = 0; 
    
    % while all of filter is inside of handwrit
    while isequal(ClustImg,Cut)
        f = f + 1;
        FiltPts = ListMask{f};
        FiltPts(:,1) = FiltPts(:,1) + Pt(1);
        FiltPts(:,2) = FiltPts(:,2) + Pt(2);
        
        % remove points outside of image
        RowIndex = and(FiltPts(:,1) > 0 , FiltPts(:,1) <= Row);
        ColIndex = and(FiltPts(:,2) > 0 , FiltPts(:,2) <= Col);
        Index = and(RowIndex,ColIndex);
        FiltPts = FiltPts(Index,:);
        
        % convert filter point to index
        FiltIndex = FiltPts(:,1) + ((FiltPts(:,2)-1).*Row);
        
        ClustImg(FiltIndex) = 1;
        Cut = and(ClustImg,CCImg);
        
    end
    
        % draw last filter
        if f > 1
            f = f - 1;
        elseif f <= 1
            f = 1;
        end
        
        if f == 0
            disp('error');
        end
        FiltPts = ListMask{f};
        FiltPts(:,1) = FiltPts(:,1) + Pt(1);
        FiltPts(:,2) = FiltPts(:,2) + Pt(2);
        
        % remove points outside of image
        RowIndex = and(FiltPts(:,1) > 0 , FiltPts(:,1) < Row);
        ColIndex = and(FiltPts(:,2) > 0 , FiltPts(:,2) < Col);
        Index = and(RowIndex,ColIndex);
        FiltPts = FiltPts(Index,:);
        
        % convert filter point to index
       FiltIndex = FiltPts(:,1) + ((FiltPts(:,2)-1).*Row);
        
        ClustImg(FiltIndex) = 1;
end
