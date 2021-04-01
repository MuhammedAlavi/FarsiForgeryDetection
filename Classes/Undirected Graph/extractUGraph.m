function UG = extractUGraph(CC,BW)
%EXTRACTUGRAPH extract array of ugraph for each connected component
% Input:
%       CC: connected component structure of image
% Output:
%       UG: undirected graph array

UG = UGraph;
UG(CC.NumObjects) = UG;

for i=1:CC.NumObjects
    ThumbNail = ind2img(CC.PixelIdxList{i},BW);
    T = bwmorph(ThumbNail,'Thin',inf);
    J1 = bwmorph(T,'endpoints');
    J2 = bwmorph(T,'branchpoints');
    J = or(J1,J2);
    UG(i) = trace_graph(T,J,ThumbNail);
end


end

