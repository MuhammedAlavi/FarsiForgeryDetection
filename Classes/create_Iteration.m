function Itr = create_Iteration(FolderAddress)
%CREATE_ITERATION create an Iteration object.
% Input:
%       FolderAddress: address of iteration images
% Output: 
%       Itr: iteration object to collect data

% capture address of images
d = dir(FolderAddress);

% predefine forgery images
ForgImg = cell(1,3);

for f=3:length(d)
    ImageFile = strcat(FolderAddress,'\',d(f).name);
    switch d(f).name
        case 'TRAIN.tiff'
            TrainImg = imread(ImageFile);
        case {'FORG_1.tiff','FORG_2.tiff','FORG_3.tiff'}
            P = strsplit(d(f).name,'_');
            P = P{2};
            % remove file extention
            P = strsplit(P,'.');
            P = str2double(P{1});
            ForgImg{P} = imread(ImageFile);
        case 'MAIN.tiff'
            MainImg = imread(ImageFile);
        otherwise
    end
end

Itr = Iteration(TrainImg,ForgImg,MainImg);

end

