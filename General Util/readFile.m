function [Train,Train2,Main,Forge1,Forge2,Forge3] = readFile(address)
%READFILE read files of image.
Main = imread(strcat(address,'\','MAIN.tiff'));
Main = Main(:,:,1:3);

Train = imread(strcat(address,'\','TRAIN.tiff'));
Train = Train(:,:,1:3);

Train2 = imread(strcat(address,'\','TRAIN2.tiff'));
Train2 = Train2(:,:,1:3);

for i=1:3
    Img = imread(strcat(address,'\','FORGE_',num2str(i),'.tiff'));
    Img = Img(:,:,1:3);
    switch i
        case 1
            Forge1 = Img;
        case 2
            Forge2 = Img;
        case 3
            Forge3 = Img;
    end
end

end

