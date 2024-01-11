%% mp4 converter

clear, clc
close all

Path="PNG_folder";
PngNames="GridSearch";
Extension="png";

writerObj = VideoWriter('GridSearch.mp4','MPEG-4');
open(writerObj);

Dir=dir(append(Path,'\*',Extension));
totLen=size(Dir,1);

reverseStr = '';
for k = 1 : totLen
    name=append(Path,"\",PngNames,"_",num2str(k),".png");
    image = imread(name);
    writeVideo(writerObj, image);

    msg = sprintf('Processed %d/%d', k, totLen);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

fprintf("\n");
close(writerObj);

