%% create video




locdir = 'C:\Users\Casa\Desktop\Rafael\Tese\calcio_automatic\inputs\first_batch_01\processed\';
bmpfiles = dir(strcat(locdir, '*.bmp'));
sortedbmpfiles = natsortfiles({bmpfiles.name});

cellAvi = VideoWriter(fullfile(locdir,'teste.avi'));
open(cellAvi)

numofframes = size(dir([locdir '*.bmp']),1);

[rows, columns, numberOfColorChannels] = size(imread(strcat(locdir, sortedbmpfiles{1})));

for frameNumber = 1:numofframes
   im = imread(strcat(locdir, sortedbmpfiles{frameNumber}));
   newIm = imresize(im, [rows, columns]);

   img = im2double(newIm);
   writeVideo(cellAvi,img)
end


close(cellAvi)
clear cellAvi img inffiles

%% play video
vid = VideoReader(fullfile(locdir,'teste.avi'));

% Set up video figure window
videofig(vid.NumberOfFrames, @(frm) redraw(frm, vid));

% Display initial frame
redraw(1, vid);

pause
clear vid
