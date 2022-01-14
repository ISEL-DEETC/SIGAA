function [excelname] = patternfindercell(localdir)
%%Read input frame and take the ROI

batchname = input('What batch to work? ','s');
batchname = ['\' batchname '\'];

locdir = strcat(localdir,batchname);

bmpfiles = dir(strcat(locdir, '*.bmp'));
sortedbmpfiles = natsortfiles({bmpfiles.name});
numofframes = size(dir([locdir '*.bmp']),1);

inffiles = dir([locdir '*.INF']);
timename=inffiles(1).name;

time = gettime(strcat(locdir,timename));

%% process ratio video frames





%% create video


cellAvi = VideoWriter(fullfile(locdir,'cellAvi.avi'));
open(cellAvi)


for frameNumber = 1:numofframes
    
   img = im2double(imread(strcat(locdir, sortedbmpfiles{frameNumber})));
   writeVideo(cellAvi,img)
end


close(cellAvi)

%% play video
vid = VideoReader(fullfile(locdir,'cellAvi.avi'));

% Set up video figure window
videofig(vid.NumberOfFrames, @(frm) redraw(frm, vid));

% Display initial frame
redraw(1, vid);

%%


firstframe = im2double(imread(strcat(locdir, sortedbmpfiles{1})));
imshow(firstframe)

str = input('How many cells to study? ')
close all;

excelname=input('Name the output excel file with extension (.xls or .xlsx)? ','s');
%excelname=[excelname '.xlsx'];


%%
I = im2double(imread(strcat(locdir, sortedbmpfiles{1})));


J={};
rect={};
output={};
output{end+1}=time';
corr_offset_before = [];


for crop=1:str
    [J{end+1}, rect{end+1}] = imcrop(I);
end
I = rgb2gray(I);
close all;
f1 = figure;
%f2 = figure;
pixelVals = [];



%%

%writematrix(time,ExcelFile);



for eachcell=1:length(J)
    %clf(f2,'reset');
    for frameNumber=1:numofframes
      
        K_color = im2double(imread(strcat(locdir, sortedbmpfiles{frameNumber})));
        K = rgb2gray(K_color);
        
       
        %% Calculate normalized cross correlation
        % This method works really well for this case if the images are not
        % blurred
        %only calculate correlation between frames and main frame once(corr_offset_total),
        %only on first cell for performance
        if(eachcell==1)
            t = normxcorr2(I, K);
            [max_t, imax] = max(abs(t(:)));
            %Offset between the images found by correlation
            [ytpeak, xtpeak] = ind2sub(size(t),imax(1));
            corr_offset_total{frameNumber} = [(xtpeak-size(I,2)) (ytpeak-size(I,1))];
        end
       
       
        c = normxcorr2(rgb2gray(J{1,eachcell}), K);
        %% find peak correlation and find the offset
       
        [max_c, imax] = max(abs(c(:)));
        %Offset between the images found by correlation
        [ypeak, xpeak] = ind2sub(size(c),imax(1));
        corr_offset = [(xpeak-size(J{1,eachcell},2)) (ypeak-size(J{1,eachcell},1))];
       
        if (frameNumber == 1)
            corr_offset_before = corr_offset;
            corr_offset_total_before = corr_offset_total{frameNumber};
        else
            normcell = norm(corr_offset-corr_offset_before);
            normframe = norm(corr_offset_total{frameNumber}-corr_offset_total_before);
            if (normcell>normframe+10)
                corr_offset = corr_offset_before;
            else
                corr_offset_before = corr_offset;
                corr_offset_total_before = corr_offset_total{frameNumber};
            end
        end
       
       
       
        %% display best match
       
        figure(f1), imshow(K_color); hold on;
        set(f1, 'NumberTitle', 'off', ...
            'Name', sprintf('Frame: %s', int2str(frameNumber)));
        rectangle('position',[corr_offset(1) corr_offset(2) rect{1,eachcell}(3) rect{1,eachcell}(4)],...
            'edgecolor','g','linewidth',1);
       
        cell = imcrop(K,[corr_offset(1) corr_offset(2) rect{1,eachcell}(3) rect{1,eachcell}(4)]);
       
        %figure(f2), imshow(cell); hold on;
        %% center crop circle
        %if all we want is to probe for the circle
        center = floor((size(cell)/2+.5));
        %viscircles(center,5)
       
        %%
        if (frameNumber == 1 )%||not(isequal(size(cell), size(cellbefore))))
            bw = imbinarize(cell);
            bw = bwareaopen(bw,50);
            se = strel('disk',5);
            bw = imclose(bw,se);
            bw = imfill(bw,'holes');
            [B,L] = bwboundaries(bw,'noholes');
            %imshow(label2rgb(L,@jet,[.5 .5 .5]))

        end
            cellbefore = cell;
       
        try

                            pixelvalue = regionprops(L,cell,'MeanIntensity');
        catch 
            clf(f1,'reset');
            %clf(f2,'reset');
            disp('The cell has moved partially out of the frame. Please crop another region.')
            [J{1,eachcell},rect{1,eachcell}] = imcrop(rgb2gray(im2double(imread(strcat(locdir, sortedbmpfiles{frameNumber})))));
            c = normxcorr2(rgb2gray(J{1,eachcell}), K);
            %% find peak correlation and find the offset

            [max_c, imax] = max(abs(c(:)));
            %Offset between the images found by correlation
            [ypeak, xpeak] = ind2sub(size(c),imax(1));
            corr_offset = [(xpeak-size(J{1,eachcell},2)) (ypeak-size(J{1,eachcell},1))];

            if (frameNumber == 1)
                corr_offset_before = corr_offset;
                corr_offset_total_before = corr_offset_total{frameNumber};
            else
                normcell = norm(corr_offset-corr_offset_before);
                normframe = norm(corr_offset_total{frameNumber}-corr_offset_total_before);
                if (normcell>normframe+10)
                    corr_offset = corr_offset_before;
                else
                    corr_offset_before = corr_offset;
                    corr_offset_total_before = corr_offset_total{frameNumber};
                end
            end

            %% display best match

            figure(f1), imshow(K); hold on;
            set(f1, 'NumberTitle', 'off', ...
                'Name', sprintf('Frame: %s', int2str(frameNumber)));
            rectangle('position',[corr_offset(1) corr_offset(2) rect{1,eachcell}(3) rect{1,eachcell}(4)],...
                'edgecolor','g','linewidth',1);

            cell = imcrop(K,[corr_offset(1) corr_offset(2) rect{1,eachcell}(3) rect{1,eachcell}(4)]);

            %figure(f2), imshow(cell); hold on;
            %% center crop circle
            %if all we want is to probe for the circle
            center = floor((size(cell)/2+.5));
            %viscircles(center,5)


            bw = imbinarize(cell);
            bw = bwareaopen(bw,50);
            se = strel('disk',5);
            bw = imclose(bw,se);
            bw = imfill(bw,'holes');
            [B,L] = bwboundaries(bw,'noholes');
            cellbefore = cell;
        end
            
        for k = 1:length(B)
   boundary = B{k};
   %plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1)
end
            %boundary = B{1};
            %plot(boundary(:,2),boundary(:,1),'g','LineWidth',1);

            brilho=[pixelvalue.MeanIntensity];
            pixelVals=[pixelVals brilho(1)];
       
        end
   % V = 1:1:numofframes-1;
   
    %figure;
    maxval = max(pixelVals);
    pixelVals = double(pixelVals)/double(pixelVals(1));
    %plot(V,pixelVals)
   
    output{end+1}=pixelVals';
   
    pixelVals = [];
%     if(eachcell~=length(J))
%         m=input('Do you want to continue to next cell, Y/N [Y]:','s')
%         if m=='N'
%             break
%         end
%     end
end


for (i=1:length(output))
    for(j=1:length(output{1,i}))
        if (~isempty(output{1,i}))
            xlsoutput(j,i) = output{1,i}(j,1);
        else
            xlsoutput(j,i) = 0;
        end
    end
end
xlswrite(excelname,xlsoutput,'Folha1');
outputfile = [pwd '\' excelname];
RemoveSheet123(outputfile);
end

