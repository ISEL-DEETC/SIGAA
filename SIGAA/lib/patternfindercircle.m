function [excelname testname] = patternfindercircle(localdir)
%%Read input frame and take the ROI

batchname = input('What batch to work? ','s');
batchname = ['\inputs\' batchname '\'];

locdir = strcat(localdir,batchname);

name=input('Name the output test name? ','s');
currDate = strrep(datestr(datetime), ':', '_');
testname = [pwd '\outputs\' name '_' currDate];
mkdir(testname)
excelname = [testname '\ca2Function.xlsx'];


%numofframes = size(dir([locdir '*.bmp']),1);

inffiles = dir([locdir '*.INF']);
timename=inffiles(1).name;

time = gettime(strcat(locdir,timename));
clear timename

%% create process ratio frames

ratiodir = [locdir 'ratio\'];
dir340 = [locdir '340\'];
dir380 = [locdir '380\'];
numofframes = size(dir([dir340 '*.bmp']),1);

bmp340files = dir(strcat(dir340, '*.bmp'));
sorted340bmpfiles = natsortfiles({bmp340files.name});
bmp380files = dir(strcat(dir380, '*.bmp'));
sorted380bmpfiles = natsortfiles({bmp380files.name});



tic
if ((~exist(ratiodir, 'dir')) | ((exist(ratiodir, 'dir')) && (length(dir([ratiodir, '\*.mat']))<length(dir([dir380, '\*.bmp'])))))
    if ~exist(ratiodir, 'dir')
        mkdir(ratiodir)
    end
    bar = waitbar(0,'Please wait...');
    for frameNumber=1:numofframes
        waitbar(frameNumber/numofframes,bar,'Processing 340/380 ratio video frames.');
        
        %frame340 = imread(strcat(dir340, sorted340bmpfiles{frameNumber}));
        
%          f1 = figure;
%          figure(f1), imshow(frame340)
         
         frame340 = rgb2gray(imread(strcat(dir340, sorted340bmpfiles{frameNumber})));
%          x = max(max (frame340bw))
%          y = min(min (frame340bw))
%          f1 = figure;
%          figure(f1), imshow(frame340)
%          teste = 255-frame340bw;
%          f4 = figure;
%          figure(f4), imshow(teste)
         
         
         
        frame380 = rgb2gray(imread(strcat(dir380, sorted380bmpfiles{frameNumber})));
% 
%         f2 = figure;
%          figure(f2), imshow(frame380)
        [r,c] = size(frame380);
        
        for i=1:r
            for j=1:c
                if frame380(i, j)>0
                    ratioframe(i, j)=(double(frame340(i,j)))./(double(frame380(i,j)));
                else
                    ratioframe(i, j)=0;
                end
            end
        end
        %ratioframe = double(frame340)./double(frame380);

%         x = min (ratioframe)
%         y = max (ratioframe)
%         ratioframe = ratioframe.*255;
%         ratioframe= uint8(ratioframe);
%          f3 = figure;
%          figure(f3), imshow(ratioframe)       

                if ~exist(ratiodir, 'dir')
                    mkdir(ratiodir)
                end

        save(fullfile(ratiodir, ['ptx' int2str(frameNumber) '.mat']),'ratioframe');
        
    end
    delete(bar);
end

%% create video

bmpfiles = dir(strcat(dir380, '*.bmp'));
sortedbmpfiles = natsortfiles({bmpfiles.name});
cellAvi = VideoWriter(fullfile(dir380,'cellAvi.avi'));
open(cellAvi)


for frameNumber = 1:numofframes
    
   img = im2double(imread(strcat(dir380, sortedbmpfiles{frameNumber})));
   writeVideo(cellAvi,img)
end


close(cellAvi)
clear cellAvi img inffiles

%% play video
vid = VideoReader(fullfile(dir380,'cellAvi.avi'));

% Set up video figure window
videofig(vid.NumberOfFrames, @(frm) redraw(frm, vid));

% Display initial frame
redraw(1, vid);

pause
clear vid

J={};
rect={};
output={};
output{end+1}=time';
ref_panel=[];

rframe= im2double(imread(strcat(dir380, sortedbmpfiles{str2num('1')})));
I = rframe;
[r c] = size(I);
[J{end+1}, rect{end+1}] = imcrop(I,[0 0 c r]);

I = rgb2gray(I);
%f2 = figure;

%% Processing Video to Generate Significant region

processeddir = [dir380 'processed\'];
processeddir_ratio = [locdir 'processed_ratio\'];
matratiofiles = dir(strcat(ratiodir, '*.mat'));
sortedmatratiofiles = natsortfiles({matratiofiles.name});

tic
%f1=figure;
if ((~exist(processeddir, 'dir')) | ((exist(processeddir, 'dir')) && (length(dir([processeddir, '\*.bmp']))<length(dir([dir380, '\*.bmp'])))))
    %fprintf('Finding Significant Region\n')
    
    bar = waitbar(0,'Please wait...');
    for iteration=1:2
       
        for frameNumber=1:numofframes
            
            % get frame number frameNumber
            if (iteration==1)
              waitbar(frameNumber/numofframes,bar,'Finding Significant Region.(step 1/2)');
            else
             waitbar(frameNumber/numofframes,bar,'Aligning the Video. (step 2/2)');   
            end
            
            if (iteration==2&&frameNumber==1)
                [r c] = size(ref_panel);
                corr_offset_total={};
                rect{1,1}(3)=c;
                rect{1,1}(4)=r;
            end
            K_color = im2double(imread(strcat(dir380, sortedbmpfiles{frameNumber})));
            K = rgb2gray(K_color);
            
            ratio_frame = load(strcat(ratiodir, sortedmatratiofiles{frameNumber}));

            if(iteration==1)
                t = normxcorr2mem(I, K);
            else
                t = normxcorr2mem(ref_panel, K);
            end
            clearAllMemoizedCaches
            [max_t, imax] = max(abs(t(:)));
            %Offset between the images found by correlation
            [ytpeak, xtpeak] = ind2sub(size(t),imax(1));
            if(iteration==1)
                corr_offset_total{frameNumber} = [(xtpeak-size(I,2)) (ytpeak-size(I,1))];
            else
                corr_offset_total{frameNumber} = [(xtpeak-size(ref_panel,2)) (ytpeak-size(ref_panel,1))];
            end
            clear max_t imax ytpeak xtpeak t
            
            %% display best match
            
%             figure(f1), imshow(K_color); drawnow;hold on;
%             set(f1, 'NumberTitle', 'off', ...
%                'Name', sprintf('Frame: %s', int2str(frameNumber)));
%             rectangle('position',[corr_offset_total{1,frameNumber}(1) corr_offset_total{1,frameNumber}(2) rect{1,1}(3) rect{1,1}(4)],...
%                'edgecolor','g','linewidth',1);drawnow;
            
            cell = imcrop(K,[corr_offset_total{1,frameNumber}(1) corr_offset_total{1,frameNumber}(2) rect{1,1}(3) rect{1,1}(4)]);
            bw_cell = cell;
            
            cell2 = imcrop(ratio_frame.ratioframe,[corr_offset_total{1,frameNumber}(1) corr_offset_total{1,frameNumber}(2) rect{1,1}(3) rect{1,1}(4)]);
                        
            if(iteration==1)
                [rows_bw, columns_bw] = size(bw_cell);
                [rows_ref, columns_ref] = size(ref_panel);
                if(frameNumber==1)
                    ref_panel = bw_cell;
                else
                    ref_panel = imcrop(bw_cell,[0 0 min(columns_bw,columns_ref) min(rows_bw,rows_ref)]);
                end
%                 clf(f2,'reset')
%                 figure(f2), imshow(ref_panel); drawnow;hold on;
            else
                if ~exist(processeddir, 'dir')
                    mkdir(processeddir)
                end
                if ~exist(processeddir_ratio, 'dir')
                    mkdir(processeddir_ratio)
                end
                %% save new video
                imwrite(cell,fullfile(processeddir, sortedbmpfiles{frameNumber}));
                save(fullfile(processeddir_ratio,sortedmatratiofiles{frameNumber}), 'cell2');
            end
        end
    end
    delete(bar);
    
%     %% resize all processed pictures
%     bar = waitbar(0,'Please wait...');
%     min_row=0;
%     min_column=0;
%     for frameNumber=1:numofframes
%         waitbar(frameNumber/(numofframes),bar,'Resizing aligned video frames. 1/2');
%         K = load(strcat(processeddir_ratio, sortedmatratiofiles{frameNumber}));
%         
%         [rows, columns, numberOfColorChannels] = size(K.cell2);
%         
%         if (frameNumber==1)
%             min_row=rows;
%             min_column=columns;
%         else
%             if (rows<min_row)
%                 min_row=rows;
%             end
%             if (columns<min_column)
%                 min_column=columns;
%             end
%         end
%     end
%     delete(bar);
%     bar = waitbar(0,'Please wait...');
%     for frameNumber=1:numofframes
%         waitbar(frameNumber/(numofframes),bar,'Resizing aligned video frames. 2/2');
%         K = load(strcat(processeddir_ratio, sortedmatratiofiles{frameNumber}));
%         cell2=imcrop(K.cell2,[0 0 min_column min_row]);
%         save(fullfile(processeddir_ratio,sortedmatratiofiles{frameNumber}), 'cell2');
%         
%     end
%     delete(bar);
%     
end
toc


%% This section to check for cells

%rframe= im2double(imread(strcat(processeddir, sortedbmpfiles{str2num('1')})));

I = imread(strcat(processeddir, sortedbmpfiles{str2num('1')}));
%imshow(I)
%I=rgb2gray(I_color);



% % gmag = imgradient(I);
%  imshow(I)
%  title('Gradient Magnitude')
%  pause;
% I2 = imcomplement(I);
% I3 = imhmin(I2,5); %20 is the height threshold for suppressing shallow minima
% L = watershed(I3);
%  Lrgb = label2rgb(L);
%  imshow(Lrgb)
%  title('Transformada Watershed')
%  pause;

se = strel('disk',7);
% Io = imopen(I,se);
% imshow(Io)
% title('Opening')
% pause;

Ie = imerode(I,se);
Iobr = imreconstruct(Ie,I);
%    imshow(Iobr)
%    title('Opening-by-Reconstruction')
%    pause;

% Ioc = imclose(Io,se);
% imshow(Ioc)
% title('Opening-Closing')
% pause;

Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
%   imshow(Iobrcbr)
%   title('Opening-Closing by Reconstruction')
%   pause;

fgm = imregionalmax(Iobrcbr);
%   imshow(fgm)
%   title('Regional Maxima of Opening-Closing by Reconstruction')
%   pause;

%I2 = labeloverlay(I,fgm);
%  imshow(I2)
%  title('Regional Maxima Superimposed on Original Image')
%  pause;

se2 = strel(ones(5,5));
fgm2 = imclose(fgm,se2);
fgm3 = imerode(fgm2,se2);


fgm4 = bwareaopen(fgm3,20);

I3 = labeloverlay(I,fgm4);
%   imshow(fgm4)
%   title('Modified Regional Maxima Superimposed on Original Image')
%   
%   pause;


[centers,radii]=imfindcircles(fgm4,[5,30],'ObjectPolarity','bright','Sensitivity',0.85);
fprintf('Found %d cells.\n',size(centers,1));
figure(1); imshow(I); hold on;


for idx=1:length(centers)
    
    %rectangle('position',[centers(idx,1)-2*radii(idx) centers(idx,2)-2*radii(idx) 4*radii(idx) 4*radii(idx)],...
    %           'edgecolor','g','linewidth',1);drawnow;
    
    cell = imcrop(fgm4,[centers(idx,1)-2*radii(idx) centers(idx,2)-2*radii(idx) 4*radii(idx) 4*radii(idx)]);
    bw = bwareaopen(cell,50);
    se = strel('disk',5);
    bw = imclose(bw,se);
    bw = imfill(bw,'holes');
    [B,L] = bwboundaries(bw,'noholes');    
    boundary = B{1};
    plot(centers(idx,1)-2*radii(idx)+boundary(:,2), centers(idx,2)-2*radii(idx)+boundary(:,1), 'g', 'LineWidth', 1)    
    text(centers(idx,1)-2*radii(idx),centers(idx,2)-2*radii(idx),num2str(idx),'Color','green','FontSize',14);drawnow;
    
end


%viscircles(centers, radii,'Color','g');
%% STOP

%pause(5);

%%ASK for the cells to study
prompt='Please specify if there are cells you wish to remove? \n ("all" to remove all the cells, "no" to keep all the cells, or  the cellnumber you want to remove. "done" when done. )';
str = input(prompt,'s');

while not(strcmp(str,'done') | strcmp(str,'no'))
    if(strcmp(str,'all'))
        radii = [];
        centers = [];        
        break;
    end    
    radii(str2double(str)) = 0; 
    centers(str2double(str),:) = 0;    
    str = input(prompt,'s') ;
end

%radii = radii(radii ~= 0);
%centers = centers(centers ~= 0);

close all

figure(1); imshow(I); hold on;

centerRoi = {};
removedidx = 0;
for idx=1:length(centers)
    
    %rectangle('position',[centers(idx,1)-2*radii(idx) centers(idx,2)-2*radii(idx) 4*radii(idx) 4*radii(idx)],...
    %    'edgecolor','g','linewidth',1);drawnow;
    %text(centers(idx,1)-2*radii(idx)+5,centers(idx,2)-2*radii(idx)+15,num2str(idx),'Color','green','FontSize',14);drawnow;
    if  centers(idx)~=0
        cell = imcrop(fgm4,[centers(idx,1)-2*radii(idx) centers(idx,2)-2*radii(idx) 4*radii(idx) 4*radii(idx)]);
        bw = bwareaopen(cell,50);
        se = strel('disk',5);
        bw = imclose(bw,se);
        bw = imfill(bw,'holes');
        [B,L] = bwboundaries(bw,'noholes');
        boundary = B{1};
        plot(centers(idx,1)-2*radii(idx)+boundary(:,2), centers(idx,2)-2*radii(idx)+boundary(:,1), 'g', 'LineWidth', 1)
        text(centers(idx,1)-2*radii(idx),centers(idx,2)-2*radii(idx),num2str(idx-removedidx),'Color','green','FontSize',14);drawnow;
        
        centerRoi{end+1} = [centers(idx,1) centers(idx,2)];
    else
        removedidx = removedidx + 1;
    end
end


if isempty(centers)
    automaticcells = 0;
else
    automaticcells = length(centers(:,1));
end


manualcells = uint8(input('How many Cells to add? '));

J = {};
rect = {};
manual=0;
quit=0;
for crop=1:manualcells
    if (crop>1)
        L=[];
        while (nnz(L) == 0)
            %rectangle('position',[rect{end}(1) rect{end}(2) rect{end}(3) rect{end}(4)],...
            %    'edgecolor','g','linewidth',1);drawnow;
            cell = imcrop(fgm4,[rect{end}(1) rect{end}(2) rect{end}(3) rect{end}(4)]);
            bw = bwareaopen(cell,50);
            se = strel('disk',5);
            bw = imclose(bw,se);
            bw = imfill(bw,'holes');
            [B,L] = bwboundaries(bw,'noholes');
            if(nnz(L) == 0)
                %acrescentei isto
                [nx,ny,d] = size(~cell);
                [X,Y] = meshgrid(1:ny,1:nx) ;
                px = round(size(cell)/2+.5);
                ang=0:0.01:2*pi;
                xp=5*cos(ang);
                yp=5*sin(ang);
                plot(rect{end}(1)+px(1)+xp,rect{end}(2)+px(2)+yp,'g');drawnow;
                x = rect{end}(1)+px(1)+xp;
                y = rect{end}(2)+px(2)+yp;
                boundary = [x(:)', y(:)'];
                L=1;
            else                
                boundary = B{1};
                plot(rect{end}(1)+boundary(:,2), rect{end}(2)+boundary(:,1), 'g', 'LineWidth', 1);drawnow;
            end
            manual = manual +1;            
            if automaticcells == 0
                text(rect{1,crop-1}(1),rect{1,crop-1}(2),num2str(crop-1-removedidx),'Color','green','FontSize',14);drawnow;
            else
                text(rect{1,crop-1}(1),rect{1,crop-1}(2),num2str(automaticcells+crop-1-removedidx),'Color','green','FontSize',14);drawnow;
            end
        end
    end
    if(quit==1)
        manualcells = manual;
        break;
    end
    [J{end+1}, rect{end+1}] = imcrop(gcf);
    %lets save the test conditions
    if (crop==manualcells)        
        %rectangle('position',[rect{end}(1) rect{end}(2) rect{end}(3) rect{end}(4)],...
        %    'edgecolor','g','linewidth',1);drawnow;
        L=[];
        while (nnz(L) == 0)
            cell = imcrop(fgm4,[rect{end}(1) rect{end}(2) rect{end}(3) rect{end}(4)]);
            bw = bwareaopen(cell,50);
            se = strel('disk',5);
            bw = imclose(bw,se);
            bw = imfill(bw,'holes');
            [B,L] = bwboundaries(bw,'noholes');
            if(nnz(L) == 0)
                %acrescentei isto
                [nx,ny,d] = size(~cell);
                [X,Y] = meshgrid(1:ny,1:nx) ;
                px = round(size(cell)/2+.5);
                ang=0:0.01:2*pi;
                xp=5*cos(ang);
                yp=5*sin(ang);
                plot(rect{end}(1)+px(1)+xp,rect{end}(2)+px(2)+yp,'g');drawnow;
                x = rect{end}(1)+px(1)+xp;
                y = rect{end}(2)+px(2)+yp;
                boundary = [x', y'];
                L=1;
                %boundary = boundary(rect{end}(1)+px(1)+xp,);
            else                
                boundary = B{1};
                plot(rect{end}(1)+boundary(:,2), rect{end}(2)+boundary(:,1), 'g', 'LineWidth', 1);drawnow;
            end
            manual = manual +1;            
            if automaticcells == 0
                text(rect{1,crop}(1),rect{1,crop}(2),num2str(crop),'Color','green','FontSize',14);drawnow;
            else
                text(rect{1,crop}(1),rect{1,crop}(2),num2str(automaticcells+crop-removedidx),'Color','green','FontSize',14);drawnow;
            end
        end
    end
end

saveas(gcf, [testname '\' name '_startcondition.bmp'])
refframe = im2double(getimage(gcf));
imwrite(refframe, [testname '\' name '_refframe.bmp']);

for i=1:manualcells
    xCentroid = rect{i}(1) + rect{i}(3)/2;
    yCentroid = rect{i}(2) + rect{i}(4)/2;
    centerRoi{end+1} = [xCentroid yCentroid];
end
xlswrite([testname '\ROIs.xlsx'],cell2mat(centerRoi'));

%%
%Create ca2+ the functions
pixelVals = [];
tic
cellindicator = 0;
for (cellnum=1:automaticcells+manualcells)
    % se for uma célula removida, salta
    
    if(cellnum<automaticcells&&centers(cellnum,1) == 0)
        continue;
    else
       cellindicator = cellindicator+1;
    end
    
    bar = waitbar(0,'Please wait...');
    pri = sprintf('Processing cell %d, please wait.',cellindicator);
%     f10 = figure;
%     f11 = figure;
    for frameNumber=1:numofframes
        waitbar(frameNumber/numofframes,bar,pri);
        %collect rgb greyscale
        K= load(strcat(processeddir_ratio, sortedmatratiofiles{frameNumber}));
%         figure(f11), imshow(K.cell2),hold on
%         rectangle('position',[rect{cellnum-automaticcells}(1) rect{cellnum-automaticcells}(2) rect{cellnum-automaticcells}(3) rect{cellnum-automaticcells}(4)],...
%          'edgecolor','g','linewidth',1);drawnow;
        %collect luminance
        %auxiliar=imread(strcat(processeddir, sortedbmpfiles{frameNumber}));
        %(0.2126*R + 0.7152*G + 0.0722*B)
        %K=im2double(auxiliar(:, :, 1).*0.2126+im2double(auxiliar(:, :, 2)).*0.7152+im2double(auxiliar(:, :, 3)).*0.0722);
        %K=im2double(auxiliar(:, :, 1)).*0.7+im2double(auxiliar(:, :, 3)).*0.3;

        if (manualcells == 0 | cellnum<automaticcells+1)
            
            % processar brilho celula automatica
            
            cell = imcrop(K.cell2,[centers(cellnum,1)-2*radii(cellnum) centers(cellnum,2)-2*radii(cellnum) 4*radii(cellnum) 4*radii(cellnum)]);

        else
            % processar brilho celula manual
            cell = imcrop(K.cell2,[rect{cellnum-automaticcells}(1) rect{cellnum-automaticcells}(2) rect{cellnum-automaticcells}(3) rect{cellnum-automaticcells}(4)]);
                   
%           figure(f10), imshow(cell),hold on

        end
        if (frameNumber==1)
            
            [nx,ny,d] = size(~cell);
            [X,Y] = meshgrid(1:ny,1:nx) ;
            px = round(size(cell)/2+.5);
            ang=0:0.01:2*pi;
            xp=10*cos(ang);
            yp=10*sin(ang);
            
            [cellr,cellc] = size(cell);
        end
%         plot(px(1)+xp,px(2)+yp,'g');drawnow;
        [ro,co] = size(cell);
        if (ro==cellr && co==cellc)
            idx = inpolygon(X(:),Y(:),(px(1)+xp)',(px(2)+yp));
            brilhotmp =  mean2(cell(idx));
        end

        pixelVals=[pixelVals brilhotmp];
%         clf(f10,'reset')
%         clf(f11,'reset')
    end
    
    pixelVals = double(pixelVals)/double(pixelVals(1));
    output{end+1}=pixelVals';
    clear pixelVals
    pixelVals = [];
    delete(bar);
    toc
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
%outputfile = [testname '\' excelname];
RemoveSheet123(excelname);
end


