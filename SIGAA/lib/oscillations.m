function oscillations(excelsheet,testname)
    close all
    file=excelsheet
    marge=2;    % limiar = average + marge*dv
    fc=1;       % 1 - compensates; 0 - do not compensates
    compensacao=2400;
    k=4;        % order of the compensation polynomial
    umaauma=1;
    baseline=300;  % 5 minutes
    ni=1200;
    nl=2400;
    duracao_min=10;   % above the margin lasting longer than
    duracao_max=100000;   % below the margin lasting longer than
    gravacao=5;      % sampling period
    noise=1;         % Clean noise
    normalizepeaks=10; % normalize the number of transients by the variable
    
    localdir = pwd;
    
    calcio=xlsread(file);
    t=calcio(:,1);
    calcio(:,1)=[];
    [maxi, kf]=max(t>baseline);   % time function array index
    [maxi, nii]=max(t>ni);
    [maxi, nll]=max(t>nl);
    [maxi, comp]=max(t>compensacao);
    nduracao_min=round(duracao_min/gravacao,0);   % peak sample number
    nduracao_max=round(duracao_max/gravacao,0);   % peak sample number
    [li, cl]=size(calcio);

    fncslopes = {};
    rises = {};
    decays = {};
    risestimes = {};
    decaystimes = {};
    betweenpeakstimes = {};
    peakstarttimes = {};
    peaktimes = {};
    
    for i=1:size(calcio)
         legendtext{i} = ['Celula ' num2str(i)];        
    end
    
    
    plot(t,calcio);
    xlabel('Time [s]') 
    ylabel('Ca2+ intensity [Ratio 340/380]')
    legend(legendtext);
    %hold on
    pause


    kf=kf-1;
    for i=1:cl  % cell by cell
        f1 = figure;
        title(['celula ' num2str(i)]);
        %      m=(calcio(comp,i)-calcio(1,i))/(comp-1)  % compensation slope
        %      declives(i)=m
        %      d=calcio(:,i)-fc*m*(1:li)';             % compensation
        if fc % compensation with k-order polynomial interpolation
            p=polyfit(t(1:comp), calcio(1:comp,i), k);
            %p=polyfit(t(nii:comp), calcio(nii:comp,i), k);
            d=calcio(:,i)-polyval(p, t)+polyval(p,t(1));
            bl = polyval(p, t);
        else
            d=calcio(:,i);
            p=polyfit(t(1:comp), calcio(1:comp,i), k);
            bl = polyval(p, t);
           % to plot ca2+ function versus the baseline
            plot(t,d);
            hold on
%             plot(t,bl,'r');
            xlim([0 3000])
%            ylim([0.9 2.8])
            xlabel('Time [s]')
            ylabel('Ca2+ intensity [Ratio 340/380]')
            legendtext{1} = ['Ca2+ function'];
            legendtext{2} = ['Baseline degree ' num2str(k)];
            legend(legendtext);

             pause
            %bl = polyval(p, t)+calcio(1,i)-polyval(p,t(1));
        end
        if noise   % MV noise reduction - moving average
            e=d;
            for j=3:length(d)-2
                e(j)=mean(d(j-2:j+2));
            end
            d=e;
        end
% to plot noise and compensation effect        
% plot(t,d);
% hold on
% %plot(t,bl,'r');
% xlabel('Time [s]') 
% ylabel('Ca2+ intensity [Ratio 340/380]')
% xlim([0 3000])
% legendtext{1} = ['Ca2+ function']
% legendtext{2} = ['Baseline degree ' num2str(k)]
% legend(legendtext)
% pause

        medcel(1,i)=mean(d(1:kf));
        medcel(2,i)=std(d(1:kf));
        medcel(3,i)=medcel(1,i)+marge*medcel(2,i);

        threshold=medcel(3, i);
        plot(t,d,'DisplayName',['Cell ' num2str(i)]);
        xlabel('Time [s]') 
        ylabel('Ca2+ intensity [Ratio 340/380]')
        xlim([0 3000])
        %ylim([-1 4])
        hold on
        %pause
        x(1:li)=threshold;
        plot(t,x,'r','DisplayName','Transient threshold');


        pico=d>threshold;   % peaks above threshold
        pico(1:nii-1)=0;   % clean below nii
        pico(nll:end)=0;   % clean above nll
        picon=pico;

        ptr=nii;
        j=2;
        while(ptr<nll)
            [maxi, niii]=max(pico(ptr:nll));
            niii=niii+ptr-1;
            [mini, nlll]=min(pico(niii:nll));
            nlll=nlll+niii-2;
            if(nlll<niii)
                if(j==2)
                    amp(:,i)=0;
                    fprintf('Warning: Célula %d sem picos!!!!!!', i);
                end
                break;
            end
            %        if( (nlll-niii+1>nduracao_min) && (nlll-niii+1<nduracao_max) && (nii~=niii) )
            if( (nlll-niii+1>nduracao_min) && (nlll-niii+1<nduracao_max)) % && (nii~=niii) )
                nlll-niii+1
                [amp(j,i), ind(j,i)]=max(calcio(niii: nlll,i));
                j=j+1;
                picon(niii:nlll)=0;
            else
                pico(niii:nlll)=0;
            end
            ptr=nlll+1;
            if(normalizepeaks>0)
                amp(1,i)=(j-2)./normalizepeaks;
            else
                amp(1,i)=j-2;
            end
        end


        %% Estimates Peaks Features

        piconotzero = find(pico); % "1" indexes


        for (ipico=1:length(piconotzero))
            for (z=piconotzero(ipico):-1:2)
                %if it is the first index of the signal, it jumps%
                %if the previous one is already a peak, jump%
                if(pico(z-1)== 1)
                    break
                end
                % if two indices before is a peak, skip to avoid peak union
                if (z > 2 &&  pico(z-2)== 1)
                    break
                end
                 % if the previous value of the signal is less than the current value, 
                 % then this index belongs to this peak
                if(d(z)>d(z-1))
                    pico(z-1)=1;
                else
                    break
                end
            end

            for (z=piconotzero(ipico):(length(pico)-1))
                %if it is the last index signal index jump if the next one is already a peak jump%
                if(pico(z+1)== 1)
                    break
                end
                % if two indices later is a peak, skip to avoid peak union
                if (z < (piconotzero(end)-1) &&  pico(z+2)== 1)
                    break
                end
                % if the next signal value is less than the current value, 
                % then this index belongs to this peak
                if(d(z)>d(z+1))
                    pico(z+1)=1;
                else
                    break
                end
            end
        end

        % for debug %
        piconotzerodepois = find(pico);

        % collect features

        peaks = pico .* d;

        fncslopes{end+1}= getslopesfivemin(calcio(:,i),t);

        [dividedpeaks,dividedtimes] = getpeaks(peaks,t);
        [rises{end+1},decays{end+1},risestimes{end+1},decaystimes{end+1},betweenpeakstimes{end+1},peakstarttimes{end+1},peaktimes{end+1}] = getpeaksfeatures(dividedpeaks, dividedtimes, threshold, 'linearreg');

        

        %%

        %plot(t,1.3*threshold*picon,'g')
        plot(t,1.3*threshold*pico,'DisplayName','Transient')
        %plot([ni ni], [0 max(d)], 'k','DisplayName','Analysis starting point')
        %plot([nl nl], [0 max(d)], 'k')        
        plot(t,bl,'DisplayName','Baseline');        
        legend;

        if(umaauma)
            fprintf('Célula %d\n', i);
            disp(amp(:,i));
            if(i~=cl)
                jj=input('Célula seguinte...', 's');
            end
        end
    end

    if(cl>0)
        medcel;
        amp;
        x=size(amp);
        %surfc(1:x(2), 1:x(1), amp, 'Edgecolor', 'none');
        %axis xy; axis tight; colormap(jet); view(0,90);

        %surf(1:cl, t, calcio)
        %axis xy; axis tight; colormap(jet); %view(0,90);
%         figure
%         plot(t,calcio)
%         %set(gca, 'XLim', [0, 3000], 'XTick', 0:100:3000,...
%    % 'XTickLabel', 0:100:3500);
%     xlabel('Time [s]') 
%     ylabel('Ca2+ intensity [Ratio 340/380]')

        ylim([-2 2])
        axis tight
    end

    % daqui
    
    %% draw arrows for peaktime analysis
    f2 = figure;
    title('Comunicação inter-celular');
    orderedpeaktimes = {};
    
    for (i=1:length(peakstarttimes))
        if (isempty(peakstarttimes{1,i}))
            continue
        end
        orderedpeaktimes{end+1}=[i peakstarttimes{1,i}{1,1}];
    end
    %sort
    for (i=1:length(orderedpeaktimes))
        for (j=i+1:length(orderedpeaktimes))
            if (orderedpeaktimes{i}(1,2)>orderedpeaktimes{j}(1,2))
                aux = orderedpeaktimes{i};
                orderedpeaktimes{i} = orderedpeaktimes{j};
                orderedpeaktimes{j} = aux;
            end
        end
    end
%    opts = detectImportOptions([localdir '\outputs\' testname '\ROIs.xlsx'],'NumHeaderLines',0);
    opts = detectImportOptions([testname '\ROIs.xlsx'],'NumHeaderLines',0);
    T= readtable([testname '\ROIs.xlsx'],opts,'ReadVariableNames',true);
    refframename=dir([testname, '\*_refframe.bmp']);
    propagationtimes = imread([testname '\' refframename(1).name]);
    imshow(propagationtimes); drawnow;hold on;
    for (i=1:length(orderedpeaktimes)-1)
    
        p1 = [T(orderedpeaktimes{1,i}(1),1) T(orderedpeaktimes{1,i}(1),2)];                         % First Point
        p1 = table2array(p1);
        p2 = [T(orderedpeaktimes{1,i+1}(1),1) T(orderedpeaktimes{1,i+1}(1),2)];                          % Second Point
        p2 = table2array(p2);

        dp = p2-p1;                         % Difference
               
        quiver(p1(1),p1(2),dp(1),dp(2),0,'color','g')
        grid
        text(p1(1),p1(2), sprintf('%.0f',i),'Color', 'g')
        %
        %text(p2(1),p2(2), sprintf('(%.0f,%.0f)',p2),'Color', 'g')
        if(i==length(orderedpeaktimes)-1) 
            text(p2(1),p2(2), sprintf('%.0f',i+1),'Color', 'g')
        end
    end
    
    saveas(gcf, [testname '\propagationtimes.bmp'])

    %here
  %% save results to disk
    s=input('Do you want to save?','s')
    if(s(1)~='n' || length(s)~=1)
        amp(1, cl+1)=marge;
        amp(2, cl+1)=fc;
        amp(3, cl+1)=300;
        amp(4, cl+1)=duracao_min;
        amp(5, cl+1)=duracao_max;
        amp(6, cl+1)=noise;
        %currentFolder = pwd
        %a=[testname '\results.xlsx'];
        xlswrite([testname '\results.xlsx'], amp, 'Peak amplitude');
        
        amp_column = amp;
        amp_column([1],:) = [];
        amp_column=amp_column(:,1:end-1);
        amp_column = amp_column(:);
        amp_column = nonzeros(amp_column);
        xlswrite([testname '\results.xlsx'], amp_column, 'Non zero amplitudes');
        
        g=1;
        for (i=1:length(rises))
            outputrisedecay(1,g)=length(rises{1,i});
            outputrisedecay(1,g+1)=length(decays{1,i});
            k=2;
            for(j=1:length(rises{1,i}))
                if(j~=1)
                    k=k+1;
                end

                if (~isempty(rises{1,i}))
                    outputrisedecay(k,g)=rises{1,i}{1,j};
                else
                    outputrisedecay(k,g) = 0;
                end
            end
            k=2;
            for(j=1:length(decays{1,i}))
                if(j~=1)
                    k=k+1;
                end

                if (~isempty(decays{1,i}))
                    outputrisedecay(k,g+1)=decays{1,i}{1,j};
                else
                    outputrisedecay(k,g+1) = 0;
                end
            end
            g=g+2;
        end
        xlswrite([testname '\results.xlsx'],outputrisedecay,'Slope (rise e decay)');

        for (i=1:length(risestimes))
            outputrisedecaytimes(1,i)=length(risestimes{1,i});
            k=2;
            for(j=1:length(risestimes{1,i}))
                if(j~=1)
                    k=k+1;
                end

                if (~isempty(risestimes{1,i}))
                    outputrisedecaytimes(k,i)=risestimes{1,i}{1,j};
                    k=k+1;
                    outputrisedecaytimes(k,i)=decaystimes{1,i}{1,j};
                else
                    outputrisedecaytimes(k,i) = 0;
                end
            end
        end
        xlswrite([testname '\results.xlsx'],outputrisedecaytimes,'Rise and decay time (s)');


        for (i=1:length(risestimes))
            if(length(betweenpeakstimes{1,i})>0)
                for(j=1:length(betweenpeakstimes{1,i}))
                    if (~isempty(betweenpeakstimes{1,i}))
                        bptoutput(j,i) = betweenpeakstimes{1,i}{1,j};
                    else
                        bptoutput(j,i) = 0;
                    end
                end
            else
                bptoutput(1,i) = 0;
            end
        end
        xlswrite([testname '\results.xlsx'],bptoutput,'Time between peaks (s)');

        for (i=1:length(risestimes))
            if(length(peakstarttimes{1,i})>0)
                for(j=1:length(peakstarttimes{1,i}))
                    if (~isempty(peakstarttimes{1,i}))
                        bptoutput(j,i) = peakstarttimes{1,i}{1,j};
                    else
                        bptoutput(j,i) = 0;
                    end
                end
            else
                bptoutput(1,i) = 0;
            end
        end
        xlswrite([testname '\results.xlsx'],bptoutput,'Start peak time (s)');
        
        for (i=1:length(risestimes))
            if(length(peaktimes{1,i})>0)
                for(j=1:length(peaktimes{1,i}))
                    if (~isempty(peaktimes{1,i}))
                        bptoutput(j,i) = peaktimes{1,i}{1,j};
                    else
                        bptoutput(j,i) = 0;
                    end
                end
            else
                bptoutput(1,i) = 0;
            end
        end
        xlswrite([testname '\results.xlsx'],bptoutput,'Peak time (s)');
        
        for (i=1:length(fncslopes))
            if(length(fncslopes{1,i})>0)
                for(j=1:length(fncslopes{1,i}))
                    if (~isempty(fncslopes{1,i}))
                        fncsoutput(j,i) = fncslopes{1,i}{1,j};
                    else
                        fncsoutput(j,i) = 0;
                    end
                end
            else
                fncsoutput(1,i) = 0;
            end
        end
        xlswrite([testname '\results.xlsx'],fncsoutput,'5 minutes slope');


    end
end
