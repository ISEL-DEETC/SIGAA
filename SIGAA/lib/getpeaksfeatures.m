function [rises,decays,risestimes,decaystimes,betweenpeakstimes,peakstarttimes,peaktimes] = getpeaksfeatures(dividedpeaks,dividedtimes,threshold,type)
% This function receives two cell arrays, one with the amplitudes of the peaks 
% and the other with the times of each peak.
rises = {};
decays = {};
risestimes = {};
decaystimes = {};
betweenpeakstimes = {};
peakstarttimes = {};
peaktimes = {};

peakstarttime = 0;
peakendtime = 0;

% if the peaks start too fast, too far from the threshold value can harm the slope measurement, 
% in this way, we use the threshold value as the first and last value of the slope estimation

for i=1:length(dividedpeaks)
    
    %plot(dividedtimes{1,i},dividedpeaks{1,i})
    [maxvalue,index] = max (dividedpeaks{1,i}); % tirar maximo do pico %
    

    %%    
    % Estimate raise and decay times %
    % raise time = max - first peak indice %    
    risestimes{end+1} = dividedtimes{1,i}(index)-dividedtimes{1,i}(1);
    
    % decay time = last peak indice - max
    dividedtimes{1,i}(length(dividedtimes{1,i}))
    dividedtimes{1,i}(index)
    decaystimes{end+1} = dividedtimes{1,i}(length(dividedtimes{1,i}))-dividedtimes{1,i}(index);
    
    %% 
    % Estimate time between peaks %
    if i>1        
        betweenpeakstimes{end+1} = dividedtimes{1,i}(1) - peakendtime;        
    end
    
    % new beginning and end of peak %
    peakstarttime = dividedtimes{1,i}(1);
    peakendtime = dividedtimes{1,i}(length(dividedtimes{1,i}));

    peakstarttimes{end+1} = peakstarttime;
    peaktimes{end+1} = dividedtimes{1,i}(index);
    
    %% 
    % Slope deltapeaks / deltatimes
    % estimate slopes for each peak
    if strcmp(type, 'deltaslope')
        if ((dividedtimes{1,i}(index)-dividedtimes{1,i}(1))==0) % se não houver rise ou seja, curva só desce
            rises{end+1} = 0;
        else
            rises{end+1} = (dividedpeaks{1,i}(index) - threshold)./(dividedtimes{1,i}(index)-dividedtimes{1,i}(1));
            if ((dividedtimes{1,i}(length(dividedtimes{1,i}))-dividedtimes{1,i}(index))==0) % se não houver decay ou seja, curva só sobe
                decays{end+1} = 0;
            else
                decays{end+1} = (threshold - dividedpeaks{1,i}(index))./(dividedtimes{1,i}(length(dividedtimes{1,i}))-dividedtimes{1,i}(index));
            end
        end
    end
    % slope meanslope
    % estimate all average slopes on the ascent and descent
    if strcmp(type, 'meanslope')
        slopes = diff (dividedpeaks{1,i})./diff(dividedtimes{1,i})
        if ((dividedtimes{1,i}(index)-dividedtimes{1,i}(1))==0) % se não houver rise ou seja, curva só desce
            rises{end+1} = 0;
        else
            rises{end+1} = mean (slopes(1:index-1));
            if ((dividedtimes{1,i}(length(dividedtimes{1,i}))-dividedtimes{1,i}(index))==0) % se não houver decay ou seja, curva só sobe
                decays{end+1} = 0;
            else
                decays{end+1} = mean (slopes(index:length(slopes)));
            end
        end
    end
    
    if strcmp(type, 'linearreg')      
        
        r=polyfit(dividedtimes{1,i}(1:index),dividedpeaks{1,i}(1:index),1);
        rises{end+1} = r(1);
        d=polyfit(dividedtimes{1,i}(index:length(dividedtimes{1,i})),dividedpeaks{1,i}(index:length(dividedpeaks{1,i})),1);
        decays{end+1} = d(1);
    end
    
end

end