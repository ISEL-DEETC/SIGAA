function [dividedpeaks,dividedtimes] = getpeaks(peaks,times)
    idxpeaks = find(peaks); % peak indices
    dividedpeaks = {};
    dividedtimes = {};
    newpeak = 1;
    
    for k=1:length(idxpeaks)
        if newpeak == 1 % if first index creates a new peak
            peak = peaks(idxpeaks(k));  
            time = times(idxpeaks(k));
            newpeak = 0;
        else                
            if idxpeaks(k)==(idxpeaks(k-1)+1) % if indices are followed, both belong to the same peak
                peak = [peak peaks(idxpeaks(k))];
                time = [time times(idxpeaks(k))];
                if(k==length(idxpeaks)) % ends flow, closes peak and adds to cell array with signal peaks
                    dividedpeaks{end+1} = peak;
                    dividedtimes{end+1} = time;
                end    
            else % if not, the definition of a new peak is finished and it is added to the cell array with the signal peaks
                dividedpeaks{end+1} = peak; 
                dividedtimes{end+1} = time;
                newpeak = 1;
            end                             
        end
    end
end