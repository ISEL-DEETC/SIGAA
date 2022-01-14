function [blslopes] = getslopesfivemin(fnc,t)
    % this function estimates the slope every 5 minutes (or 300 secconds)
    fivemint=[1];    
    for(i=2:length(t))        
        if(t(i)-t(fivemint(end))>300)
            fivemint=[fivemint i-1];
        end
    end 
    blslopes = num2cell(diff(fnc(fivemint)')./diff(t(fivemint)'));
end