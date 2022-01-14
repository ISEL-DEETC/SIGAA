function [time] = gettime(datfile)
fid = fopen(datfile);
tline = fgetl(fid);
time=[];
while ischar(tline)
%    if (contains(tline,'Clock reset to 0.0'))
    if (~isempty(strfind(tline,'Clock reset to 0.0')))

        tline = fgetl(fid);
        while ischar(tline)
            tmp = regexp(tline,' ','split');
            time = [time str2double(tmp{1,length(tmp)})/double(100)];
            tline = fgetl(fid);
        end
        break;
    end
    tline = fgetl(fid);
end
fclose(fid);
end