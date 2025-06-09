optype = 'XM';
mkdir test2X/BAND
for ind = 1:128

    [band_n,band_logn] = readmat(ind,optype);

    if ~anynan(band_logn)
        writematrix(band_logn, ...
            append('./test2X/BAND/MATBAND_LOGLEG_',optype,sprintf('%03d.dat',ind)));
    end

    if ~anynan(band_n)
        writematrix(band_n, ...
            append('./test2X/BAND/MATBAND_LEG_',optype,sprintf('%03d.dat',ind)));
    end
end

function [band_n,band_logn] = readmat(ind,optype)

recycle('off')

formatSpec = "%f, ";
dim = 350+1-ind;

filename = append('./test2X/MATX_LEG_',optype,sprintf('%03d.dat',ind));
fileID = fopen(filename,'r');
if fileID == -1
    band_n = nan;
else
    band_n = textscan(fileID,formatSpec);
    band_n = double(reshape(band_n{1,1},dim,[]));
    band_n(:,dim+1:end) = 0;
    band_n = spdiags(band_n);
    fclose(fileID);
    if anynan(band_n)
        delete(filename);
    end
end

filename = append('./test2X/MATX_LOGLEG_',optype,sprintf('%03d.dat',ind));
fileID = fopen(filename,'r');
if fileID == -1
    band_logn = nan;
else
    band_logn = textscan(fileID,formatSpec);
    band_logn = double(reshape(band_logn{1,1},dim,[]));
    band_logn(:,dim+1:end) = 0;
    band_logn = spdiags(band_logn);
    fclose(fileID);
    if anynan(band_logn)
        delete(filename);
    end
end

recycle('on')

end