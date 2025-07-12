% TEST NORM AND LOGNORM
% TFM%NORM(NRCHOPDIM+14,NTCHOPDIM): normalization factors of P^m_n
% TFM%LOGNORM(NRCHOPDIM+14,NTCHOPDIM): exp(lognorm) = norm

% READ NORM AND LOGNORM
fileID = fopen('NORMS.dat','r');
formatSpec = "%f, ";
norms = textscan(fileID,formatSpec);
norms = double(reshape(norms{1,1},370,128)');
fclose(fileID);

fileID = fopen('LOGNORMS.dat','r');
formatSpec = "%f, ";
lognorms = textscan(fileID,formatSpec);
lognorms = double(reshape(lognorms{1,1},370,128)');
fclose(fileID);

% NORM(J,I+1)/NORM(J,I)
ratio_norms = norms(:,2:end)./norms(:,1:end-1);
ratio_lognorms = exp(lognorms(:,2:end)-lognorms(:,1:end-1));

format longEng
for i = 1:size(ratio_lognorms,1)
    disp(sum(ratio_norms(i,:)-ratio_lognorms(i,:)));
end

