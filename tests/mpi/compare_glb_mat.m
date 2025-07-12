mat0 = read_mat('./test5X/UZ_GLB_FFF.dat','FFF');
mat1 = read_mat('./test2X/UZ_GLB_FFF.dat','FFF');

% NOTE: in PPP space, if the numbers of the collocation points in some
% direction are different, the two 3d matrix files cannot be compared
% directly. We need to find the matching physical points and compare the
% values there.
% For example, in PPP-space, the azimuthal direction has NTH complex data
% points, which correspond to 2*(NTH+1) actual azimuthal angles stored in
% the real and imaginary parts of the complex data points:
%     TFM%THR = (/ (2*PI/NTH*(I-1), I=1,NTH+1) /)
%     TFM%THI = (/ (2*PI/NTH*(I-0.5D0), I=1,NTH+1) /)
% which is equivalent to: th = 2*pi/nth.*(0:nth) +
% 2*pi/nth.*((1:nth+1)-0.5).*1i. Therefore, if the 2nd run has exactly two
% times NTH of the first run, the real part of the complex data points of
% the 2nd run should match all the entries of the first run, i.e.
%     1st run: (1,2i), (3,4i), (5,6i), (7,8i), (9,  X)
%                 |       |       |       |       |
%     2nd run: (1,3 ), (5,7 ), (9,11),(13,15), (17, X)
% Sample code:
% mat1r = real(mat1);
% mat1c = mat1r(:,:,1:2:end-1)+mat1r(:,:,2:2:end-1).*1i;
% mat0c = mat0(:,:,1:end-1);

dmat = abs(mat0-mat1);
[maxdiff,loc] = max(dmat,[],'all')
[loc_1,loc_2,loc_3] = ind2sub(size(dmat),loc)

function mat3d = read_mat(filename,sp_flg)

formatSpec = "%f, ";
% filename = './test3X/DCHI_GLB_FFF.dat';

fileID = fopen(filename,'r');
mat3d = textscan(fileID,formatSpec);
fclose(fileID);

if nargin < 1
    sp_flg = 'fff';
end

if lower(sp_flg) == 'fff'
    % 356 x 64 x 255
    mat3d = reshape(mat3d{1,1},356,255,[]);
elseif lower(sp_flg) == 'ppp'
    % 373 x 129 x 256
    mat3d = reshape(mat3d{1,1},373,256,[]);
end

end