function p=SNR(I1,I2,opt)
% PSNR: calculate peak signal-to-noise ratio
% I1 is the corrupted image with respect to I2
if nargin<3
    p=20*log10(norm(I2,'fro')/norm(I1-I2,'fro'));
else
    p=10*log10(norm(I2,'fro')^2/(numel(I1)*MSE(I1,I2,opt)));
end
end