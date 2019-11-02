function p=PSNR(I1,I2,range,opt)
% PSNR: calculate peak signal-to-noise ratio
% I1 is the corrupted image with respect to I2
if nargin<3
    range=max(abs(I2(:)));
    mse=MSE(I1,I2);
elseif nargin<4
    mse=MSE(I1,I2);
else
    mse=MSE(I1,I2,opt);
end
p=10*log10(range^2/mse);
end