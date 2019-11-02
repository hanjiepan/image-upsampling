% use the locaized mask function for image up-sampling
clear;
close all;
addpath('./image/');
addpath('tools/');

save_opt=0;
im_name='butterfly_256_256.tif';
'chip_200_244.tif';
'flowers_362_500.tif';
'house_256_256.tif';
'mandrill_512_512.tif';
'mit_256_256.tif';
'peppers_512_512.tif';
'ppt3_656_529.tif';
'temple_370_370.tif';
'zebra_391_586.tif';
[~,basename]=fileparts(im_name);
extension='circular';
downSamp_ratio=3;

if usejava('desktop')
    show_fig=1;
else
    show_fig=0;
end

im=imread(im_name);

color=(length(size(im))==3);
% crop ground truth image so that it is integer multiples of the
% up-sampling ratio
s_im=size(im)-mod(size(im),downSamp_ratio);
if color
    im=cat(3,im(1:s_im(1),1:s_im(2),1),...
        im(1:s_im(1),1:s_im(2),2),...
        im(1:s_im(1),1:s_im(2),3));
else
    im=im(1:s_im(1),1:s_im(2));
end

if color
    YCbCr=rgb2ycbcr(im);
    im=double(YCbCr(:,:,1));
else
    im=double(im);
end
% taper boundary to avoid artefacts due to periodic extension
im=double(edgetaper(im,fspecial('gaussian',3,1)));

% down-sampling
cubicFunc=@(x) (1.5*abs(x).^3 - 2.5*abs(x).^2 + 1) .* (abs(x) <= 1) + ...
                (-0.5*abs(x).^3 + 2.5*abs(x).^2 - 4*abs(x) + 2) .* ...
                ((1 < abs(x)) & (abs(x) <= 2));
h_hori=cubicFunc(linspace(-2,2,4*downSamp_ratio+1));
h_hori=h_hori(2:end-1);
h_bicubic=kron(h_hori,h_hori.');
h_bicubic=h_bicubic./sum(h_bicubic(:));
im_lp=imfilter(im,h_bicubic,extension,'conv');
% down-sampling
offset=floor(downSamp_ratio/2); % in order to be consistent with Matlab's imresize results
samp_grid=zeros(size(im_lp));
samp_grid(offset+1:downSamp_ratio:offset+end-(downSamp_ratio-1),...
    offset+1:downSamp_ratio:offset+end-(downSamp_ratio-1))=1;

lowRes_img=im_lp(offset+1:downSamp_ratio:offset+end-(downSamp_ratio-1),...
        offset+1:downSamp_ratio:offset+end-(downSamp_ratio-1));

%% up-sampling
[I_sup0,t_I0,obj0]=upSamp_ell2_admm(h_bicubic,samp_grid,lowRes_img,300);
I_sup0=max(min(I_sup0,255),0);

psnr0=PSNR(double(I_sup0(downSamp_ratio+1:end-downSamp_ratio,...
    downSamp_ratio+1:end-downSamp_ratio)),...
    im(downSamp_ratio+1:end-downSamp_ratio,...
    downSamp_ratio+1:end-downSamp_ratio),255);
ssim0=ssim_index(double(I_sup0(downSamp_ratio+1:end-downSamp_ratio,...
    downSamp_ratio+1:end-downSamp_ratio)),...
    im(downSamp_ratio+1:end-downSamp_ratio,...
    downSamp_ratio+1:end-downSamp_ratio));

if show_fig
    figure,imshow(I_sup0,[0,255])
end

%% obtain edge model by combining localized linear descriptor
[mask,coefs,pos_e,pos_smooth,sinTheta_all,cosTheta_all,...
    r_x_rescal_all,r_y_rescal_all,r_x,r_y,L_rescal_all,L_smooth,t_mask]=...
    edge_global_linearCons(I_sup0,[2,2],extension,show_fig,0.28);

%% up-sampling
ell=15;
h_data=h_bicubic;
[I_sup,t_up,obj]=upSamp_ell1_admm_direct(h_data,samp_grid,...
    lowRes_img,ell,mask,2000,I_sup0);

psnr_ell1=PSNR(double(I_sup(downSamp_ratio+1:end-downSamp_ratio,...
    downSamp_ratio+1:end-downSamp_ratio)),...
    im(downSamp_ratio+1:end-downSamp_ratio,...
    downSamp_ratio+1:end-downSamp_ratio),255);
ssim1=ssim_index(double(I_sup(downSamp_ratio+1:end-downSamp_ratio,...
    downSamp_ratio+1:end-downSamp_ratio)),...
    im(downSamp_ratio+1:end-downSamp_ratio,...
    downSamp_ratio+1:end-downSamp_ratio));

if show_fig
    figure,imshow(I_sup,[0,255])
end
save(['resUpSamp_',basename,'.mat'],'-v7.3')

%%
bicubic_res=imresize(lowRes_img,size(im),'bicubic');

psnr_bicubic=PSNR(double(bicubic_res(...
    downSamp_ratio+1:end-downSamp_ratio,...
    downSamp_ratio+1:end-downSamp_ratio)),...
    im(downSamp_ratio+1:end-downSamp_ratio,...
    downSamp_ratio+1:end-downSamp_ratio),255);
ssim_bicubic=ssim_index(double(bicubic_res(...
    downSamp_ratio+1:end-downSamp_ratio,...
    downSamp_ratio+1:end-downSamp_ratio)),...
    im(downSamp_ratio+1:end-downSamp_ratio,...
    downSamp_ratio+1:end-downSamp_ratio));

fprintf(['\n',basename,', ell=%.2e\n'],ell);
fprintf('PSNR_bicubic: %.2fdB,\tSSIM_bicubic: %.4f\n',psnr_bicubic,ssim_bicubic);
fprintf('PSNR0: %.2fdB,\tSSIM0: %.4f\n',psnr0,ssim0);
fprintf('PSNR1: %.2fdB,\tSSIM1: %.4f\n',psnr_ell1,ssim1);
