% use the locaized mask function for bind image up-sampling
clear;
close all;
addpath('./image/');
addpath('tools/');
save_opt=1;
im_name='sculpture.png';
[~,basename]=fileparts(im_name);
extension='circular';
downSamp_ratio=3;
% print resolution
dpi='300';
frame=[293,123,60,40];
zoom_factor=0.7;

if usejava('desktop')
    show_fig=1;
else
    show_fig=0;
end

im=imread(im_name);
if length(size(im))>2
    YCbCr=rgb2ycbcr(im);
    lowRes_img=double(YCbCr(:,:,1));
else
    lowRes_img=double(im);
end
% taper boundary to avoid artefacts due to periodic extension
lowRes_img=double(edgetaper(lowRes_img,fspecial('gaussian',2*downSamp_ratio+1,1)));
% assume bicubic kernel
cubicFunc=@(x) (1.5*abs(x).^3 - 2.5*abs(x).^2 + 1) .* (abs(x) <= 1) + ...
                (-0.5*abs(x).^3 + 2.5*abs(x).^2 - 4*abs(x) + 2) .* ...
                ((1 < abs(x)) & (abs(x) <= 2));
h_hori=cubicFunc(linspace(-2,2,4*downSamp_ratio+1));
h_hori=h_hori(2:end-1);
h_bicubic=kron(h_hori,h_hori.');
h_bicubic=h_bicubic./sum(h_bicubic(:));
% sampling grid
samp_grid=zeros(size(lowRes_img).*downSamp_ratio);
offset=floor(downSamp_ratio/2);
samp_grid(offset+1:downSamp_ratio:offset+end-(downSamp_ratio-1),...
    offset+1:downSamp_ratio:offset+end-(downSamp_ratio-1))=1;

%% up-sampling
if ~exist([basename,'_blind.mat'],'file')
I_sup0=upSamp_ell2_admm(h_bicubic,samp_grid,lowRes_img,200);
I_sup0=max(min(I_sup0,255),0);
else
    load([basename,'_blind.mat'],'I_sup0');
end

%%
bicubic_res=imresize(lowRes_img,size(samp_grid),'bicubic');

%% obtain edge model by combining localized linear descriptor
if ~exist([basename,'_blind.mat'],'file')
[mask,coefs,pos_e,pos_smooth,sinTheta_all,cosTheta_all,...
    r_x_rescal_all,r_y_rescal_all,r_x,r_y,L_rescal_all,L_smooth,t_mask]=...
    edge_global_linearCons(I_sup0,[2,2],extension,show_fig,0.28);%
else
    load([basename,'_blind.mat'],'mask');
end
%% up-sampling
ell=15;
h_data=h_bicubic;
snr_allowed=35;
epsilon=norm(lowRes_img(:),2)/(10^(snr_allowed/20));
if ~exist([basename,'_blind.mat'],'file')
I_sup=upSamp_ell1_noisy(h_data,samp_grid,...
    lowRes_img,ell,mask,200,epsilon);
save([basename,'_blind.mat']);
else
    load([basename,'_blind.mat'],'I_sup');
end

%% zoom-in plot
if show_fig
    figure(1)
    % enlarge for display purposes
    imshow(imresize(lowRes_img,size(samp_grid),'nearest'),[0,255])
    set(gcf,'Name','Low Resolution Image')
    if save_opt
        set(gcf,'paperpositionmode','auto');
        print(gcf,['/Volumes/RamDisk/',basename,'lowRes.png']...
            ,'-dpng',['-r',dpi])
    end
    
    figure(2),imshow(bicubic_res,[0,255])
    set(gcf,'Name','Bicubic Interpolation')
    if save_opt
        set(gcf,'paperpositionmode','auto');
        print(gcf,['/Volumes/RamDisk/',basename,'bicubic.png']...
            ,'-dpng',['-r',dpi])
    end
    
    figure(3)
    imshow(I_sup,[0,255])
    set(gcf,'Name','Annihilation-driven Approach')
    if save_opt
        set(gcf,'paperpositionmode','auto');
        print(gcf,['/Volumes/RamDisk/',basename,'AnniRes.png']...
            ,'-dpng',['-r',dpi])
    end

end
