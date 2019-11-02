function [mask,coefs,pos_e,pos_smooth,sinTheta_all,cosTheta_all,...
    r_x_rescal_all,r_y_rescal_all,r_x,r_y,L_rescal_all,L_smooth,t_consume]=...
    edge_global_linearCons(im,sigmas,extension,showfig,th_canny)
t0=tic;
if nargin<3
    extension='circular';
end
if nargin<4
    showfig=1;
end
if nargin<5
    [~,th_canny]=edge(im,'canny');
    th_canny=max(1.8*th_canny(2),0.28);
end
r_x=sigmas(1);
r_y=sigmas(2);
s_im=size(im);
plt_size=s_im;
max_factor=1.7;

% compute the derivative image
wd=floor(s_im./2);
[wx,wy]=meshgrid(-wd(2):s_im(2)-1-wd(2),-wd(1):s_im(1)-1-wd(1));
wx=ifftshift(wx).*sqrt(2*pi./prod(s_im));
wy=ifftshift(wy).*sqrt(2*pi./prod(s_im));
Dr_ft=1j.*wy;
% when the size is even, then Dr and Di are pre-processed such that
% applying Dr/Di to an image I is equivalent to 
% ifft2(Dr.*fft2(I),'symmetric')
if rem(s_im(1),2)==0 && rem(s_im(2),2)==0
    Dr_ft(wd(1)+1,:)=[0,-wd(1)*1j.*ones(1,wd(2)-1),0,wd(1)*1j.*ones(1,wd(2)-1)].*sqrt(2*pi./prod(s_im));
elseif rem(s_im(1),2)==0 && rem(s_im(2),2)~=0
    Dr_ft(wd(1)+1,:)=[0,-wd(1)*1j.*ones(1,wd(2)),wd(1)*1j.*ones(1,wd(2))].*sqrt(2*pi./prod(s_im));
end

Di_ft=-1j.*wx;
if rem(s_im(2),2)==0
    Di_ft(:,wd(2)+1)=0;
end

d_im_sq=ifft2(Di_ft.*fft2(im)).^2+ifft2(Dr_ft.*fft2(im)).^2;
% each block is from -L to L
L=round(4*[r_y,r_x]);

[edgePt_subset,c_e_ref,lambda_e_ref,]=edge_linear(im,sigmas,extension,0,1,size(im),th_canny);

% remove block centers that are on the boundary of the image
boundary_remove=false(s_im);
boundary_remove(4:end-3,4:end-3)=true;

% disqualify points that have too small gradients as edge points
abs_d_im=sqrt(d_im_sq);
[y_cord_qualify,x_cord_qualify]=find(abs_d_im>th_canny*max(abs_d_im(:)) & boundary_remove);
keep=false(size(edgePt_subset,1),1);
for loop=1:numel(y_cord_qualify)
    idx_loop=(max(abs(edgePt_subset(:,1)-y_cord_qualify(loop)),...
        abs(edgePt_subset(:,2)-x_cord_qualify(loop)))<1);
    keep=keep | idx_loop;
end
not_keep=~keep;
edgePt_subset(not_keep,:)=[];
c_e_ref(not_keep,:)=[];
lambda_e_ref(not_keep,:)=[];
% select a subset of edge points that are at least 'sep'-apart
sep_edge=sqrt(2*log(2))*sqrt(r_x*r_y);%
stop=0;
bg_idx=1;
while ~stop
    blk_ctr=edgePt_subset(bg_idx,:);
    dist_loop=(edgePt_subset(:,1)-blk_ctr(1)).^2+...
        (edgePt_subset(:,2)-blk_ctr(2)).^2;
    idx_remove=dist_loop<sep_edge^2 & dist_loop~=0;
    edgePt_subset(idx_remove,:)=[];
    c_e_ref(idx_remove,:)=[];
    lambda_e_ref(idx_remove,:)=[];
    bg_idx=bg_idx+1;
    stop=(bg_idx>size(edgePt_subset,1));
end
L_smooth=round(4*[r_y,r_x]);
sep_smooth=sqrt(r_x*r_y);
[smooth_x_cord,smooth_y_cord]=meshgrid(-L_smooth(2):r_x:s_im(2)+L_smooth(2)+1,-L_smooth(1):r_y:s_im(1)+L_smooth(1)+1);
smoothPt_subset=[smooth_y_cord(:),smooth_x_cord(:)];
% further remove the smooth points that are in the neighborhood of the edge
% points
for loop=1:size(edgePt_subset,1)
    blk_ctr=edgePt_subset(loop,:);
    dist_loop=(smoothPt_subset(:,1)-blk_ctr(1)).^2+...
        (smoothPt_subset(:,2)-blk_ctr(2)).^2;
    idx_remove=dist_loop<sep_smooth^2;
    smoothPt_subset(idx_remove,:)=[];
end

y_cord=edgePt_subset(:,1);
x_cord=edgePt_subset(:,2);
N_pt=size(edgePt_subset,1);
rescal_factor_all=min((lambda_e_ref(:,2)./lambda_e_ref(:,1)).^0.25,max_factor);
r_x_rescal_all=r_x.*rescal_factor_all;
r_y_rescal_all=r_y./rescal_factor_all;
norm_c_e_all=sqrt(c_e_ref(:,1).*c_e_ref(:,1)+c_e_ref(:,2).*c_e_ref(:,2));
sinTheta_all=c_e_ref(:,1)./norm_c_e_all;
cosTheta_all=c_e_ref(:,2)./norm_c_e_all;
cri=(abs(sinTheta_all)<sin(pi/4));
L_rescal_all=zeros(N_pt,2);
L_rescal_all(cri,1)=[max(round(L(1)./rescal_factor_all(cri)),L(1))];
L_rescal_all(cri,2)=[max(round(L(2).*rescal_factor_all(cri)),L(2))];
L_rescal_all(~cri,1)=[max(round(L(1).*rescal_factor_all(~cri)),L(1))];
L_rescal_all(~cri,2)=[max(round(L(2)./rescal_factor_all(~cri)),L(2))];

smooth_y_cord=smoothPt_subset(:,1)+L_smooth(1)+1;
smooth_x_cord=smoothPt_subset(:,2)+L_smooth(2)+1;
% build smooth area mask function
[X_smooth,Y_smooth]=meshgrid(-L_smooth(2):L_smooth(2),-L_smooth(1):L_smooth(1));
% create Gaussian window
W_smooth=window(X_smooth,Y_smooth,r_x,r_y);
mask_smooth=zeros(plt_size+2.*(2.*L_smooth+1));
for count=1:numel(smooth_y_cord)
    mask_smooth(smooth_y_cord(count):smooth_y_cord(count)+2*L_smooth(1),...
        smooth_x_cord(count):smooth_x_cord(count)+2*L_smooth(2))=...
        mask_smooth(smooth_y_cord(count):smooth_y_cord(count)+2*L_smooth(1),...
        smooth_x_cord(count):smooth_x_cord(count)+2*L_smooth(2))+W_smooth;
end
mask_smooth=mask_smooth(2*L_smooth(1)+1+1:end-(2*L_smooth(1)+1),2*L_smooth(2)+1+1:end-(2*L_smooth(2)+1));

cons_val=-1;
r_x_rescal_sq_all=r_x_rescal_all.*r_x_rescal_all;
r_y_rescal_sq_all=r_y_rescal_all.*r_y_rescal_all;

[Mtx_ab,w1]=buildM(d_im_sq,x_cord,y_cord,sinTheta_all,cosTheta_all,r_x_rescal_sq_all,r_y_rescal_sq_all,L_rescal_all);
mask_dIm2=mask_smooth.*d_im_sq;
w=build_ObjOffset(mask_dIm2,x_cord,y_cord,sinTheta_all,cosTheta_all,r_x_rescal_all,r_y_rescal_all,L_rescal_all);
w_ab=w(1:2*N_pt);
ab=Mtx_ab\(w1.*(-cons_val)-w_ab);
coefs=[ab(1:N_pt),ab(N_pt+1:2*N_pt),cons_val.*ones(N_pt,1)];

pos_e=[y_cord,x_cord];
pos_smooth=[smoothPt_subset(:,1),smoothPt_subset(:,2)];

mask_edge=build_edgeMask(coefs,x_cord,y_cord,sinTheta_all,cosTheta_all,...
            r_x_rescal_sq_all,r_y_rescal_sq_all,L_rescal_all,plt_size,s_im);

mask=mask_smooth+mask_edge;
mask=min(abs(mask),max(mask_smooth(:)));
if showfig
    figure,imshow(abs(mask),[])
end
t_consume=toc(t0);
end
