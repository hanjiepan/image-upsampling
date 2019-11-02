function [pos_e,c_e,lambda_e,c_all,pos_iso,pos_canny]=edge_linear(im,sigmas,extension,showfig,refinement,plt_size,th_canny)
if nargin<3
    extension='replicate';
end
if nargin<4
    showfig=1;
end
if nargin<5
    refinement=0;
end
if nargin<6
    plt_size=size(im);
end
if nargin<7
    [~,th_canny]=edge(im,'canny');
    th_canny=max(1.8*th_canny(2),0.25);
end
r_x=sigmas(1);
r_y=sigmas(2);
s_im=size(im);

wd=floor(s_im./2);
[X,Y]=meshgrid(-wd(2):s_im(2)-1-wd(2),-wd(1):s_im(1)-1-wd(1));

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

Di_ft=1j.*wx;%-
if rem(s_im(2),2)==0
Di_ft(:,wd(2)+1)=0;
end

canny_edge=edge(im,'canny',th_canny);%

% compute the derivative image
dx_im=ifft2(Di_ft.*fft2(im));
dy_im=ifft2(Dr_ft.*fft2(im));
d_im_sq=dx_im.^2+dy_im.^2;

L=[5*r_y,5*r_x];
[X_w,Y_w]=meshgrid(-L(2):L(2),-L(1):L(1));
% create Gaussian window
W=window(X_w,Y_w,r_x,r_y);

M11=imfilter(d_im_sq,X_w.^2.*W,extension,'conv');
M22=imfilter(d_im_sq,Y_w.^2.*W,extension,'conv');
M33=imfilter(d_im_sq,W,extension,'conv');
M12=imfilter(d_im_sq,X_w.*Y_w.*W,extension,'conv');
M13=imfilter(d_im_sq,-X_w.*W,extension,'conv');
M23=imfilter(d_im_sq,-Y_w.*W,extension,'conv');

alpha_all=M33+eps;
beta_all=M13.^2+M23.^2-M22.*M33-M11.*M33;
gamma_all=M11.*M22.*M33-M11.*M23.^2-M12.^2.*M33-M13.^2.*M22+2.*M13.*M12.*M23;

Delta=beta_all.^2-4.*alpha_all.*gamma_all;
lambda1=(-beta_all-sqrt(Delta))./(2*alpha_all);
lambda2=(-beta_all+sqrt(Delta))./(2*alpha_all);
lambda=min(lambda1,lambda2).*(abs(Delta)>1e-10)+...
    (-beta_all./(2*alpha_all)).*(abs(Delta)<=1e-10);
c1=M13.*(M22-lambda)-M12.*M23;
c2=(M11-lambda).*M23-M13.*M12;
c3=-(M11-lambda).*(M22-lambda)+M12.^2;
norm_c1_c2=sqrt(c1.^2+c2.^2);
c1=c1./norm_c1_c2;
c2=c2./norm_c1_c2;
c3=c3./norm_c1_c2;

% extract points around Canny edge point
[y_cord,x_cord]=find(canny_edge==1);

c1_e=c1(y_cord+(x_cord-1)*s_im(1));
c2_e=c2(y_cord+(x_cord-1)*s_im(1));
c3_e=c3(y_cord+(x_cord-1)*s_im(1));

x_off=-c1_e.*c3_e;
y_off=-c2_e.*c3_e;

x_e=X(y_cord+(x_cord-1)*s_im(1))+x_off;
y_e=Y(y_cord+(x_cord-1)*s_im(1))+y_off;

lambda1_e=lambda1(y_cord+(x_cord-1)*s_im(1));
lambda2_e=lambda2(y_cord+(x_cord-1)*s_im(1));

pos_e=[y_e+wd(1)+1,x_e+wd(2)+1];
pos_canny=[Y(y_cord+(x_cord-1)*s_im(1))+wd(1)+1,...
            X(y_cord+(x_cord-1)*s_im(1))+wd(2)+1];

c_e=[c1_e,c2_e,c3_e];
lambda_e=[lambda1_e,lambda2_e];
c_all=[c1(:),c2(:),c3(:)];

% refine detected edge points' parameters with adaptived window size
if refinement
    pos_iso=pos_e;
    for further=1:5
        [pos_e,c_e,lambda_e]=edge_refine_linear(d_im_sq,sigmas,pos_e,c_e,lambda_e,extension);
    end
    % sign corrections
    x_pos_round=round(pos_e(:,2));
    y_pos_round=round(pos_e(:,1));
    idx_remove=x_pos_round>s_im(2) | x_pos_round<1 | y_pos_round>s_im(1) | y_pos_round<1;
    x_pos_round(idx_remove)=[];
    y_pos_round(idx_remove)=[];
    dx_normalised=dx_im./sqrt(dx_im.^2+dy_im.^2);
    dy_normalised=dy_im./sqrt(dx_im.^2+dy_im.^2);
    img_idx=y_pos_round+(x_pos_round-1)*s_im(1);
    c_e(idx_remove,:)=[];
    pos_e(idx_remove,:)=[];
    lambda_e(idx_remove,:)=[];
    inner_product_all=sign(c_e(:,1).*dx_normalised(img_idx)+...
        c_e(:,2).*dy_normalised(img_idx));
    c_e=bsxfun(@times,c_e,inner_product_all);

    if showfig
        iptsetpref('ImshowBorder','tight');
        figure,imshow(im,[0,255])
        hold on
        plot(pos_iso(:,2),pos_iso(:,1),'g+','MarkerSize',10,'LineWidth',2)
        plot(pos_e(:,2),pos_e(:,1),'y.','MarkerSize',20)
        plot(pos_canny(:,2),pos_canny(:,1),'w^','MarkerSize',10,'LineWidth',2)
        hold off
        set(gcf,'Name','Linear Model (Refined)')
    end
else
    if showfig
        iptsetpref('ImshowBorder','tight');
        figure,imshow(im,[0,255])
        hold on
        plot(pos_e(:,2),pos_e(:,1),'y.','MarkerSize',20)
        plot(pos_canny(:,2),pos_canny(:,1),'w^','MarkerSize',10,'LineWidth',2)
        hold off
        set(gcf,'Name','Linear Model')
    end
end
return

