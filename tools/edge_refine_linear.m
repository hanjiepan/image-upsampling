function [pos_refine,c_refine,lambda_refine]=edge_refine_linear(d_im_sq,sigmas,pos_e,c_e,lambda_e,extension)%#codegen
pos_refine=zeros(size(pos_e));
c_refine=c_e;%zeros(size(c_e));
lambda_refine=zeros(size(lambda_e));

r_x=sigmas(1);
r_y=sigmas(2);

s_im=size(d_im_sq);

rescal_factor_all=min((lambda_e(:,2)./lambda_e(:,1)).^0.25,1.5);
sigma_x_all=r_x.*rescal_factor_all;
sigma_y_all=r_y./rescal_factor_all;
sin_theta_all=c_e(:,1)./sqrt(c_e(:,1).*c_e(:,1)+c_e(:,2).*c_e(:,2));
cos_theta_all=c_e(:,2)./sqrt(c_e(:,1).*c_e(:,1)+c_e(:,2).*c_e(:,2));
w_blk_all=[max(round(2*r_y.*sin_theta_all),2),max(round(2*r_x.*cos_theta_all),2)];
x_ctr_all=round(pos_e(:,2));
y_ctr_all=round(pos_e(:,1));
pad_size=max(max(max([y_ctr_all,x_ctr_all]+w_blk_all,[],1)-s_im,...
    1+max(w_blk_all-[y_ctr_all,x_ctr_all],[],1)),[0,0]);
d_im_sq_pad=padarray(d_im_sq,pad_size,extension,'both');
x_round_error_all=pos_e(:,2)-x_ctr_all;
y_round_error_all=pos_e(:,1)-y_ctr_all;

for loop=1:numel(x_ctr_all)
    sigma_x=sigma_x_all(loop);
    sigma_y=sigma_y_all(loop);
    sin_theta=sin_theta_all(loop);
    cos_theta=cos_theta_all(loop);
    
    w_blk=w_blk_all(loop,:);
    [X_w,Y_w]=meshgrid(-w_blk(2):w_blk(2),-w_blk(1):w_blk(1));
    
    x_ctr=x_ctr_all(loop);
    y_ctr=y_ctr_all(loop);
    x_round_error=x_round_error_all(loop);
    y_round_error=y_round_error_all(loop);
	dImg_blk=d_im_sq_pad(pad_size(1)+y_ctr-w_blk(1):pad_size(1)+y_ctr+w_blk(1),...
                        pad_size(2)+x_ctr-w_blk(2):pad_size(2)+x_ctr+w_blk(2));
    
    W=window(X_w.*cos_theta-Y_w.*sin_theta+x_round_error,...
        X_w.*sin_theta+Y_w.*cos_theta+y_round_error,...
        sigma_x,sigma_y);
    dImg_flipW=dImg_blk(2*w_blk(1)+1:-1:1,2*w_blk(2)+1:-1:1).*W;
    dImg_flipWX=dImg_flipW.*X_w;
    dImg_flipWY=dImg_flipW.*Y_w;
    M11=sum2(dImg_flipWX.*X_w);
    M22=sum2(dImg_flipWY.*Y_w);
    M33=sum2(dImg_flipW);
    M12=sum2(dImg_flipWX.*Y_w);
    M13=sum2(-dImg_flipWX);
    M23=sum2(-dImg_flipWY);
    
    alpha=M33+eps;
    beta=M13*M13+M23*M23-M22*M33-M11*M33;
    gamma=M11*M22*M33-M11*M23*M23-M12*M12*M33-M13*M13*M22+2*M13*M12*M23;

    Delta=beta*beta-4*alpha*gamma;
    lambda1=(-beta-sqrt(Delta))/(2*alpha);
    lambda2=(-beta+sqrt(Delta))/(2*alpha);
    lambda=min(lambda1,lambda2)*(abs(Delta)>1e-10)+...
        (-beta/(2*alpha))*(abs(Delta)<=1e-10);
    
    c1=M13*(M22-lambda)-M12*M23;
    c2=(M11-lambda)*M23-M13*M12;
    c3=-(M11-lambda)*(M22-lambda)+M12*M12;
    norm_c1_c2=sqrt(c1*c1+c2*c2);
    c1=c1/norm_c1_c2;
    c2=c2/norm_c1_c2;
    c3=c3/norm_c1_c2;
    
    x_off=-c1*c3;
    y_off=-c2*c3;
    x_e=x_ctr+x_off;
    y_e=y_ctr+y_off;
    
    pos_refine(loop,:)=[y_e,x_e];
    c_refine(loop,:)=[c1,c2,c3];
    lambda_refine(loop,:)=[lambda1,lambda2];
end
end

function output=sum2(patch)
output=sum(patch(:));
end
