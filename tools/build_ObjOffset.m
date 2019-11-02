function w=build_ObjOffset(mask_dIm2,x_cord,y_cord,sinTheta_all,cosTheta_all,r_x_rescal_all,r_y_rescal_all,L_rescal_all)
% build the linear offset term in the objective function
s_im=size(mask_dIm2);
N_pt=length(x_cord);
x_w_dIm2=zeros(N_pt,1);
y_w_dIm2=zeros(N_pt,1);
w_dIm2=zeros(N_pt,1);
x1_ctr_all=min(max(round(x_cord),1),s_im(2));
y1_ctr_all=min(max(round(y_cord),1),s_im(1));
for loop=1:N_pt
    x1_ctr_subpixel=x_cord(loop);
    y1_ctr_subpixel=y_cord(loop);
    x1_ctr=x1_ctr_all(loop);
    y1_ctr=y1_ctr_all(loop);
    x1_round_error=x1_ctr-x1_ctr_subpixel;
    y1_round_error=y1_ctr-y1_ctr_subpixel;
    
    L_loop=L_rescal_all(loop,:);
    r_x_rescal=r_x_rescal_all(loop);
    r_y_rescal=r_y_rescal_all(loop);
    sinTheta=sinTheta_all(loop);
    cosTheta=cosTheta_all(loop);
    [X_s,Y_s]=meshgrid(-L_loop(2):L_loop(2),-L_loop(1):L_loop(1));
    
    ver_sel=mod(y1_ctr-L_loop(1):y1_ctr+L_loop(1),s_im(1));
    ver_sel(ver_sel==0)=s_im(1);
    hor_sel=mod(x1_ctr-L_loop(2):x1_ctr+L_loop(2),s_im(2));
    hor_sel(hor_sel==0)=s_im(2);
    Img_blk=mask_dIm2(ver_sel,hor_sel);

    W1_s=window((X_s+x1_round_error).*cosTheta...
            -(Y_s+y1_round_error).*sinTheta,...
            (X_s+x1_round_error).*sinTheta...
            +(Y_s+y1_round_error).*cosTheta,...
            r_x_rescal,r_y_rescal,L_loop);

    x_w_dIm2(loop)=sum2(Img_blk.*(X_s+x1_round_error).*W1_s);
    y_w_dIm2(loop)=sum2(Img_blk.*(Y_s+y1_round_error).*W1_s);
    w_dIm2(loop)=sum2(Img_blk.*W1_s);
end
% w=[x_w_dIm2;y_w_dIm2];
w=[x_w_dIm2;y_w_dIm2;w_dIm2];
end

function output=sum2(patch)
output=sum(patch(:));
end