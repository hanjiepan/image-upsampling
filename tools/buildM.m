function [Mtx,w]=buildM(d_im_sq,x_cord,y_cord,sinTheta_all,cosTheta_all,r_x_rescal_sq_all,r_y_rescal_sq_all,L_rescal_all)%#codegen
s_im=size(d_im_sq);
N_pt=length(x_cord);
M_x_y_w2=zeros(N_pt,N_pt);
M_x_w2=zeros(N_pt,N_pt);
M_y_w2=zeros(N_pt,N_pt);
M_x2_w2=zeros(N_pt,N_pt);
M_y2_w2=zeros(N_pt,N_pt);
M_w2=zeros(N_pt,N_pt);

x1_ctr_all=min(max(round(x_cord),1),s_im(2));
y1_ctr_all=min(max(round(y_cord),1),s_im(1));
for loop_out=1:N_pt
    x1_ctr_subpixel=x_cord(loop_out);
    y1_ctr_subpixel=y_cord(loop_out);
    x1_ctr=x1_ctr_all(loop_out);
    y1_ctr=y1_ctr_all(loop_out);
    
    x1_round_error=x1_ctr-x1_ctr_subpixel;
    y1_round_error=y1_ctr-y1_ctr_subpixel;
    
    L_out=L_rescal_all(loop_out,:);
    r_x_rescal_sq_out=r_x_rescal_sq_all(loop_out);
    r_y_rescal_sq_out=r_y_rescal_sq_all(loop_out);
    sinTheta_out=sinTheta_all(loop_out);
    cosTheta_out=cosTheta_all(loop_out);

    x_offset_all=x1_ctr_subpixel-x_cord;
    y_offset_all=y1_ctr_subpixel-y_cord;
    inrange_all=find((abs(x_offset_all)<=L_rescal_all(:,2)+L_out(2)) & ...
        (abs(y_offset_all)<=L_rescal_all(:,1)+L_out(1)));
    total_idx_num=numel(inrange_all);
    
    for count_inner=1:total_idx_num
        loop_inner=inrange_all(count_inner);
        x_offset=x_offset_all(loop_inner);
        y_offset=y_offset_all(loop_inner);
        L_inner=L_rescal_all(loop_inner,:);

        r_x_rescal_sq_inner=r_x_rescal_sq_all(loop_inner);
        r_y_rescal_sq_inner=r_y_rescal_sq_all(loop_inner);
        sinTheta_inner=sinTheta_all(loop_inner);
        cosTheta_inner=cosTheta_all(loop_inner);

        Lout_2Lin=L_out+2.*L_inner;
        [X_s,Y_s]=meshgrid(-Lout_2Lin(2):Lout_2Lin(2),-Lout_2Lin(1):Lout_2Lin(1));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if y1_ctr-Lout_2Lin(1)<=0 || x1_ctr-Lout_2Lin(2)<=0 || ...
                y1_ctr+Lout_2Lin(1)>s_im(1) || x1_ctr+Lout_2Lin(2)>s_im(2)
            ver_sel=mod(y1_ctr-Lout_2Lin(1):y1_ctr+Lout_2Lin(1),s_im(1));
            ver_sel(ver_sel==0)=s_im(1);
            hor_sel=mod(x1_ctr-Lout_2Lin(2):x1_ctr+Lout_2Lin(2),s_im(2));
            hor_sel(hor_sel==0)=s_im(2);
            dImg_blk=d_im_sq(ver_sel,hor_sel);
        else
            dImg_blk=d_im_sq(y1_ctr-Lout_2Lin(1):y1_ctr+Lout_2Lin(1),...
                x1_ctr-Lout_2Lin(2):x1_ctr+Lout_2Lin(2));
        end
        X_s_x1e=X_s+x1_round_error;
        Y_s_y1e=Y_s+y1_round_error;
        X_s_x1e_x_offset=X_s_x1e+x_offset;
        Y_s_y1e_y_offset=Y_s_y1e+y_offset;
        
        W1_s=window(X_s_x1e.*cosTheta_out...
            -Y_s_y1e.*sinTheta_out,...
            X_s_x1e.*sinTheta_out...
            +Y_s_y1e.*cosTheta_out,...
            r_x_rescal_sq_out,r_y_rescal_sq_out,L_out);
        W2_s=window(X_s_x1e_x_offset.*cosTheta_inner...
            -Y_s_y1e_y_offset.*sinTheta_inner,...
            X_s_x1e_x_offset.*sinTheta_inner...
            +Y_s_y1e_y_offset.*cosTheta_inner,...
            r_x_rescal_sq_inner,r_y_rescal_sq_inner,L_inner);
        W1sW2s=W1_s.*W2_s;

        dImg_blkXW2=dImg_blk.*X_s_x1e.*W1sW2s;
        dImg_blkYW2=dImg_blk.*Y_s_y1e.*W1sW2s;
        M_x2_w2(loop_out,loop_inner)=...
            sum2(dImg_blkXW2.*...
            X_s_x1e_x_offset);
        M_y2_w2(loop_out,loop_inner)=...
            sum2(dImg_blkYW2.*...
            Y_s_y1e_y_offset);
        M_w2(loop_out,loop_inner)=...
            sum2(dImg_blk.*W1sW2s);
        M_x_y_w2(loop_out,loop_inner)=...
            sum2(dImg_blkXW2.*...
            Y_s_y1e_y_offset);
        M_x_w2(loop_out,loop_inner)=...
            sum2(dImg_blkXW2);
        M_y_w2(loop_out,loop_inner)=...
            sum2(dImg_blkYW2);
    end
end

if nargout==1
Mtx=[M_x2_w2,M_x_y_w2,M_x_w2;
    M_x_y_w2',M_y2_w2,M_y_w2;
    M_x_w2',M_y_w2',M_w2];
else
    Mtx=[M_x2_w2,M_x_y_w2;
        M_x_y_w2',M_y2_w2];
    w=sum([M_x_w2;M_y_w2],2);
end
end

function output=sum2(patch)
output=sum(patch(:));
end

function w_val=window(x,y,sigma_xsq,sigma_ysq,w_wd)
% create Gaussian window
idx=(abs(x)<=w_wd(2) & abs(y)<=w_wd(1));
w_val=zeros(size(x));
x_sub=x(idx);
y_sub=y(idx);
w_val(idx)=exp(-x_sub.*x_sub./(2*sigma_xsq)-y_sub.*y_sub./(2*sigma_ysq));
end
