function mask_edge=build_edgeMask(coefs,x_cord,y_cord,sinTheta_all,cosTheta_all,r_x_rescal_sq_all,r_y_rescal_sq_all,L_rescal_all,plt_size,s_im)
% build edge mask
N_pt=length(x_cord);
mask_edge=zeros(plt_size);
for count=1:N_pt
    x1_ctr_subpixel=x_cord(count);
    y1_ctr_subpixel=y_cord(count);

    x1_ctr=round(x1_ctr_subpixel);
    y1_ctr=round(y1_ctr_subpixel);
    
    if x1_ctr<=0
        continue;
    end
    if x1_ctr>s_im(2)
        continue;
    end
    if y1_ctr<=0
        continue;
    end
    if y1_ctr>s_im(1)
        continue;
    end
    
    x1_round_error=x1_ctr-x1_ctr_subpixel;
    y1_round_error=y1_ctr-y1_ctr_subpixel;
    
    L_loop=L_rescal_all(count,:);
    r_x_rescal_sq=r_x_rescal_sq_all(count);
    r_y_rescal_sq=r_y_rescal_sq_all(count);
    sinTheta=sinTheta_all(count);
    cosTheta=cosTheta_all(count);
    [X_s,Y_s]=meshgrid(-L_loop(2):L_loop(2),-L_loop(1):L_loop(1));
    
    W1_s=window((X_s+x1_round_error).*cosTheta...
            -(Y_s+y1_round_error).*sinTheta,...
            (X_s+x1_round_error).*sinTheta...
            +(Y_s+y1_round_error).*cosTheta,...
            r_x_rescal_sq,r_y_rescal_sq,L_loop);%
    
    mask_edge=padarray(mask_edge,L_loop,0,'both');
    mask_edge(y1_ctr:y1_ctr+2*L_loop(1),x1_ctr:x1_ctr+2*L_loop(2))=...
        mask_edge(y1_ctr:y1_ctr+2*L_loop(1),x1_ctr:x1_ctr+2*L_loop(2))+...
        (coefs(count,1).*(X_s+x1_round_error)+...
        coefs(count,2).*(Y_s+y1_round_error)+coefs(count,3)).*W1_s;
    mask_edge=mask_edge(L_loop(1)+1:end-L_loop(1),L_loop(2)+1:end-L_loop(2));
end
end

function w_val=window(x,y,sigma_xsq,sigma_ysq,w_wd)
% create Gaussian window
idx=(abs(x)<=w_wd(2) & abs(y)<=w_wd(1));
w_val=zeros(size(x));
x_sub=x(idx);
y_sub=y(idx);
w_val(idx)=exp(-x_sub.*x_sub./(2*sigma_xsq)-y_sub.*y_sub./(2*sigma_ysq));
end
