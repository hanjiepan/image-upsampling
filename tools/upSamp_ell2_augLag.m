function [x,t_consume,obj]=upSamp_ell2_augLag(LapFilter,DataFilter,grid,y,max_iter,x_ini)
% Implement the up-sampling problem subject to data-fidelity constraint
% only efficiently with augmented Lagrangian.
% 
% Input:    LapFilter: the Laplacian filter used in the smoothness
%               regularization
%           DataFilter: the sampling kernel that links the high-resolution
%               image with the given low-resolution image
%           grid: down-sampling grid
%           y: the given low-resolution image
%           max_iter: (optional) maximum number of iterations used in the
%               augmented Lagranian iteration.
t0=tic;
warning('off','MATLAB:nearlySingularMatrix');
if nargout<3
    eval_obj=false;
else
    eval_obj=true;
    obj=zeros(max_iter,2);
end
grid=logical(grid);
s_x=size(grid);
s_y=size(y);

% preprocess the filter if it is of even size
[LapFilter,s_lap]=modiFilter(LapFilter);
[DataFilter,s_dataKernel]=modiFilter(DataFilter);

Lap_ft=fft2(circshift(padarray(LapFilter,s_x-s_lap,0,'post'),...
    -floor(s_lap./2)));
Lap_t_Lap_ft=conj(Lap_ft).*Lap_ft;
Lap_t_Lap=conv2(LapFilter,LapFilter);

A_ss=@(im) reshape(downSamp(imfilter(im,DataFilter,...
    'circular','conv'),grid),s_y);
A_sf=@(im) fft2(reshape(downSamp(imfilter(im,DataFilter,...
    'circular','conv'),grid),s_y));

DataFilter_dual=conj(rot90(DataFilter,2));
At_ss=@(im) imfilter(upSamp(im,grid,s_x),DataFilter_dual,...
    'circular','conv');
At_fs=@(im_ft) imfilter(upSamp(ifft2(im_ft),grid,s_x),DataFilter_dual,...
    'circular','conv');

Phi_PhiT=conv2(circshift(padarray(DataFilter,s_x-s_dataKernel,0,'post'),...
    -floor(s_dataKernel./2)),...
    circshift(padarray(DataFilter_dual,s_x-s_dataKernel,0,'post'),...
    -floor(s_dataKernel./2)),'same');
AAt_ft=fft2(ifftshift(reshape(downSamp(Phi_PhiT,grid),s_y)));
I_lamb_AtA_inv=@(im,reg) ...
    im-reg.*(At_fs(A_sf(im)./(1+reg.*AAt_ft)));

At_y=At_ss(y);

z=ones(s_y);

lambda=20;
if nargin>5
    x=x_ini;
else
    x=At_y;
end

inner_max=5;
num_basis=8;
F=zeros(numel(x),num_basis);
DF=zeros(numel(x),num_basis);
AF=zeros(numel(y),num_basis);

Ax=A_ss(x);
for count=1:max_iter
    At_z=At_ss(z);
    for inner=1:inner_max
        % build LET basis
        F1=x;
        F2_1=imfilter(x,Lap_t_Lap,'circular','conv')+At_z;
%         F2_2=At_sIN_sOUT(A_sIN_sOUT(x))-At_y;
        F2_2=At_ss(Ax)-At_y;
        F2=F2_1+lambda.*F2_2;
        
        F3=ifft2(fft2(F2_1+0.1*lambda.*F2_2)./(Lap_t_Lap_ft+0.1*lambda));
        F4=ifft2(fft2(F2)./(Lap_t_Lap_ft+lambda));
        F5=ifft2(fft2(F2_1+10*lambda.*F2_2)./(Lap_t_Lap_ft+10*lambda));
        
        F6=I_lamb_AtA_inv(F2_1+0.1*lambda.*F2_2,0.1*lambda);
        F7=I_lamb_AtA_inv(F2,lambda);
        F8=I_lamb_AtA_inv(F2_1+10*lambda.*F2_2,10*lambda);
        
        F(:,1)=F1(:);F(:,2)=F2(:);F(:,3)=F3(:);F(:,4)=F4(:);
        F(:,5)=F5(:);F(:,6)=F6(:);F(:,7)=F7(:);F(:,8)=F8(:);
        
        F_12=F1+1j.*F2;F_34=F3+1j.*F4;F_56=F5+1j.*F6;F_78=F7+1j.*F8;
        
        DF12=reshape(imfilter(F_12,LapFilter,'circular','conv'),[],1);
        DF34=reshape(imfilter(F_34,LapFilter,'circular','conv'),[],1);
        DF56=reshape(imfilter(F_56,LapFilter,'circular','conv'),[],1);
        DF78=reshape(imfilter(F_78,LapFilter,'circular','conv'),[],1);
        DF(:,1)=real(DF12);DF(:,2)=imag(DF12);
        DF(:,3)=real(DF34);DF(:,4)=imag(DF34);
        DF(:,5)=real(DF56);DF(:,6)=imag(DF56);
        DF(:,7)=real(DF78);DF(:,8)=imag(DF78);

        AF12=reshape(A_ss(F_12),[],1);
        AF34=reshape(A_ss(F_34),[],1);
        AF56=reshape(A_ss(F_56),[],1);
        AF78=reshape(A_ss(F_78),[],1);
        AF(:,1)=real(AF12);AF(:,2)=imag(AF12);
        AF(:,3)=real(AF34);AF(:,4)=imag(AF34);
        AF(:,5)=real(AF56);AF(:,6)=imag(AF56);
        AF(:,7)=real(AF78);AF(:,8)=imag(AF78);
        
        M_coef=(DF'*DF)+lambda.*(AF'*AF);
        coef=M_coef\(AF'*(lambda.*y(:)-z(:)));
        x=reshape(F*coef,s_x);
        Ax=reshape(AF*coef,s_y);
    end
    z=z+lambda.*(Ax-y);
    if eval_obj
        obj(count,:)=[norm(imfilter(x,LapFilter,'circular','conv'),'fro')^2,...
            norm(y-Ax,'fro')];
    end
end
t_consume=toc(t0);
end

function output=downSamp(input,grid)
output=input(grid);
end

function output=upSamp(input,grid,s_grid)
output=zeros(s_grid);
output(logical(grid))=input;
end

function [h,h_size]=modiFilter(h)
s_h=size(h);
if mod(s_h(1),2)==0
    h=padarray(h,[1,0],0,'post');
end
if mod(s_h(2),2)==0
    h=padarray(h,[0,1],0,'post');
end
h_size=size(h);
end
