function [I,t_consume,obj]=upSamp_ell2_admm(DataFilter,grid,I0,max_iter,I_ini)
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
    obj=zeros(max_iter,3);
end
grid=logical(grid);
s_I=size(grid);
s_I0=size(I0);
numel_I=numel(grid);
wd=floor(s_I./2);
[wx,wy]=meshgrid(-wd(2):s_I(2)-1-wd(2),-wd(1):s_I(1)-1-wd(1));
wx=ifftshift(wx).*sqrt(2*pi./numel_I);
wy=ifftshift(wy).*sqrt(2*pi./numel_I);

% preprocess the filter if it is of even size
[DataFilter,s_dataKernel]=modiFilter(DataFilter);
DataFilter_dual=conj(rot90(DataFilter,2));

Lap_ft=wx.^2+wy.^2;

Phit_I0=imfilter(upSamp(I0,grid,s_I),DataFilter_dual,'circular','conv');
lambda=1/max(DataFilter(:));

gamma=0.3*lambda;
LtL_ft_gamma=conj(Lap_ft).*Lap_ft+gamma;
if nargin>5
    I=I_ini;
else
    I=Phit_I0;
end
w=I;
v=ones(s_I);
z=ones(s_I0);

Phi_PhiT=conv2(circshift(padarray(DataFilter,s_I-s_dataKernel,0,'post'),...
    -floor(s_dataKernel./2)),...
    circshift(padarray(DataFilter_dual,s_I-s_dataKernel,0,'post'),...
    -floor(s_dataKernel./2)),'same');
PhiPhit_ft=fft2(ifftshift(reshape(downSamp(Phi_PhiT,grid),s_I0)));
gamma_PhiPhiT_ft=gamma+lambda.*PhiPhit_ft;

inv_NM=1/numel_I;

for count=1:max_iter
    % I update
    I_ft=(gamma.*fft2(w+v))./LtL_ft_gamma;
    I=ifft2(I_ft,'symmetric');
    % w update
    PhiT_I0z_Ir=lambda.*imfilter(upSamp(I0-z,grid,s_I),DataFilter_dual,...
        'circular','conv')+gamma.*(I-v);
    w=1/gamma.*(PhiT_I0z_Ir-lambda.*imfilter(upSamp(ifft2(fft2(reshape(downSamp(...
        imfilter(PhiT_I0z_Ir,DataFilter,'circular','conv'),...
        grid),s_I0))./gamma_PhiPhiT_ft),grid,s_I),DataFilter_dual,'circular','conv'));
    
    % update vector Lagrange multipliers
    Phi_w=reshape(downSamp(imfilter(w,DataFilter,'circular','conv'),grid),s_I0);
    z=z+(Phi_w-I0);
    v=v+w-I;
    if eval_obj
        obj(count,:)=[norm(I_ft.*Lap_ft,'fro')^2*inv_NM,...
            norm(I0-Phi_w,'fro'),norm(w-I,'fro')];
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
