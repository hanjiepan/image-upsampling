function [I,t_consume,obj]=upSamp_ell1_admm_direct(DataFilter,grid,I0,ell,mask,max_iter,I_ini)
% Implement the up-sampling problem with annihilation constraint (ell-1)
% with alternating direction method of multipliers (admm).
% 
% Input:    LapFilter: the Laplacian filter used in the smoothness
%               regularization
%           DataFilter: the sampling kernel that links the high-resolution
%               image with the given low-resolution image
%           grid: down-sampling grid
%           y: the given low-resolution image
%           ell: annihilation constraint regularization weight
%           mu: mask function, which contains the modelisation of the edges
%           max_iter: (optional) maximum number of iterations used in the
%               augmented Lagranian iteration.
% 
% Here we have implemented the ADMM version that involves three auxillary
% variables
% 
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

Dr_ft=1j.*wy;
% when the size is even, then Dr and Di are pre-processed such that
% applying Dr/Di to an image I is equivalent to 
% ifft2(Dr.*fft2(I),'symmetric')
if rem(s_I(1),2)==0 && rem(s_I(2),2)==0
    Dr_ft(wd(1)+1,:)=[0,-wd(1)*1j.*ones(1,wd(2)-1),0,wd(1)*1j.*ones(1,wd(2)-1)].*sqrt(2*pi./numel_I);
elseif rem(s_I(1),2)==0 && rem(s_I(2),2)~=0
    Dr_ft(wd(1)+1,:)=[0,-wd(1)*1j.*ones(1,wd(2)),wd(1)*1j.*ones(1,wd(2))].*sqrt(2*pi./numel_I);
end
Dr_conj_ft=conj(Dr_ft);

Di_ft=-1j.*wx;
if rem(s_I(2),2)==0
Di_ft(:,wd(2)+1)=0;
end
Di_conj_ft=conj(Di_ft);

Dr_jDi=Dr_ft+1j.*Di_ft;

DtD_ft=Dr_conj_ft.*Dr_ft+Di_conj_ft.*Di_ft;

Lap_ft=wx.^2+wy.^2;
Lap_conj_ft=conj(Lap_ft);
LtL_ft=Lap_ft.*Lap_conj_ft;

% preprocess the filter if it is of even size
DataFilter=modiFilter(DataFilter);
DataFilter_dual=conj(rot90(DataFilter,2));

% initialisation
if nargin>6
    I=I_ini;
else
    I=imfilter(upSamp(I0,grid,s_I),DataFilter_dual,'circular','conv');
end
I_ft=fft2(I);
inv_NM=1/numel_I;
z=ones(s_I0);
u_12=ifft2(I_ft.*Dr_jDi);
u1=real(u_12);
u2=imag(u_12);
v1=ones(s_I);
v2=v1;
w=I;
r=ones(s_I);

lambda=0.1*ell*max(mask(:))/max(DataFilter(:));20;%;
rho=.5*max(mask(:))*ell;%
gamma=rho;
% mask=mask(:);
th_level=ell/(2*rho).*abs(mask);

% preprocess the filter if it is of even size
[DataFilter,s_dataKernel]=modiFilter(DataFilter);
Phi_PhiT=conv2(circshift(padarray(DataFilter,s_I-s_dataKernel,0,'post'),...
    -floor(s_dataKernel./2)),...
    circshift(padarray(DataFilter_dual,s_I-s_dataKernel,0,'post'),...
    -floor(s_dataKernel./2)),'same');
PhiPhit_ft=fft2(ifftshift(reshape(downSamp(Phi_PhiT,grid),s_I0)));
gamma_PhiPhiT_ft=gamma+lambda.*PhiPhit_ft;

LtL_rhoDtD_gamma_ft=LtL_ft+rho.*DtD_ft+gamma;

for count=1:max_iter
    % I update
    I_ft=(rho.*(fft2(u1+v1).*Dr_conj_ft+fft2(u2+v2).*Di_conj_ft)...
        +gamma.*fft2(w+r))./LtL_rhoDtD_gamma_ft;
    I=ifft2(I_ft,'symmetric');
    
    % u update
    Dri_I=ifft2(I_ft.*Dr_jDi);
    Dr_I=real(Dri_I);
    Di_I=imag(Dri_I);
    dr=Dr_I-v1;
    di=Di_I-v2;
    norm_d=sqrt(dr.*dr+di.*di);
    u1=(dr./norm_d).*max(norm_d-th_level,0);
    u2=(di./norm_d).*max(norm_d-th_level,0);
    
    % w update
    PhiT_I0z_Ir=lambda.*imfilter(upSamp(I0-z,grid,s_I),DataFilter_dual,...
        'circular','conv')+gamma.*(I-r);
    w=1/gamma.*(PhiT_I0z_Ir-lambda.*imfilter(upSamp(ifft2(fft2(reshape(downSamp(...
        imfilter(PhiT_I0z_Ir,DataFilter,'circular','conv'),...
        grid),s_I0))./gamma_PhiPhiT_ft),grid,s_I),DataFilter_dual,'circular','conv'));

    % update vector Lagrange multipliers
    Phi_w=reshape(downSamp(imfilter(w,DataFilter,'circular','conv'),grid),s_I0);
    z=z+Phi_w-I0;
    v1=v1+u1-Dr_I;
    v2=v2+u2-Di_I;
    r=r+w-I;
    
    if eval_obj
        obj(count,:)=[norm(I_ft.*Lap_ft,'fro')^2*inv_NM+...
            ell*sum(sum(sqrt((mask.*Dr_I).^2+(mask.*Di_I).^2))),...
            norm(I0-Phi_w,'fro'),...
            norm(w-I,'fro')];
    end
end
t_consume=toc(t0);
end

function output=downSamp(input,grid)
output=input(grid);
end

function output=upSamp(input,grid,s_grid)
output=zeros(s_grid);
output(grid)=input;
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

