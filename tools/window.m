function w_val=window(x,y,sigma_x,sigma_y,w_wd)%#codegen
% create Gaussian window
if nargin==5
    idx=(abs(x)<=w_wd(2) & abs(y)<=w_wd(1));
    w_val=zeros(size(x));
    x_sub=x(idx);
    y_sub=y(idx);
    w_val(idx)=exp(-x_sub.*x_sub./(2*sigma_x*sigma_x)-y_sub.*y_sub./(2*sigma_y*sigma_y));
else
    w_val=exp(-x.*x./(2*sigma_x*sigma_x)-y.*y./(2*sigma_y*sigma_y));
end
return