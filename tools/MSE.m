function e=MSE(I1,I2,opt)
% MSE: calculate mean-square-error
s1=size(I1);
s2=size(I2);
if s1~=s2
    error('Size of two images does not match');
elseif nargin>2
    I1=[I1(:),ones(numel(I1),1)];
    coef=(I1'*I1)\(I1'*I2(:));
    I1=I1*coef;
    e=(I1-I2(:))'*(I1-I2(:))/prod(s1);
else
    e=(I1(:)-I2(:))'*(I1(:)-I2(:))/prod(s1);
end