%my method of moments attempt for solving the deconvolution problem



function [output,coeff]=mom_2d_cart(x,y,t,B, delta_k, method)
%2d implementation: x,t,b, should be m*n arrays representing a grid
%square grids for now, without considering complex meshing

%input:
%x: x position of each grid point
%y: y position of each grid point
%t: 2d psf function
%B: 2d image function
%delta_k: discretized k domain step
% method: basis function type. currently supports "point" and 'cos'

%output:
%output: the recovered 2d object function
% coeff: the set of coefficients calculated for the basis functions


%flatten B into array
b=B(:);

max_dim=length(x(:));

%basis matrix:
basis_mat=zeros(max_dim);

for col=1:max_dim     
    basis_func=cos(delta_k*(col-1)*x).*cos(delta_k*(col-1)*y);

    %the following line may help with conditioning:
%    basis_func=basis_func/max(max(basis_func));
    
%convolution F(g(x))      
    
    basis_mat_temp=conv2(basis_func, t, 'same'); 

    basis_mat(:,col)=basis_mat_temp(:);
end

%weighting for different weight methods
basis_2=zeros(size(basis_mat));
b_2=zeros(size(b));
if strcmp(method, 'point')==1
    %point matching
    basis_2=basis_mat;
    b_2=b;
elseif strcmp(method, 'cos')==1
%cosine weights
    for row=1:max_dim
    %weighting the output 
        
        weight_func=cos(delta_k*(row-1)*x).*cos(delta_k*(row-1)*y);
        weight_func=weight_func(:)';
        
        b_2(row)=weight_func*b;
        
        %weighting the basis matrix
        for col=1:max_dim
            basis_2(row,col)=weight_func*basis_mat(:,col);
        end
    end


else
    error('this method is currently not supported')
end
cond(basis_2)
rank(basis_2)

coeff=basis_2\b_2;
coeff=coeff/max(max(coeff),-min(coeff)); %normalization

%construct the output using basis
output=zeros(size(x));

for index=1:max_dim
    output_temp=coeff(index)*cos(delta_k*(index-1)*x).*cos(delta_k*(index-1)*y);
    output=output+output_temp;
end

end