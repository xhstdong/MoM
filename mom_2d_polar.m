%my method of moments attempt for solving the deconvolution problem

function [output,coeff]=mom_2d_polar(r,t,B, delta_k, method)
%2d implementation: r,b, should be m*n arrays representing a grid
%square grids for now, without considering complex meshing

%we have assumed radial symmetry for now for simplicity

%input:
%r: r position of each grid point
%t: 2d psf function
%B: 2d image function
%delta_k: discretized k domain step
% method: basis function type. currently supports "point" and 'cos'

%output:
%output: the recovered 2d object function
% coeff: the set of coefficients calculated for the basis functions

%flatten B into array
b=B(:);

max_dim=length(r(:));
%basis matrix:
basis_mat=zeros(max_dim);

for col=1:max_dim
    %various possible basis functions:
    %basis_func=cos(delta_k*(col-1)*r); %could change here
    basis_func=besselj(0,delta_k*(col-1/2)*r);
    

    %scaling may help conditioning
%    basis_func=basis_func/max(max(basis_func));
    
%convolution F(g(x))      
    
%strictly speaking this is incorrect; however, we seem to get decent
%results, will be improved in the future:
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
elseif strcmp(method, 'bes')==1
%bessel weights
    for row=1:max_dim

    %weighting the output  

        %weight_func=cos(delta_k*(row-1)*r);
        weight_func=besselj(0,delta_k*(row-1/2)*r);
%         weight_func=weight_func/weight_func(round(size(r,1)/2), round(size(r,2)/2)); %scaling
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

output=zeros(size(r));
for index=1:max_dim
    %output_temp=coeff(index)*cos(delta_k*(index-1)*r);
    output_temp=coeff(index)*besselj(0,delta_k*(index-1/2)*r);
    %output_temp=output_temp/output_temp(round(size(r,1)/2),round(size(r,2)/2));  %scaling
    
    output=output+output_temp;
end

end