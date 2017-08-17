%my method of moments attempt for solving the deconvolution problem

function [output,coeff]=mom_1d_cart(x,t,b, delta_k, method)
%1d implementation
%input:
%x: x position of each grid point
%t: 1d psf function
%b: 1d image function
%delta_k: discretized k domain step
% method: basis function type. currently supports "point" and 'cos'

%output:
%output: the recovered 1d object function
% coeff: the set of coefficients calculated for the basis functions


max_dim=length(x);


% %basis matrix:
basis_mat=zeros(max_dim); 

for col=1:max_dim
   basis_func=cos(delta_k*(col-1)*x);
   
   %convolution F(g(x))      
    basis_mat(:,col)=conv(basis_func,t, 'same')';    
%     if mod(col,50)==0
%         figure, 
%         hold on 
%         plot(abs(fft(basis_func)))
%         hold off
%     end
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
        weight_func=cos(delta_k*(row-1)*x);
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
    output_temp=coeff(index)*cos(delta_k*(index-1)*x);

    output=output+output_temp;

end

end