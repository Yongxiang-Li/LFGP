function B = basismatrix(d, p, knots, U) 
%  
% Parameters:  
%   d	: Degree of the B-Spline. 
%   p   : Degree freedom of basis matrix.
%   knots	: Knot sequence, row vector of size nk. 
%   U	: Parametric evaluation points, row vector of size nu. 
% output: 
%   B	: Bspline basis matrix (nu, p) 
nu = numel(U); 
B = sparse(nu, p); 
for col=1:nu 
    s = findspan(p-1, d, U(col), knots);
    N = basisfun(s,U(col),d,knots);
    B(col, s-d+1:s+1) = N;
end
end

function s = findspan(n,p,u,U)                  
%  INPUT: 
%    n - number of control points - 1 
%    p - spline degree 
%    u - parametric point 
%    U - knot sequence 
%  
%  RETURN: 
%    s - knot span                                          
if (u==U(n+2)), s=n; return,  end                                                          
low = p;                                        
high = n + 1;                                   
mid = floor((low + high) / 2);                  
while (u < U(mid+1) || u >= U(mid+2))           
    if (u < U(mid+1))                           
        high = mid;                             
    else                                        
        low = mid;                                    
    end  
    mid = floor((low + high) / 2);              
end                                                                                      
s = mid;
end

function N = basisfun(i,u,p,U)                 
%    INPUT: 
%      i - knot span  ( from FindSpan() ) 
%      u - parametric point 
%      p - spline degree 
%      U - knot sequence 
%    OUTPUT: 
%      N - Basis functions vector[p+1]                                      
i = i + 1;                                                
left = zeros(p+1,1);                              
right = zeros(p+1,1);                                                                        
N(1) = 1;                                         
for j=1:p                                         
    left(j+1) = u - U(i+1-j);                     
    right(j+1) = U(i+j) - u;                      
    saved = 0;                                    
    for r=0:j-1                                   
        temp = N(r+1)/(right(r+2) + left(j-r+1)); 
        N(r+1) = saved + right(r+2)*temp;         
        saved = left(j-r+1)*temp;                 
    end                                           
    N(j+1) = saved;                               
end
end