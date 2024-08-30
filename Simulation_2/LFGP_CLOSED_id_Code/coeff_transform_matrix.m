function [ Transform ] = coeff_transform_matrix( knots, p, d )
    Transform = zeros(p-1, p);
    for i = 1:(p-1)
        Transform(i,i) = -d/(knots(d+i+1)-knots(i+1));
        Transform(i,i+1) = d/(knots(d+i+1)-knots(i+1));
    end
end