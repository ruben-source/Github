function L = gradientfunc(x,lambda,alpha,gamma,s,n)
    x_T = x(1,:); x_M1 = x(2,:); x_M2 = x(3,:);
    if n == 1
        L = gamma(1)*(alpha(s,1,:) - alpha(1,1,:)) + lambda(1,:).*x_M1.*x_T;
    elseif n == 2
        L = gamma(2)*(alpha(s,2,:) - alpha(1,2,:)) - lambda(1,:).*x_M2.*x_T;
    elseif n == 3
        L = gamma(3)*(alpha(s,3,:) - alpha(1,3,:)) - (lambda(2,:).*x_M1.*x_T.*(1-(x_M1+x_T)/beta_M));
    elseif n == 4
        L = gamma(4)*(alpha(s,4,:) - alpha(1,4,:)) - (lambda(3,:).*x_M2.*x_T.*(1-(x_M1+x_T)/beta_M));
    else 
        L = gamma(5)*(alpha(s,5,:) - alpha(1,4,:)) + x_M1.*x_T.*(lambda(2,:) - lambda(3,:));
    end
end