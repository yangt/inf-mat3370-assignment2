function [C,x]=interior(A,b,c,n,m)
    itr = 0;
    maxitr = 1000;
	x = ones(n,1);
    w = ones(m,1);
    y = w;
    z = x;
    X1 = 1;
    eps = 0.000000001;
    while(abs(b'*y - c'*x) > eps & itr <= maxitr)
        rho = b - A*x - w;
        sigma = c - A'*y + z;
        gamma = z'*x + y'*w;
        delta = 1/10;
        mu = delta*(gamma/(n+m));
        X = zeros(n);
        W = zeros(m);
        Y = W;
        Z = X;
        for i=1:n
            X(i,i) = x(i);
            Z(i,i) = z(i);
        end
        for j=1:m
            Y(j,j) = y(j);
            W(j,j) = w(j);
        end
        e1 = ones(m,1);
        e2 = ones(n,1);
        iYW = Y\W;
        iZAt = Z\A';
        iYe = Y\e1;
        iXe = X\e2;
        y1 = Z\(c - A'*y + mu*iXe);
        dy = (-(iYW + A*X*iZAt))\(b - A*x - mu*iYe - A*X*y1);
        dx = X*(Z\(c - A'*y + mu*iXe - A'*dy));
        dz = X\(mu*e2 - X*Z*e2 - Z*dx);
        dw = Y\(mu*e1 - Y*W*e1 - W*dy);
        r = 0.99;
        theta = min((r/max([max(-dx./x),...
                    max(-dw./w),max(-dy./y),...
                    max(-dz./z)])),1);
        x = x + theta*dx;
        y = y + theta*dy;
        w = w + theta*dw;
        z = z + theta*dz;
        itr = itr + 1;
        C = c'*x;
        fprintf('Interior-point:Iteration=%d, C =%f\n',itr,C)
    end
    C = c'*x;
    fprintf('Interior-point:Total Iterations=%d, C =%f\n',itr,C)
    
