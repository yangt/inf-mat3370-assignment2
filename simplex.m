function [C, iterations] = simplex_new(A, b, c)
    n = size(A,2);
    m = size(b,1);
    iterations = 0;
    N_indices = 1:1:n;
    B_indices = n+1:1:m+n;
    size(B_indices)
    B = eye(m);
    N = A;
    xs_B = b;
    zs_N = -c;
    while(min(zs_N) < 0)
        % step 2
        [tmp j] = min(zs_N);
        % step 3
        e_j = zeros(n, 1);
        e_j(j) = 1;
        dx_B = B\(N*e_j);
        % step 4 and 5
        [maximum i] = max(dx_B./xs_B);
        if(maximum <= 0)
            printf('unbounded')
            break
        end
        t = 1./maximum;
        % step 6
        e_i = zeros(m, 1);
        e_i(i) = 1;
        dz_N = -(B\N)'*e_i;
        % step 7
        s = zs_N(j)/dz_N(j);
        % step 8
        xs_B = xs_B - t*dx_B;
        xs_B(i) = t;
        zs_N = zs_N - s*dz_N;
        zs_N(j) = s; 
        % step 9
        tmp = B(:,i);
        B(:,i) = N(:, j);
        N(:,j) = tmp;
        tmp = B_indices(i);
        B_indices(i) = N_indices(j);
        N_indices(j) = tmp;
        iterations = iterations + 1;
    end
    C = 0;
    for i=1:m
        if (B_indices(i) <= n)
            C = C + c(B_indices(i))*xs_B(i);
        end
    end
end 
