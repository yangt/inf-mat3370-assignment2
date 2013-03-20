function [C,xs]=simplex(A,b,c,n,m)
	itr = 0;
    Nb = 1:1:n;
	Bb = n+1:1:n+m;
	B = eye(m);
	N = A;
	xBs = b;
	zNs = -c;
	while (min(zNs) < 0)
		% step 2
		j = 1;
		zj = zNs(j,1);
		while zj >= 0
			zj = zNs(j+1,1);
			j = j+1;
		end
		% step 3
		ej = zeros(n,1);
		ej(j) = 1;
		dxB = B\(N*ej);
		% step 4
		[t i] = max(dxB./xBs);
		t = 1/t;
		% step 5
		if (length(i) > 1)
			i = i(1);
			t = t(1);
        end
		% step 6
		ei = zeros(m,1);
		ei(i) = 1;
        BiN = B\N;
        dzN  = -(BiN)'*ei;
		% step 7
		s = zNs(j)/dzN(j);
		% step 8
		xjs = t;
		xBs = xBs - t*dxB;
		xBs(i) = xjs;
		zis = s;
		zNs = zNs - s*dzN;
		zNs(j) = zis;
		% step 9
		colN = N(1:m,j);
		colB = B(1:m,i);
		N(1:m,j) = colB;
		B(1:m,i) = colN;	
		Ntmp = Nb(j);	
		Nb(j) = Bb(i);
		Bb(i) = Ntmp;
        % iteration and objective function
        C = 0;
        xs = zeros(n,1);
        for i=1:m
            if (Bb(i) <= n)
                 C = C + c(Bb(i))*xBs(i);
                 xs(Bb(i))=xBs(i);
            end
        end
        itr = itr + 1;
        fprintf('Simplex: Iteration=%d, C=%f\n',itr,C)
    end
    C = 0;
    xs = zeros(n,1);
    for i=1:m
        if (Bb(i) <= n)
            C = C + c(Bb(i))*xBs(i);
            xs(Bb(i))=xBs(i);
        end
    end
    fprintf('Simplex: Total Iterations=%d, C=%f\n',itr,C)
    
	
