function triData = computeTriangles( wavData )

s = length(wavData.m);
npairs = s*(s-1)/2;
gamma = zeros(npairs,1); 
displacement = zeros(npairs,1);
sym_wmp = zeros(npairs,1);
sym_wgt = zeros(npairs,1);
sym_hst = zeros(npairs,1);
sym_clr = zeros(npairs,1);
pntP = cell(npairs,1);
pntQ = cell(npairs,1);
rp = randperm(s);
itr = 1;
for jj = 1:s-1
    for kk = jj+1:s
        j = rp(jj);
        k = rp(kk);
        
        p = [wavData.x(j); wavData.y(j)];
        q = [wavData.x(k); wavData.y(k)];
        taup = [cos(wavData.a(j)); sin(wavData.a(j))];
        tauq = [cos(wavData.a(k)); sin(wavData.a(k))];
        T_pq = (q-p)/norm(q-p);
        T_perp_pq = [0 -1; 1 0]*T_pq;
        theta = ang(T_perp_pq(1),T_perp_pq(2));
        S = [cos(2*theta) sin(2*theta); sin(2*theta) -cos(2*theta)]; % reflection matrix
        wmp = abs(tauq'*S*taup); % weak mirror potential (the closer to 1, the more mirror symmetric)
        
        d = T_pq;
        if d(2) < 0
            d = -d;
        end
        gamma(itr) = atan2(d(2), d(1));
        displacement(itr) = dot((p+q)/2, d);
        sym_wmp(itr) = wmp;
        sym_wgt(itr) = wavData.m(j)*wavData.m(k);

        localV1 = wavData.v{j};
        indMaxOrV1 = wavData.vi(j);
        localV2 = wavData.v{k};
        indMaxOrV2 = wavData.vi(k);
        localV1 = localV1([indMaxOrV1:end 1:indMaxOrV1-1]);
        localV2_rev = localV2([indMaxOrV2:-1:1 end:-1:indMaxOrV2+1]);    
        histW_rev = 0;
        for z=1:numel(localV1)
            histW_rev = histW_rev + min([localV1(z),localV2_rev(z)]);
        end
        sym_hst(itr) = histW_rev;
        pntP{itr} = p;
        pntQ{itr} = q;

        hc1 = wavData.hC{j};
        hc2 = wavData.hC{k};
        histIntrC = 0;
        for z=1:numel(hc1)
            histIntrC = histIntrC + min([hc1(z),hc2(z)]);
        end
        sym_clr(itr) = histIntrC;

        itr = itr+1;
    end
end

triData.gamma = gamma;
triData.displacement = displacement;
triData.sym_wmp = sym_wmp;
triData.sym_wgt = sym_wgt;
triData.sym_hst = sym_hst;
triData.sym_clr = sym_clr;
triData.p = pntP;
triData.q = pntQ;

end

