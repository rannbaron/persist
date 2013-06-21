function [ matching cost ] = kWass( a, b, k)

W=distMatrix(a,b).^k;
[m c]=m_hungarian(W);
matching=m;
cost=c^(1/k);

end
