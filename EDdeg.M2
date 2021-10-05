

needsPackage "Quasidegrees"

A=  matrix(
   {{0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1},
   {0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1},
   {1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1},
   {1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1},
   {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1},
   {3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0}});

   A = matrix {{1, 1, 1, 0, 0, 0},
        {1, 0, 0, 1, 1, 0},
        {0, 1, 0, 1, 0, 1},
        {0, 0, 1, 0, 1, 1}};

u =

R = QQ[a..f];
I = toricIdeal(A,R);
sing = I + minors(codim I,jacobian(I));
M = (matrix{apply(# gens R,i->(gens R)_i-u_i)})||(transpose(jacobian I));
J = saturate(I + minors((codim I)+1,M), sing);
dim J, degree J
