\ENCODING=ISO-8859-1
\Problem name: IloCplex

Minimize
 obj1: - x1 - 2 x2 - 3 x3 + [ 33 x1 ^2 - 12 x1 * x2 + 22 x2 ^2 - 23 x2 * x3
       + 11 x3 ^2 ] / 2
Subject To
 c1: - x1 + x2 + x3 <= 20
 c2: x1 - 3 x2 + x3 <= 30
 c3:  = 0
 q1: [ 0.5 x1 ^2 + 0.5 x2 ^2 + 0.5 x3 ^2 ] <= 1
Bounds
-1000000 <= x1 <= 1000000
-1000000 <= x2 <= 1000000
-1000000 <= x3 <= 1000000
Generals
 x1  x2 
End
