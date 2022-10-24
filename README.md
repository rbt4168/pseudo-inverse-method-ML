# pseudo-inverse-method-ML

Written in C17.

Theory is finding " Moore–Penrose inverse of A " and times it to b to slove x\*.
In liner system : Ax\* = b.
x\* = ( transposA \* A )^-1 \* transposA \* b
Ref : [https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse](https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse "https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse")

And the best trainning result.
```c
[ 7783 / 10000 ] validation , accuracy = 77.830000 %
> Judge matrix :
        0       1       2       3       4       5       6       7       8       9
0       938     2       2       1       0       13      15      1       7       1
1       4       171     137     184     20      42      34      34      466     43
2       17      6       875     15      14      0       41      18      41      5
3       6       1       42      851     2       40      8       20      24      16
4       0       5       11      2       859     10      15      2       18      60
5       20      1       4       32      7       731     24      10      51      12
6       22      2       11      0       11      24      882     0       6       0
7       5       2       22      15      18      3       4       894     3       62
8       18      6       17      24      24      77      16      12      758     22
9       21      4       6       14      46      6       1       73      14      824
```

Data get from [http://yann.lecun.com/exdb/mnist/](http://yann.lecun.com/exdb/mnist/ "http://yann.lecun.com/exdb/mnist/")

Just for fun.

rbt4168 
-- 2022/10/24.
