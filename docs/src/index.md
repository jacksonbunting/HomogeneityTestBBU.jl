# HomogeneityTestBBU.jl

HomogeneityTestBBU implements the permutation test from Bugni, Bunting and Ura (2023). By default, the test function takes the data ```X``` and a test statistic ```test_stat_fn``` and returns the p-value computed as Equation (4), Bugni, Bunting and Ura (2023). The test rejects the null of homogeneity for all significance levels greater than the p-value.

```@docs
homogeneity_test_fn
```

```@docs
euler_algorithm
```