# HomogeneityTestBBU.jl

HomogeneityTestBBU implements the permutation test from Bugni, Bunting and Ura (2023). By default, the test function takes the data ```X``` and a test statistic ```test_stat_fn``` and returns the p-value computed as Equation (4), Bugni, Bunting and Ura (2023). The test rejects the null of homogeneity for all significance levels greater than the p-value.

```@docs
homogeneity_test_fn
```

```@docs
euler_algorithm
```


## Empirical vignette

The following code replicates the empirical application in Bugni, Bunting and Ura (2023).

```@setup init
using Random, HomogeneityTestBBU
S = [ 13  13  13  13  11   8   8   8   7   7   6   6   9   9   9   8   8   9   9;
 15  20  19  19  20  20  14  15  21  20  20  20  21  21  19  19  19  20  20;
 14  14  14  14  14  15  15  15  15  15  15  15  13  13  12  12  12  12  12;
 12  12  12  12  12  12  12  12  12  12  13  13  13  13  13  13  12  13  13;
 32  36  36  37  38  38  38  36  32  33  33  35  36  35  33  33  33  35  34;
 13  18  16  16  16  17  17  17  14  15  15  14  14  15  13  13  13  13  13;
  7  10  10  10  10   8   8   9  11  11  11  10  10  10   9  10  10  10  10;
 17  17  18  18  18  18  18  16  16  16  16  15  15  15  14  14  14  14  14;
 14  13  13  11  11  11  11  12  12  12  12  11  11  11  10  11  11  11  10;
 23  23  22  24  24  24  24  24  23  22  22  20  20  20  19  20  20  20  20;
  8  10  10  10  10  10  10  10  10   8   8   8   8   8   8   8   9   9  10;
 13  13  13  12  12   7   7   9  12  12  12  13  13  13  12  12  12  12  12;
 14  14  14  14  14  14  14  14  14  14  14  14  14  14  13  12  12  13  13;
 11  11  11  11  11  11  11  11  11  11   9   9   9   9   8   8   8   8   8;
 11  11  11   9   9   9   8   9   9   6   6   6   6   6   5   6   6   6   6;
 17  18  18  18  18  17  18  18  18  18  18  18  18  18  16  16  16  16  17;
 28  25  23  25  25  23  23  23  23  23  23  23  24  24  22  22  22  20  20;
 22  21  21  22  21  21  22  22  22  21  22  21  21  22  20  20  20  20  21;
 19  19  15  17  17  17  17  17  17  17  17  17  18  18  14  14  14  14  15;
 11  12  12  12  11  11  11   9   8   8   8   8   8   5   5   5   5   5   5;
 33  32  27  30  30  30  30  31  31  31  31  33  33  34  30  31  30  31  31;
 12  12  12  11  12  12  12  12  12  12  12  13  12  12  11  11  11  11  12;
 47  46  45  51  48  47  47  44  43  39  39  39  39  39  36  36  36  36  36]

function test_fn(S_data)
    S_support = sort(unique(S_data[:,1:(end-1)]))
    outcome_support = sort(unique(S_data[:,2:end]))
    
    S_cols = size(S_data, 2);

    T_P = 0;
    T_P_star = 0;
    
    for j2 in eachindex(S_support)
        tmp0 = S_data[:,1:(end-1)] .== S_support[j2];
        tmp1 = unique(S_data[:,2:end][tmp0])
        tmp_i = findall(vec(sum(tmp0, dims=2) .> 0))
        for j3 in eachindex(tmp1)
            P_d_tmp = sum((S_data[:,2:end][tmp0] .== tmp1[j3])) / sum(tmp0);
            for j1 in tmp_i
                P_j_d_tmp = sum((S_data[j1,2:end][tmp0[j1,:]] .== tmp1[j3])) / sum(tmp0[j1,:]);
                W_j_d_tmp = sum(tmp0[j1,:]) / P_d_tmp
                T_P = T_P + W_j_d_tmp * (P_j_d_tmp - P_d_tmp)^2
                if P_j_d_tmp > 0
                    T_P_star = T_P_star + 2 *  W_j_d_tmp * P_d_tmp * P_j_d_tmp * log(P_j_d_tmp / P_d_tmp)
                end
            end
        end
    end
    return T_P, T_P_star
end

homogeneity_test_fn(X = (S[:,12:end], []),
    test_stat_fn = test_fn,
    K=3)
Random.seed!(1234)
```

```@example init
Before = S[:,1:11];
After = S[:,12:end];

@time tAfter = homogeneity_test_fn(X = (After, []),
    test_stat_fn = test_fn,
    K=50000)
println(tAfter.Pvalue)

@time tBefore = homogeneity_test_fn(X = (Before, []),
    test_stat_fn = test_fn,
    K=50000)
println(tBefore.Pvalue)
```


## Data and test statistic for empirical vignette
```@example
S = [ 13  13  13  13  11   8   8   8   7   7   6   6   9   9   9   8   8   9   9;
 15  20  19  19  20  20  14  15  21  20  20  20  21  21  19  19  19  20  20;
 14  14  14  14  14  15  15  15  15  15  15  15  13  13  12  12  12  12  12;
 12  12  12  12  12  12  12  12  12  12  13  13  13  13  13  13  12  13  13;
 32  36  36  37  38  38  38  36  32  33  33  35  36  35  33  33  33  35  34;
 13  18  16  16  16  17  17  17  14  15  15  14  14  15  13  13  13  13  13;
  7  10  10  10  10   8   8   9  11  11  11  10  10  10   9  10  10  10  10;
 17  17  18  18  18  18  18  16  16  16  16  15  15  15  14  14  14  14  14;
 14  13  13  11  11  11  11  12  12  12  12  11  11  11  10  11  11  11  10;
 23  23  22  24  24  24  24  24  23  22  22  20  20  20  19  20  20  20  20;
  8  10  10  10  10  10  10  10  10   8   8   8   8   8   8   8   9   9  10;
 13  13  13  12  12   7   7   9  12  12  12  13  13  13  12  12  12  12  12;
 14  14  14  14  14  14  14  14  14  14  14  14  14  14  13  12  12  13  13;
 11  11  11  11  11  11  11  11  11  11   9   9   9   9   8   8   8   8   8;
 11  11  11   9   9   9   8   9   9   6   6   6   6   6   5   6   6   6   6;
 17  18  18  18  18  17  18  18  18  18  18  18  18  18  16  16  16  16  17;
 28  25  23  25  25  23  23  23  23  23  23  23  24  24  22  22  22  20  20;
 22  21  21  22  21  21  22  22  22  21  22  21  21  22  20  20  20  20  21;
 19  19  15  17  17  17  17  17  17  17  17  17  18  18  14  14  14  14  15;
 11  12  12  12  11  11  11   9   8   8   8   8   8   5   5   5   5   5   5;
 33  32  27  30  30  30  30  31  31  31  31  33  33  34  30  31  30  31  31;
 12  12  12  11  12  12  12  12  12  12  12  13  12  12  11  11  11  11  12;
 47  46  45  51  48  47  47  44  43  39  39  39  39  39  36  36  36  36  36];

function test_fn(S_data)
    S_support = sort(unique(S_data[:,1:(end-1)]))
    outcome_support = sort(unique(S_data[:,2:end]))
    
    S_cols = size(S_data, 2);

    T_P = 0;
    T_P_star = 0;
    
    for j2 in eachindex(S_support)
        tmp0 = S_data[:,1:(end-1)] .== S_support[j2];
        tmp1 = unique(S_data[:,2:end][tmp0])
        tmp_i = findall(vec(sum(tmp0, dims=2) .> 0))
        for j3 in eachindex(tmp1)
            P_d_tmp = sum((S_data[:,2:end][tmp0] .== tmp1[j3])) / sum(tmp0);
            for j1 in tmp_i
                P_j_d_tmp = sum((S_data[j1,2:end][tmp0[j1,:]] .== tmp1[j3])) / sum(tmp0[j1,:]);
                W_j_d_tmp = sum(tmp0[j1,:]) / P_d_tmp
                T_P = T_P + W_j_d_tmp * (P_j_d_tmp - P_d_tmp)^2
                if P_j_d_tmp > 0
                    T_P_star = T_P_star + 2 *  W_j_d_tmp * P_d_tmp * P_j_d_tmp * log(P_j_d_tmp / P_d_tmp)
                end
            end
        end
    end
    return T_P, T_P_star
end
```
