"""
    homogeneity_test_fn

This function implements the homogeneity test of Bugni, Bunting and Ura (2023).


## Arguments
- `X::Tuple`: The data in the form of tuple X=(S, A)
- `test_stat_fn::Function`: A function that takes the arguments S and A and returns a (vector of) test statistics.
- `K::Int`: An integer indicating the length of the MCMC chain. Defaults to 10,000.
- `verbose::Boolean`: A Boolean indicating if additional output is to be returned. Defaults to false.  

## Values
- `Pvalue`: The test's p-value (Bugni, Bunting and Ura (2023), equation (4)).
- `test_stat`: The test statistic computed on X (optional).
- `Pvalue_chain`: The p-value computed at each of K steps (optional).
- `test_stat_chain`: The test statistic computed at each of K steps (optional).
"""
function homogeneity_test_fn(;
    X::Tuple, ## The data
    test_stat_fn::Function,  ## The test statistic
    K::Int=10000, ## The length of the MCMC chain to generate
    verbose::Bool=false 
    )

    S_data, A_data = X
    if !all(size(S_data) .== size(A_data) )
        println("***S and A are not the same dimension***")
    end
    
    if !(K>0)
        println("***Choose K > 0***")
    end

    test_stat = test_stat_fn(S_data, A_data)
    numb_statistics = length(test_stat)
    
    local_n, local_t = size(S_data);

    SS_list =    hcat(S_data[:,1],S_data[:,2])     
    for t = 2:(local_t-1)
        SS_list = vcat(SS_list, hcat(S_data[:,t],S_data[:,t+1]) );
    end
    SS_list = vcat(SS_list, hcat(S_data[:,end], zeros(Int, local_n)) )
    
    S_data_permuted = copy(S_data);juli
    A_data_permuted = copy(A_data);

    perm_TestStat = fill(NaN, K, numb_statistics);
    for k = 1:K
        rand_markets = rand(1:local_n, 2);
        if rand_markets[1]!=rand_markets[2]
            S_data_mixed = [S_data_permuted[rand_markets[1],:] ;0 ;
                                        S_data_permuted[rand_markets[2],:];0];
            S_data_mixed  = euler_algorithm(S_data_mixed);
            deleteat!(S_data_mixed, [local_t+1 ;length(S_data_mixed)])
            S_data_permuted[rand_markets[1],:] = S_data_mixed[1:local_t];
            S_data_permuted[rand_markets[2],:] = S_data_mixed[local_t.+(1:local_t)];
        end

        for i in 1:local_n
            S_data_permuted[i,:] = euler_algorithm(S_data_permuted[i,:]);
        end

        
        for ind in axes(unique(SS_list, dims=1), 1)
            ss_ind = unique(SS_list, dims=1)[ind, :]'
            ss_ind = all(SS_list .== ss_ind, dims= 2)
            ss_ind = reshape(ss_ind, size(S_data))
            A_data_permuted[ss_ind] = shuffle(A_data[ss_ind])
        end
        
        perm_TestStat[k,:] .= test_stat_fn(S_data_permuted, A_data_permuted);
    end

    pvalue = fill(NaN, numb_statistics);
    pvalue_eachk = fill(NaN, K, numb_statistics)

    for ind in 1:numb_statistics
        pvalue[ind] = sum(perm_TestStat[:,ind] .>= test_stat[ind])/K; 
        if verbose
            for ind2 in 1:K
                pvalue_eachk[ind2, ind] = sum(perm_TestStat[1:ind2, ind] .>= test_stat[ind]) / ind2
            end
        end
    end


    if verbose
        return (Pvalue = pvalue, 
                test_stat = test_stat, 
                Pvalue_chain = pvalue_eachk,
                test_stat_chain = perm_TestStat)
    else
        return Pvalue = pvalue
    end
end