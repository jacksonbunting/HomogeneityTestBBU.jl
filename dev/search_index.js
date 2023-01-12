var documenterSearchIndex = {"docs":
[{"location":"#HomogeneityTestBBU.jl","page":"Home","title":"HomogeneityTestBBU.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"HomogeneityTestBBU implements the permutation test from Bugni, Bunting and Ura (2023). By default, the test function takes the data X and a test statistic test_stat_fn and returns the p-value computed as Equation (4), Bugni, Bunting and Ura (2023). The test rejects the null of homogeneity for all significance levels greater than the p-value.","category":"page"},{"location":"","page":"Home","title":"Home","text":"homogeneity_test_fn","category":"page"},{"location":"#HomogeneityTestBBU.homogeneity_test_fn","page":"Home","title":"HomogeneityTestBBU.homogeneity_test_fn","text":"homogeneity_test_fn\n\nThis function implements the homogeneity test of Bugni, Bunting and Ura (2023).\n\nArguments\n\nX::Tuple: The data in the form of tuple X=(S, A)\ntest_stat_fn::Function: A function that takes the arguments S and A and returns a (vector of) test statistics.\nK::Int: An integer indicating the length of the MCMC chain. Defaults to 10,000.\nverbose::Boolean: A Boolean indicating if additional output is to be returned. Defaults to false.  \n\nValues\n\nPvalue: The test's p-value (Bugni, Bunting and Ura (2023), equation (4)).\ntest_stat: The test statistic computed on X (optional).\nPvalue_chain: The p-value computed at each of K steps (optional).\ntest_stat_chain: The test statistic computed at each of K steps (optional).\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"euler_algorithm","category":"page"},{"location":"#HomogeneityTestBBU.euler_algorithm","page":"Home","title":"HomogeneityTestBBU.euler_algorithm","text":"euler_algorithm\n\nThis function applies a Euler algorithm to the First Order Markov chain x\n\nArguments\n\nx::vector\n\n\n\n\n\n","category":"function"}]
}
