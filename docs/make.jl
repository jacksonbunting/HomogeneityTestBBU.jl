push!(LOAD_PATH,"../src/")
using HomogeneityTestBBU, Documenter
makedocs(
         sitename = "HomogeneityTestBBU.jl"
         , 
        #  modules  = [HomogeneityTestBBU],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(
    repo="github.com/bunting-econ/HomogeneityTestBBU.jl.git",
)
