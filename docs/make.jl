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
    repo="github.com/jacksonbunting/HomogeneityTestBBU.jl.git",
)
