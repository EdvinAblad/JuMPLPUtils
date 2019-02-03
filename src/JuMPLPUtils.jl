module JuMPLPUtils
using SparseArrays
using LinearAlgebra
using MathOptInterface
using JuMP
import Combinatorics
include("lp_course_utils.jl")
export get_slack
export get_opt_ranges
export get_opt_range

end # module
