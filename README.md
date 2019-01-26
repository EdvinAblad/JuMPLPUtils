# JuMPLPUtils
Some minor functions to study the sensitivity of an LP solution, e.g., optimallity range for rhs and cost perturbations

The functions included are:
s = getslack(constraint::ConstraintRef) -> gets the slack of the given constraint

c_int, rhs_int = getoptranges(model::Model) -> computes perturbations of c and b for which the given solution remain optimal

c_int  = getoptrange(var::Variable) -> wraps getoptranges and only outputs the interval for var

rhs_int = getoptrange(constraint::ConstraintRef) -> wraps getoptranges and only outputs the interval for constraint

