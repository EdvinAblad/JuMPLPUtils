"""
 Gets the current slack of the constraint, for feasible solution it's always positive.
 For double sided inequalities it's the least slack that is given.
"""
function getslack(constraint::ConstraintRef)
  i = constraint.idx
  row_val = MathProgBase.getconstrsolution(internalmodel(constraint.m))[i]
  row_lb = MathProgBase.getconstrLB(internalmodel(constraint.m))[i]
  row_ub = MathProgBase.getconstrUB(internalmodel(constraint.m))[i]
  row_basis = MathProgBase.getbasis(internalmodel(constraint.m))[2][i]
  if row_basis != :Basic
    return 0.0
  end
  return min(row_ub - row_val, row_val - row_lb)
end

"""
  Computes the perturbation range for costs and rhs, in which the current basis remains optimal.
"""
function getoptranges(model::Model)
  col_basis, row_basis = MathProgBase.getbasis(internalmodel(model))
  cost_int = zeros(model.numCols,2)  # intervals for non-basic costs
  cost_int[col_basis .== :NonbasicAtLower, 1] = -model.redCosts[col_basis .== :NonbasicAtLower]
  cost_int[col_basis .== :NonbasicAtLower, 2] .= Inf

  cost_int[col_basis .== :NonbasicAtUpper, 1] .= -Inf
  cost_int[col_basis .== :NonbasicAtUpper, 2] = model.redCosts[col_basis .== :NonbasicAtUpper]

  if model.objSense == :Max
    cost_int[:,1], cost_int[:,2] = cost_int[:,2], cost_int[:,1]
    cost_int[isinf.(cost_int)] *= -1.0;
  end
  # intervals for basic cost costcoefficient
  #c_N - c_B B^-1 N >= 0 => c_N - pi N >= dc_B B^-1 N
  #
  A = MathProgBase.getconstrmatrix(internalmodel(model))
  Ident = spdiagm(0=>ones(A.m))
  B = hcat(A[:,col_basis .== :Basic], Ident[:,row_basis .== :Basic])
  N = hcat(A[:,col_basis .== :NonbasicAtLower], A[:,col_basis .== :NonbasicAtUpper],
    Ident[:,row_basis .== :NonbasicAtLower], Ident[:,row_basis .== :NonbasicAtUpper])
  y = model.linconstrDuals

  c_red = vcat(model.redCosts[col_basis .== :NonbasicAtLower],
   model.redCosts[col_basis .== :NonbasicAtUpper],
   -y[row_basis .== :NonbasicAtLower], -y[row_basis .== :NonbasicAtUpper])

  cost_int[col_basis .== :Basic, 1] .= -Inf
  cost_int[col_basis .== :Basic, 2] .= Inf
  for i = 1:sum(col_basis .== :Basic)
    j = (1:model.numCols)[col_basis .== :Basic][i]
    dcb = zeros(A.m)
    dcb[i] = 1.0
    A_red = ((B'\dcb)'*N)'
    neg_ind = ((A_red .< -1e-6) .& (c_red .> 1e-8)) .|  ((A_red .> 1e-6) .& (c_red .< -1e-8))
    cost_int[j, 1] = maximum([cost_int[j, 1]; c_red[neg_ind]./A_red[neg_ind]])
    pos_ind = ((A_red .> 1e-6) .& (c_red .> 1e-8)) .|  ((A_red .< -1e-6) .& (c_red .< -1e-8))
    cost_int[j, 2] = minimum([cost_int[j, 2]; c_red[pos_ind]./A_red[pos_ind]])
  end
  #Intervals for rhs: U >= B^-1(b+td) >= L  =>  t B^-1d >= L-B^-1b && t B^-1d <= U -B^-1b
  B = hcat(A[:,col_basis .== :Basic], -Ident[:,row_basis .== :Basic])
  x_B = model.colVal[col_basis .== :Basic]
  x_lb = model.colLower[col_basis .== :Basic]
  x_ub = model.colUpper[col_basis .== :Basic]
  row_val = MathProgBase.getconstrsolution(internalmodel(model))[row_basis .== :Basic]
  row_lb = MathProgBase.getconstrLB(internalmodel(model))[row_basis .== :Basic]
  row_ub = MathProgBase.getconstrUB(internalmodel(model))[row_basis .== :Basic]
  x_lb_d = [x_lb - x_B; row_lb - row_val]
  x_ub_d = [x_ub - x_B; row_ub - row_val]

  rhs_int = zeros(A.m,2)
  rhs_int[:,1] .= -Inf
  rhs_int[:,2] .= Inf
  rhs_int[row_basis .== :Basic, 1] = row_val - row_ub
  rhs_int[row_basis .== :Basic, 2] = row_val - row_lb
  for i = (1:A.m)[row_basis .!= :Basic]
    ei = zeros(A.m)
    ei[i] = 1.0
    db = B\ei
    neg_ind = db .< -1e-7
    rhs_int[i, 1] = maximum([rhs_int[i, 1]; x_ub_d[neg_ind]./db[neg_ind]])
    rhs_int[i, 2] = minimum([rhs_int[i, 2]; x_lb_d[neg_ind]./db[neg_ind]])
    pos_ind = db .> 1e-7
    rhs_int[i, 1] = maximum([rhs_int[i, 1]; x_lb_d[pos_ind]./db[pos_ind]])
    rhs_int[i, 2] = minimum([rhs_int[i, 2]; x_ub_d[pos_ind]./db[pos_ind]])
  end
  return cost_int, rhs_int;
end
"""
  Computes the perturbation range for var, in which the current basis remains optimal.
"""
function getoptrange(var::Variable)
  dc, db = getoptranges(var.m)
  return dc[var.col,:]
end

"""
  Computes the perturbation range for constraint, in which the current basis remains optimal.
"""
function getoptrange(constraint::ConstraintRef)
  dc, db = getoptranges(constraint.m)
  return db[constraint.idx,:]
end
