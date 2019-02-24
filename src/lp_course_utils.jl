const MOI = MathOptInterface
"""
 Gets the current slack of the constraint, for feasible solution it's always positive.
 For double sided inequalities it's the least slack that is given.
"""
function get_slack(constraint::ConstraintRef)::Float64
  constr_obj = constraint_object(constraint)
  row_val = value(constr_obj.func)
  row_range = MOI.Interval(constr_obj.set)
  return max(0.0, min(row_range.upper - row_val, row_val - row_range.lower))
end

function get_slack_lb(constraint::ConstraintRef)::Float64
  constr_obj = constraint_object(constraint)
  row_val = value(constr_obj.func)
  row_range = MOI.Interval(constr_obj.set)
  return row_val - row_range.lower
end

function get_slack_ub(constraint::ConstraintRef)::Float64
  constr_obj = constraint_object(constraint)
  row_val = value(constr_obj.func)
  row_range = MOI.Interval(constr_obj.set)
  return row_range.upper - row_val
end

function get_var_lb(var::VariableRef)
  lb = -Inf
  if has_lower_bound(var)
    lb = lower_bound(var)
  end
  if is_fixed(var)
    lb = fix_value(var)
  end
  return lb
end

function get_var_ub(var::VariableRef)
  ub = Inf
  if has_upper_bound(var)
    ub = upper_bound(var)
  end
  if is_fixed(var)
    ub = fix_value(var)
  end
  return ub
end

"""
Builds the constraint matrix in std form, rhs, and  the upper bound on the slack varaibles
"""
function get_std_matrix(model::Model)
  vars = all_variables(model)
  eq_constr = all_constraints(model, AffExpr, MOI.EqualTo{Float64})
  lq_constr = all_constraints(model, AffExpr, MOI.LessThan{Float64})
  gr_constr = all_constraints(model, AffExpr, MOI.GreaterThan{Float64})
  int_constr = all_constraints(model, AffExpr, MOI.Interval{Float64})
  # Construct A matrix
  col_j = Array{Int64,1}()
  row_i = Array{Int64,1}()
  coeffs = Array{Float64,1}()
  b = Array{Float64,1}()
  slack_upper = Array{Float64,1}()
  col_ofset = length(vars)

  for (i, constr) in pairs(eq_constr)
    constr_obj = constraint_object(constr)
    constr_terms = constr_obj.func.terms
    for (j, var) in pairs(vars)
      if haskey(constr_terms, var)
        push!(col_j, j)
        push!(row_i, i)
        push!(coeffs, constr_terms[var])
      end
    end
    push!(b, constr_obj.set.value)
    push!(col_j, col_ofset + i)
    push!(row_i, i)
    push!(coeffs, 1.0)
    push!(slack_upper, 1e-6)
  end

  row_ofset = length(eq_constr)
  col_ofset += length(eq_constr)

  for (i, constr) in pairs(lq_constr)
    constr_obj = constraint_object(constr)
    constr_terms = constr_obj.func.terms
    for (j, var) in pairs(vars)
      if haskey(constr_terms, var)
        push!(col_j, j)
        push!(row_i, row_ofset + i)
        push!(coeffs, constr_terms[var])
      end
    end
    push!(b, constr_obj.set.upper)
    push!(col_j, col_ofset + i)
    push!(row_i, row_ofset + i)
    push!(coeffs, 1.0)
    push!(slack_upper, Inf)
  end
  row_ofset += length(lq_constr)
  col_ofset += length(lq_constr)

  for (i, constr) in pairs(gr_constr)
    constr_obj = constraint_object(constr)
    constr_terms = constr_obj.func.terms
    for (j, var) in pairs(vars)
      if haskey(constr_terms, var)
        push!(col_j, j)
        push!(row_i, row_ofset + i)
        push!(coeffs, constr_terms[var])
      end
    end
    push!(b, constr_obj.set.lower)
    push!(col_j, col_ofset + i)
    push!(row_i, row_ofset + i)
    push!(coeffs, -1.0)
    push!(slack_upper, Inf)
  end
  row_ofset += length(gr_constr)
  col_ofset += length(gr_constr)

  for (i, constr) in pairs(int_constr)
    constr_obj = constraint_object(constr)
    constr_terms = constr_obj.func.terms
    for (j, var) in pairs(vars)
      if haskey(constr_terms, var)
        push!(col_j, j)
        push!(row_i, row_ofset + i)
        push!(coeffs, constr_terms[var])
      end
    end
    push!(b, constr_obj.set.lower)
    push!(col_j, col_ofset + i)
    push!(row_i, row_ofset + i)
    push!(coeffs, -1.0)
    push!(slack_upper, constr_obj.set.upper - constr_obj.set.lower)
  end


  A = sparse(row_i,col_j,coeffs)

  fixed_vars = (1:length(vars))[is_fixed.(vars)]
  b -= A[:,fixed_vars]*value.(vars)[fixed_vars]

  A = A[:, setdiff(1:size(A)[2],fixed_vars)]

  objective = objective_function(model).terms
  c = zeros(size(A)[2])
  for (j, var) in pairs(vars)
    if haskey(objective, var)
      c[j] = objective[var]
    end
  end
  return  A, b, c, slack_upper
end
"""
Get dual solution
"""
function get_std_dual(model::Model)::Array{Float64,1}
  eq_constr = all_constraints(model, AffExpr, MOI.EqualTo{Float64})
  lq_constr = all_constraints(model, AffExpr, MOI.LessThan{Float64})
  gr_constr = all_constraints(model, AffExpr, MOI.GreaterThan{Float64})
  int_constr = all_constraints(model, AffExpr, MOI.Interval{Float64})
  return ((objective_sense(model) == MOI.MAX_SENSE) ? -1.0 : 1.0 ) .* dual.([eq_constr;lq_constr;gr_constr;int_constr])
end
"""
Get basic variables
"""
function get_col_basis(model::Model, vars::Vector{VariableRef})::Vector{Int64}
  eq_constr = all_constraints(model, AffExpr, MOI.EqualTo{Float64})
  lq_constr = all_constraints(model, AffExpr, MOI.LessThan{Float64})
  gr_constr = all_constraints(model, AffExpr, MOI.GreaterThan{Float64})
  int_constr = all_constraints(model, AffExpr, MOI.Interval{Float64})
  col_basis = Array{Int64,1}()
  try
    for (j, var) in pairs(vars)
      if MOI.get(backend(model), MOI.VariableBasisStatus(), var.index) == MOI.BASIC
        push!(col_basis, j)
      end
    end
    off = length(vars) + length(eq_constr)
    for (j, con) in pairs(lq_constr)
      if MOI.get(backend(model), MOI.ConstraintBasisStatus(), con.index) == MOI.BASIC
        push!(col_basis, j + off)
      end
    end
    off += length(lq_constr)
    for (j, con) in pairs(gr_constr)
      if MOI.get(backend(model), MOI.ConstraintBasisStatus(), con.index) == MOI.BASIC
        push!(col_basis, j + off)
      end
    end
    off += length(gr_constr)
    for (j, con) in pairs(int_constr)
      if MOI.get(backend(model), MOI.ConstraintBasisStatus(), con.index) == MOI.BASIC
        push!(col_basis, j + off)
      end
    end
  catch
    @warn("VariableBasisStatus not supported by solver, using brute force method")
    A, b, c, slack_upper = get_std_matrix(model)
    n_vars = length(vars)
    n_slack = length(slack_upper)
    n = n_slack + n_vars
    m = length(b)
    #Identify the basic and non-basic variables
    x_val = value.(vars)
    s_val = A[:,(n_vars+1):n]*(b - A[:,1:n_vars]*x_val)
    x_lb =  get_var_lb.(vars)
    x_ub =  get_var_ub.(vars)
    bound_dist = min.(x_val - x_lb, x_ub - x_val)
    slack_dist = min.(s_val, slack_upper- s_val)
    x_basis = (1:n_vars)[bound_dist .> 1e-6]
    col_basis = [x_basis; (n_vars+1:n)[slack_dist .> 1e-6]]
    if length(col_basis) < m
      y = get_std_dual(model)
      c_red = c - A'*y
      B = A[:,col_basis]
      basis_cand = (1:n)[abs.(c_red) .< 1e-6]
      new_basis_cand = setdiff(basis_cand, col_basis)
      if length(col_basis) + length(new_basis_cand) < m
        @warn "Degenerate primal solution, sensitiviy might be too restrictive, try perturbing rhs to get a non-degenerate solution"
        return col_basis
      end
      m_missing = m - length(col_basis)
      miss_range = (m - m_missing + 1):m
      alt_new_basis = collect(Combinatorics.combinations(new_basis_cand, m_missing))
      append!(col_basis, alt_new_basis[1])
      B = A[:, col_basis]
      max_det = abs(det(B_cand))
      for i in 2:min(length(alt_new_basis), 100)
          B[:,miss_range] = A[:,alt_new_basis[i]]
          cur_det = abs(det(B_cand))
          if cur_det > max_det
            max_det = cur_det
            col_basis[miss_range] = alt_new_basis[i]
          end
      end
    end
  end
  return col_basis
end

"""
  Computes the perturbation range for costs and rhs, in which the current basis remains optimal.
"""
function get_opt_ranges(model::Model)
  A, b, c, slack_upper = get_std_matrix(model)
  vars = all_variables(model)
  vars = vars[.!is_fixed.(vars)]
  n_vars = length(vars)
  n_slack = length(slack_upper)
  n = n_slack + n_vars
  m = length(b)
  #Identify the basic and non-basic variables
  x_val = value.(vars)
  s_val = A[:,(n_vars+1):n]*(b - A[:,1:n_vars]*x_val)
  col_basis = get_col_basis(model, vars)
  x_lb =  get_var_lb.(vars)
  x_ub =  get_var_ub.(vars)
  bound_dist = min.(x_val - x_lb, x_ub - x_val)
  x_basis = (1:n_vars)[bound_dist .> 1e-6]
  x_at_lower = (1:n_vars)[x_val - x_lb .<= 1e-6]
  x_at_upper = (1:n_vars)[x_ub - x_val .<= 1e-6]

  y = get_std_dual(model)
  c_red = c - A'*y
  # do the sensitiviy analysis
  cost_int = zeros(n_vars, 2)  # intervals for non-basic costs
  cost_int[x_at_lower, 1] = -c_red[x_at_lower]
  cost_int[x_at_lower, 2] .= Inf

  cost_int[x_at_upper, 1] .= -Inf
  cost_int[x_at_upper, 2] = c_red[x_at_upper]

  if objective_sense(model) == MOI.MAX_SENSE
    cost_int[:,1], cost_int[:,2] = cost_int[:,2], cost_int[:,1]
    cost_int[isinf.(cost_int)] *= -1.0;
  end
  B = A[:,col_basis]

  N = A[:,setdiff(1:n,col_basis)]
  # intervals for basic cost costcoefficient
  #c_N - c_B B^-1 N >= 0 => c_red = c_N - pi N >= dc_B B^-1 N = t*a_red[i]
  dcb = Matrix{Float64}(LinearAlgebra.I,m,m)[:,1:length(x_basis)]
  A_red = (B'\dcb)'*N
  c_red = c_red[setdiff(1:n,col_basis)]
  cost_int[x_basis, 1] .= -Inf
  cost_int[x_basis, 2] .= Inf
  for i = 1:length(x_basis)
    j = x_basis[i]
    neg_ind = ((A_red[i,:] .< -1e-6) .& (c_red .> 1e-8)) .|  ((A_red[i,:] .> 1e-6) .& (c_red .< -1e-8))
    cost_int[j, 1] = maximum([cost_int[j, 1]; c_red[neg_ind]./A_red[i,neg_ind]])
    pos_ind = ((A_red[i,:] .> 1e-6) .& (c_red .> 1e-8)) .|  ((A_red[i,:] .< -1e-6) .& (c_red .< -1e-8))
    cost_int[j, 2] = minimum([cost_int[j, 2]; c_red[pos_ind]./A_red[i,pos_ind]])
  end
  #Intervals for rhs: U >= B^-1(b+td) >= L  =>  t B^-1d >= L-B^-1b && t B^-1d <= U -B^-1b
  x_B = [x_val;s_val][col_basis]
  x_lb = [x_lb;zeros(length(s_val))][col_basis]
  x_ub = [x_ub;slack_upper][col_basis]
  x_lb_d = x_lb - x_B
  x_ub_d = x_ub - x_B

  rhs_int = zeros(A.m,2)
  rhs_int[:,1] .= -Inf
  rhs_int[:,2] .= Inf
  row_basis = col_basis[col_basis .> n_vars] .- n_vars
  eq_constr = all_constraints(model, AffExpr, MOI.EqualTo{Float64})
  lq_constr = all_constraints(model, AffExpr, MOI.LessThan{Float64})
  gr_constr = all_constraints(model, AffExpr, MOI.GreaterThan{Float64})
  int_constr = all_constraints(model, AffExpr, MOI.Interval{Float64})
  all_constr = [eq_constr;lq_constr;gr_constr;int_constr]
  basis_constraints = all_constr[row_basis]
  rhs_int[row_basis, 1] = -get_slack_ub.(basis_constraints)
  rhs_int[row_basis, 2] = get_slack_lb.(basis_constraints)

  non_row_basis = setdiff(1:m, row_basis)
  ei = Matrix{Float64}(LinearAlgebra.I,m,m)[:,non_row_basis]
  db = B\ei
  for ii = 1:length(non_row_basis)
    i = non_row_basis[ii]
    neg_ind = db[:,ii] .< -1e-7
    rhs_int[i, 1] = maximum([rhs_int[i, 1]; x_ub_d[neg_ind]./db[neg_ind,ii]])
    rhs_int[i, 2] = minimum([rhs_int[i, 2]; x_lb_d[neg_ind]./db[neg_ind,ii]])
    pos_ind = db[:,ii] .> 1e-7
    rhs_int[i, 1] = maximum([rhs_int[i, 1]; x_lb_d[pos_ind]./db[pos_ind,ii]])
    rhs_int[i, 2] = minimum([rhs_int[i, 2]; x_ub_d[pos_ind]./db[pos_ind,ii]])
  end

  dict_cost_int = Dict(vars[j]=>(cost_int[j,1],cost_int[j,2]) for j in 1:n_vars)
  dict_rhs_int = Dict(all_constr[i]=>(rhs_int[i,1],rhs_int[i,2]) for i in 1:m)

  return dict_cost_int, dict_rhs_int
end
"""
  Computes the perturbation range for var, in which the current basis remains optimal.
"""
function getoptrange(var::VariableRef)
  dc, db = getoptranges(var.model)
  return dc[var.index,:]
end

"""
  Computes the perturbation range for constraint, in which the current basis remains optimal.
"""
function getoptrange(constraint::ConstraintRef)
  dc, db = getoptranges(constraint.model)
  return db[constraint.index,:]
end
