using JuMP
using Clp
using Test
using JuMPLPUtils
import Random
Random.seed!(1234)
A = rand(40, 10)
c = rand(10)
b = rand(40)
m = Model(with_optimizer(Clp.Optimizer, LogLevel = 0))
J = 1:10
I = 1:40
@variable(m, x[J] >= 0)
@variable(m, b_var[I])
for i in I
    fix(b_var[i],b[i])
end
#minimize the cost
z = @objective(m, Min, sum(c[j]*x[j] for j in J))
#Satisfy the constraints

@constraint(m, constraint[i in I], sum(A[i, j]*x[j] for j in J) >= b_var[i])
@constraint(m, int_constr, .2 <= x[4] + x[5] <= 2)
@constraint(m, eq_constr,  x[1] + x[2] + x[3] + x[4] == 1.0)
optimize!(m)

N = [1,3,4,6,7,8,9]
B = [2,5,10]
x_val = value.(x).data

@test sum(x_val[N]) <= 1e-6
@test abs(sum(x_val[B]) -2.53198) <= 1e-3
@test abs(get_slack(constraint[10])-0.51868) <= 1e-3

s = A*x_val - b
@test sum(abs(s[i] -get_slack(constraint[i])) for i in I) <= 1e-10

c_int, b_int = get_opt_ranges(m)
z_val = value(z)
pi = dual.(constraint).data
for j in J
    dc = .8*max(c_int[x[j]][1],-1.0)
    c[j] += dc
    @objective(m, Min, sum(c[j]*x[j] for j in J))
    optimize!(m)
    @test sum(abs.(x_val - value.(x).data)) <= 1e-10
    c[j] -= dc
    dc = .8*min(c_int[x[j]][2],1.0)
    c[j] += dc
    @objective(m, Min, sum(c[j]*x[j] for j in J))
    optimize!(m)
    @test sum(abs.(x_val - value.(x).data)) <= 1e-10
    c[j] -= dc
end
@objective(m, Min, sum(c[j]*x[j] for j in J))

for i in I
    db = .9*max(b_int[constraint[i]][1],-1.0)
    fix(b_var[i], b[i] + db)
    optimize!(m)
    @test sum(value.(x).data[N]) <= 1e-10
    @test abs(z_val - value(z) + db*pi[i]) <= 1e-10
    db= .9*min(b_int[constraint[i]][2], 1.0)
    fix(b_var[i], b[i] + db)
    optimize!(m)
    @test sum(value.(x).data[N]) <= 1e-10
    @test abs(z_val - value(z) + db*pi[i]) <= 1e-10
    fix(b_var[i], b[i])
end


m = Model(with_optimizer(Clp.Optimizer, LogLevel = 0))
@variable(m, x[J] >= 0)
@variable(m, b_var[I])
for i in I
    fix(b_var[i], b[i])
end
#minimize the cost
z = @objective(m, Max, sum(c[j]*x[j] for j in J))
#Satisfy the constraints
@constraint(m, constraint[i in I], sum(A[i, j]*x[j] for j in J) <= b_var[i])
optimize!(m)

N = [1,2,4,5,6,7,8,9,10]
B = [3]
x_val =value.(x).data

@test sum(x_val[N]) <= 1e-6
@test abs(sum(x_val[B]) - .063688) <= 1e-3
@test abs(get_slack(constraint[10])-.5704) <= 1e-3

s = b- A*x_val
@test sum(abs(s[i] -get_slack(constraint[i])) for i in I) <= 1e-10

c_int, b_int = get_opt_ranges(m)
z_val = value(z)
pi = -dual.(constraint).data
for j in J
    dc = .99*max(c_int[x[j]][1],-1.0)
    c[j] += dc
    @objective(m, Max, sum(c[j]*x[j] for j in J))
    optimize!(m)
    @test sum(abs.(x_val - value.(x).data)) <= 1e-10
    c[j] -= dc
    dc = .99*min(c_int[x[j]][2],1.0)
    c[j] += dc
    @objective(m, Max, sum(c[j]*x[j] for j in J))
    optimize!(m)
    @test sum(abs.(x_val - value.(x).data)) <= 1e-10
    c[j] -= dc
end
@objective(m, Max, sum(c[j]*x[j] for j in J))

for i in I
    db = .99*max(b_int[constraint[i]][1],-1.0)
    fix(b_var[i], b[i] + db)
    optimize!(m)
    @test sum(value.(x).data[N]) <= 1e-10
    @test abs(z_val - value(z) + db*pi[i]) <= 1e-10
    db= .99*min(b_int[constraint[i]][2], 1.0)
    fix(b_var[i], b[i] + db)
    optimize!(m)
    @test sum(value.(x).data[N]) <= 1e-10
    @test abs(z_val - value(z) + db*pi[i]) <= 1e-10
    fix(b_var[i], b[i] )
end
