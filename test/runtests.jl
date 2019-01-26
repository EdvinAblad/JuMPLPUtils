using JuMP
using Clp
using Test
using JuMPLPUtils
import Random
Random.seed!(1234)
A = rand(40, 10)
c = rand(10)
b = rand(40)
m = Model()
J = 1:10
I = 1:40
@variable(m, x[J] >= 0)

#minimize the cost
z = @objective(m, Min, sum(c[j]*x[j] for j in J))
#Satisfy the constraints
@constraint(m, constraint[i in I], sum(A[i, j]*x[j] for j in J) >= b[i])
setsolver(m, ClpSolver(LogLevel = 0))
solve(m)

N = [1,3,6,7,8,9]
B = [2,4,5,10]
x_val =getvalue(x).innerArray

@test sum(x_val[N]) <= 1e-6
@test abs(sum(x_val[B]) -3.179) <= 1e-3
@test abs(getslack(constraint[10])-.775) <= 1e-3

s = A*x_val - b
@test sum(abs(s[i] -getslack(constraint[i])) for i in I) <= 1e-10

c_int, b_int = getoptranges(m)
z_val = getvalue(z)
pi = getdual(constraint)
for j in J
    dc = .8*max(c_int[j,1],-1.0)
    c[j] += dc
    @objective(m, Min, sum(c[j]*x[j] for j in J))
    solve(m)
    @test sum(abs(x_val[j] - getvalue(x)[j]) for j in J) <= 1e-10
    c[j] -= dc
    dc = .8*min(c_int[j,2],1.0)
    c[j] += dc
    @objective(m, Min, sum(c[j]*x[j] for j in J))
    solve(m)
    @test sum(abs(x_val[j] - getvalue(x)[j]) for j in J) <= 1e-10
    c[j] -= dc
end
@objective(m, Min, sum(c[j]*x[j] for j in J))

for i in I
    db = .8*max(b_int[i,1],-1.0)
    JuMP.setRHS(constraint[i], b[i] + db)
    solve(m)
    # @test sum(getvalue(x).innerArray[N]) <= 1e-10
    @test abs(z_val - getvalue(z) + db*pi[i]) <= 1e-10
    db= .8*min(b_int[i,2], 1.0)
    JuMP.setRHS(constraint[i], b[i] + db)
    solve(m)
    # @test sum(getvalue(x).innerArray[N]) <= 1e-10
    @test abs(z_val - getvalue(z) + db*pi[i]) <= 1e-10
    JuMP.setRHS(constraint[i], b[i])
end


m = Model()
@variable(m, x[J] >= 0)
#minimize the cost
z = @objective(m, Max, sum(c[j]*x[j] for j in J))
#Satisfy the constraints
@constraint(m, constraint[i in I], sum(A[i, j]*x[j] for j in J) <= b[i])
setsolver(m, ClpSolver(LogLevel = 0))
solve(m)

N = [1,2,4,5,6,7,8,9,10]
B = [3]
x_val =getvalue(x).innerArray

@test sum(x_val[N]) <= 1e-6
@test abs(sum(x_val[B]) - .063688) <= 1e-3
@test abs(getslack(constraint[10])-.5704) <= 1e-3

s = b- A*x_val
@test sum(abs(s[i] -getslack(constraint[i])) for i in I) <= 1e-10

c_int, b_int = getoptranges(m)
z_val = getvalue(z)
pi = getdual(constraint)
for j in J
    dc = .8*max(c_int[j,1],-1.0)
    c[j] += dc
    @objective(m, Max, sum(c[j]*x[j] for j in J))
    solve(m)
    @test sum(abs(x_val[j] - getvalue(x)[j]) for j in J) <= 1e-10
    c[j] -= dc
    dc = .8*min(c_int[j,2],1.0)
    c[j] += dc
    @objective(m, Max, sum(c[j]*x[j] for j in J))
    solve(m)
    @test sum(abs(x_val[j] - getvalue(x)[j]) for j in J) <= 1e-10
    c[j] -= dc
end
@objective(m, Max, sum(c[j]*x[j] for j in J))

for i in I
    db = .8*max(b_int[i,1],-1.0)
    JuMP.setRHS(constraint[i], b[i] + db)
    solve(m)
    # @test sum(getvalue(x).innerArray[N]) <= 1e-10
    @test abs(z_val - getvalue(z) + db*pi[i]) <= 1e-10
    db= .8*min(b_int[i,2], 1.0)
    JuMP.setRHS(constraint[i], b[i] + db)
    solve(m)
    # @test sum(getvalue(x).innerArray[N]) <= 1e-10
    @test abs(z_val - getvalue(z) + db*pi[i]) <= 1e-10
    JuMP.setRHS(constraint[i], b[i])
end
