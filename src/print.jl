string(p::PrintVariable) = p.str
raw(p::PrintVariable) = p.str
raw(p::PrintTerm) = p.str
string(p::PrintTerm) = paren(p.str)
string(p::PrintSource) = p.str
paren(str) = "(" * str * ")"
getindex(p::PrintSource,n) = PrintVariable(string(p)*"[$n]")

function string(c::Constraint) 
    gl = parent(c).gl[index(c)]
    gu = parent(c).gu[index(c)]
    str = string(parent(c).cons[index(c)])
    return gl==gu ? str * " == $gl" : (gl > -Inf ? "$gl <= " : "") * str * (gl < Inf ? " <= $gu" : "")
end
print(io::IO,e::Constraint) = print(io, string(e))
show(io::IO,::MIME"text/plain",e::Constraint) = print(io,e)


string(m::Model) = "NLP model with $(m.n) variables and $(m.m) constraints"
print(io::IO,e::Model) = print(io, string(e))
show(io::IO,::MIME"text/plain",e::Model) = print(io,e)

string(e::Variable) = string(func(e)(PrintSource(parent(e))...))
string(e::Parameter) = string(func(e)(PrintSource(parent(e))...))
string(e::Term) = raw(func(e)(PrintSource(parent(e))...))
print(io::IO,e::Expression) = print(io, string(e))
string(e::Source) = e.str
print(io::IO,e::Source) = print(io, string(e))

show(io::IO,::MIME"text/plain",e::Expression) = print(io,e)
show(io::IO,::MIME"text/plain",e::Source) = print(io,e)


PrintSource(::Nothing) = (PrintSource(DEFAULT_VAR_STRING),PrintSource(DEFAULT_PAR_STRING))
PrintSource(e::Source{Variable}) = (PrintSource(e.str),nothing)
PrintSource(e::Source{Parameter}) = (nothing,PrintSource(e.str))
PrintSource(e::Tuple) = (PrintSource(e[1].str),PrintSource(e[2].str))
PrintSource(m::Model) = (m.vars,m.pars)

