"""
    value(x::ModelVariable{T,RT}) where {T,RT}
Return the value of variable `x`.
"""
value(x::ModelVariable{T,RT}) where {T,RT} = x.parent.x[x.index]

"""
    value(p::ModelParameter{T,RT}) where {T,RT}
Return the value of parameter `p`.
"""
value(p::ModelParameter{T,RT}) where {T,RT}= p.parent.p[p.index]

"""
    setvalue(x::ModelVariable{T,RT},val) where {T,RT}
Set the value of variable 'x' to `val`.
"""
setvalue(x::ModelVariable{T,RT},val) where {T,RT} = x.parent.x[x.index] = val

"""
    setvalue(p::ModelParameter{T,RT},val) where {T,RT}
Set the value of parameter 'p' to `val`.
"""
setvalue(p::ModelParameter{T,RT},val) where {T,RT} = p.parent.p[p.index] = val

"""
    set_lower_bound(x::ModelVariable{T,RT},val) where {T,RT}
Set the lower bound of variable 'x' to `val`.
"""
set_lower_bound(x::ModelVariable{T,RT},val) where {T,RT} = x.parent.xl[x.index] = val

"""
    set_upper_bound(x::ModelVariable{T,RT},val) where {T,RT}
Set the upper bound of variable 'x' to `val`.
"""
set_upper_bound(x::ModelVariable{T,RT},val) where {T,RT} = x.parent.xu[x.index] = val

"""
    lower_bound(x::ModelVariable{T,RT}) where {T,RT}
Retrun the lower bound of variable `x`.
"""
lower_bound(x::ModelVariable{T,RT}) where {T,RT} = x.parent.xl[x.index]

"""
    upper_bound(x::ModelVariable{T,RT}) where {T,RT}
Retrun the upper bound of variable `x`.
"""
upper_bound(x::ModelVariable{T,RT}) where {T,RT} = x.parent.xu[x.index]

"""
    set_lower_bound(c::Constraint,val)
Set the lower bound of constraint `c` to `val`.
"""
set_lower_bound(c::Constraint,val) = c.parent.gl[c.index] = val

"""
    set_upper_bound(c::Constraint,val)
Set the upper bound of constraint `c` to `val`.
"""
set_upper_bound(c::Constraint,val) = c.parent.gu[c.index] = val

"""
    lower_bound(c::Constraint)
Retrun the lower bound of constraint `c`.
"""
lower_bound(c::Constraint) = c.parent.gl[c.index]

"""
    upper_bound(c::Constraint)
Retrun the upper bound of constraint `c`.
"""
upper_bound(c::Constraint) = c.parent.gu[c.index]

"""
    dual(c::Constraint)
Retrun the dual of constraint `c`.
"""
dual(c::Constraint) = c.parent.l[c.index]

"""
    dual(c::Constraint)
Retrun the objective value of MadDiffModel `m`.
"""

for field in fieldnames(NLPModels.NLPModelMeta)
    meth = Symbol("get_", field)
    @eval begin
        import NLPModels: $meth
        """
            $($meth)(m::MadDiffModel)
        Return the value $($(QuoteNode(field))) from MadDiffModel.
        """
        $meth
    end
end
