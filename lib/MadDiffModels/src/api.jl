"""
    value(x::ModelComponent{V}) where V <: MadDiffCore.Variable
Return the value of variable `x`.
"""
value(x::ModelComponent{V}) where V <: MadDiffCore.Variable = x.parent.x[x.inner.index]

"""
    value(p::ModelComponent{P}) where P <: MadDiffCore.Parameter
Return the value of parameter `p`.
"""
value(p::ModelComponent{P}) where P <: MadDiffCore.Parameter= p.parent.p[p.inner.index]

"""
    setvalue(x::ModelComponent{V},val) where V <: MadDiffCore.Variable
Set the value of variable 'x' to `val`.
"""
setvalue(x::ModelComponent{V},val) where V <: MadDiffCore.Variable = x.parent.x[x.inner.index] = val

"""
    setvalue(p::ModelComponent{P},val) where P <: MadDiffCore.Parameter
Set the value of parameter 'p' to `val`.
"""
setvalue(p::ModelComponent{P},val) where P <: MadDiffCore.Parameter = p.parent.p[p.inner.index] = val

"""
    set_lower_bound(x::ModelComponent{V},val) where V <: MadDiffCore.Variable
Set the lower bound of variable 'x' to `val`.
"""
set_lower_bound(x::ModelComponent{V},val) where V <: MadDiffCore.Variable = x.parent.xl[x.inner.index] = val

"""
    set_upper_bound(x::ModelComponent{V},val) where V <: MadDiffCore.Variable
Set the upper bound of variable 'x' to `val`.
"""
set_upper_bound(x::ModelComponent{V},val) where V <: MadDiffCore.Variable = x.parent.xu[x.inner.index] = val

"""
    lower_bound(x::ModelComponent{V}) where V <: MadDiffCore.Variable
Retrun the lower bound of variable `x`.
"""
lower_bound(x::ModelComponent{V}) where V <: MadDiffCore.Variable = x.parent.xl[x.inner.index]

"""
    upper_bound(x::ModelComponent{V}) where V <: MadDiffCore.Variable
Retrun the upper bound of variable `x`.
"""
upper_bound(x::ModelComponent{V}) where V <: MadDiffCore.Variable = x.parent.xu[x.inner.index]

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
