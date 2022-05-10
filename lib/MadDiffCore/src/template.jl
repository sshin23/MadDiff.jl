Field1(inner,::Vector{LagrangianEntry{HessianNull}}) = inner

Hessian11a(::HessianNull,::MyRef{Float64}) = HESSIAN_NULL
Hessian22a(::HessianNull,::HessianNull,::MyRef{Float64},::MyRef{Float64}) = HESSIAN_NULL
Hessian22a(::HessianNull,h2,::MyRef{Float64},ref2::MyRef{Float64}) = Hessian11a(h2,ref2)
Hessian22a(h1,::HessianNull,ref1::MyRef{Float64},::MyRef{Float64}) = Hessian11a(h1,ref1)

# Hessian11(::FZero1,h1,h11,::MyRef{Float64},ref::MyRef{Float64}) = Hessian11a(h1,ref)
# Hessian111(ddf22::FZero2,e,d,indexer) = Hessian11a(Hessian(e.e2,d.d1,indexer),ref(d))
# Hessian112(ddf11::FZero2,e,d,indexer) = Hessian11a(Hessian(e.e1,d.d1,indexer),ref(d))
