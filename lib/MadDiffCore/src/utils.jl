struct Dummy end
const DUMMY = Dummy()
getindex(::Dummy,key::Int) = 0.
    getindex(::Tuple{Int,Dummy},key::Int) = 0.
    setindex!(::Dummy,val,key) = nothing
setindex!(::Tuple{Int,Dummy},val,key) = nothing

index(e) = e.index
index1(e) = e.index1
index2(e) = e.index2
inner(e) = e.inner
ref(e) = e.ref
ref1(e) = e.ref1
ref2(e) = e.ref2
refval(e) = e.ref.x
brefval(e) = e.bref.x
frefval(e) = e.fref.x
frefval1(e) = e.fref1.x
frefval2(e) = e.fref2.x
refval1(e) = e.ref1.x
refval2(e) = e.ref2.x
refval11(e) = e.ref11.x
refval12(e) = e.ref12.x
refval21(e) = e.ref21.x
refval22(e) = e.ref22.x
setrefval(e,val) = ref(e)[] = val
setrefval1(e,val) = ref1(e)[] = val
setrefval2(e,val) = ref2(e)[] = val
setbrefval(e,val) = e.bref[] = val
addrefval(e,val) = ref(e)[] += val
islower(h) = h.islower



