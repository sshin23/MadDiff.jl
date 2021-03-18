abstract type PrintExpression end

struct PrintSource
    str::String
end

struct PrintVariable <: PrintExpression
    str::String
end

struct PrintTerm <: PrintExpression
    str::String
end

-(p::PrintExpression) = PrintTerm("-"*string(p))
+(p1::PrintExpression,p2::PrintExpression) = PrintTerm(raw(p1)*" + "*raw(p2))
-(p1::PrintExpression,p2::PrintExpression) = PrintTerm(raw(p1)*" - "*string(p2))
+(p1::T,p2::PrintExpression) where T <: Real = PrintTerm(string(p1)*" + "*raw(p2))
+(p1::PrintExpression,p2::T) where T <: Real = PrintTerm(raw(p1)*" + "*string(p2))
-(p1::T,p2::PrintExpression) where T <: Real = PrintTerm(string(p1)*" - "*raw(p2))
-(p1::PrintExpression,p2::T) where T <: Real = PrintTerm(raw(p1)*" - "*string(p2))

*(p1::PrintExpression,p2::PrintExpression) = PrintVariable(string(p1)*"*"*string(p2))
/(p1::PrintExpression,p2::PrintExpression) = PrintVariable(string(p1)*"/"*string(p2))
^(p1::PrintExpression,p2::PrintExpression) = PrintVariable(string(p1)*"^"*string(p2))

for o in [:*,:/,:^]
    @eval begin
        $o(p1::T,p2::PrintExpression) where T <: Real = PrintVariable(""*string(p1)*string($o)*string(p2))
        $o(p1::PrintExpression,p2::T) where T <: Real = PrintVariable(string(p1)*string($o)*string(p2)*"")
    end
end

for (M,f,nargs) in diffrules
    if nargs == 1
        @eval begin
            $f(p::PrintExpression) = PrintVariable(string($f)*"("*raw(p)*")")
        end
    elseif nargs == 2
        @eval begin
            $f(a::T,p::PrintExpression) where T <: Real = PrintVariable(string($f)*"("*string(a)*" , "*raw(p)*")")
            $f(p::PrintExpression,a::T) where T <: Real = PrintVariable(string($f)*"("*raw(p)*" , "*string(a)*")")
            $f(p1::PrintExpression,p2::PrintExpression) = PrintVariable(string($f)*"("*raw(p1)*" , "*raw(p2)*")")
        end
    end
end
