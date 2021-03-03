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
for T in Reals
    @eval begin
        +(p1::$T,p2::PrintExpression) = PrintTerm(string(p1)*" + "*raw(p2))
        +(p1::PrintExpression,p2::$T) = PrintTerm(raw(p1)*" + "*string(p2))
        -(p1::$T,p2::PrintExpression) = PrintTerm(string(p1)*" - "*raw(p2))
        -(p1::PrintExpression,p2::$T) = PrintTerm(raw(p1)*" - "*string(p2))
    end
end

*(p1::PrintExpression,p2::PrintExpression) = PrintVariable(string(p1)*"*"*string(p2))
/(p1::PrintExpression,p2::PrintExpression) = PrintVariable(string(p1)*"/"*string(p2))
^(p1::PrintExpression,p2::PrintExpression) = PrintVariable(string(p1)*"^"*string(p2))

for o in [:*,:/,:^]
    for T in Reals
        @eval begin
            $o(p1::$T,p2::PrintExpression) = PrintVariable(""*string(p1)*string($o)*string(p2))
            $o(p1::PrintExpression,p2::$T) = PrintVariable(string(p1)*string($o)*string(p2)*"")
        end
    end
end

for (M,f,nargs) in diffrules
    if f != :^
        if nargs == 1
            @eval begin
                $f(p::PrintExpression) = PrintVariable(string($f)*"("*raw(p)*")")
            end
        elseif nargs == 2
            for T in Reals
                @eval begin
                    $f(a::$T,p::PrintExpression) = PrintVariable(string($f)*"("*string(a)*" , "*raw(p)*")")
                    $f(p::PrintExpression,a::$T) = PrintVariable(string($f)*"("*raw(p)*" , "*string(a)*")")
                end
            end
            @eval begin
                $f(p1::PrintExpression,p2::PrintExpression) = PrintVariable(string($f)*"("*raw(p1)*" , "*raw(p2)*")")
            end
        end
    end
end
