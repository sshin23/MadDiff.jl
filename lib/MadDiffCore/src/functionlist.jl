@inline mone(x)=-one(x) 
@inline one(x1,x2)=one(x1)
@inline zero(x1,x2)=zero(x1)
@inline mone(x1,x2)=-one(x1)
@inline x1(x1,x2)=x1
@inline x2(x1,x2)=x2
@inline and(x::Bool,y::Bool) = x && y
@inline or(x::Bool,y::Bool) = x || y
@inline and(x,y::Bool) = x==1 && y
@inline or(x,y::Bool) = x==1 || y
@inline and(x::Bool,y) = x && y==1
@inline or(x::Bool,y) = x || y==1
@inline and(x,y) = x==1 && y==1
@inline or(x,y) = x==1 || y==1

const f_nargs_1 = [
    (
        :+,
        +,
        one,
        zero
    ),
    (
        :-,
        -,
        mone,
        zero
    ),
    (
        :inv,
        inv,
        x -> -abs2(inv(x)),
        x -> -(2*inv(x))*(-abs2(inv(x)))
    ),
    (
        :abs,
        abs,
        x->(ifelse(x >= 0, one(x), -one(x))),
        zero
    ),
    (
        :sqrt,
        sqrt,
        x->(0.5 / sqrt(x)),
        x->((0.5 * -(0.5 / sqrt(x))) / sqrt(x) ^ 2)
    ),
    (
        :(Base.sqrt),
        sqrt,
        x->(0.5 / sqrt(x)),
        x->((0.5 * -(0.5 / sqrt(x))) / sqrt(x) ^ 2)
    ),
    (
        :cbrt,
        cbrt,
        x->(0.3333333333333333 / cbrt(x) ^ 2),
        x->((0.3333333333333333 * -(2 * (0.3333333333333333 / cbrt(x) ^ 2) * cbrt(x))) / (cbrt(x) ^ 2) ^ 2)
    ),
    (
        :abs2,
        abs2,
        x -> 2x,
        x -> 2
    ),
    (
        :exp,
        exp,
        exp,
        exp
    ),
    (
        :exp2,
        exp2,
        x -> exp2(x)*0.69314718055994528622676398299518041312694549560546875,
        x -> exp2(x)*0.69314718055994528622676398299518041312694549560546875*0.69314718055994528622676398299518041312694549560546875
    ),
    (
        :exp10,
        exp10,
        x -> exp10(x)*2.30258509299404590109361379290930926799774169921875,
        x -> exp10(x)*2.30258509299404590109361379290930926799774169921875*2.30258509299404590109361379290930926799774169921875
    ),
    (
        :log,
        log,
        inv,
        x -> -abs2(inv(x))
    ),
    (
        :(Base.log),
        log,
        inv,
        x -> -abs2(inv(x))
    ),
    (
        :log2,
        log2,
        x -> inv(x)/0.69314718055994528622676398299518041312694549560546875,
        x -> (-abs2(inv(x)))/0.69314718055994528622676398299518041312694549560546875
    ),
    (
        :(Base.log2),
        log2,
        x -> inv(x)/0.69314718055994528622676398299518041312694549560546875,
        x -> (-abs2(inv(x)))/0.69314718055994528622676398299518041312694549560546875
    ),
    (
        :log1p,
        log1p,
        x->(1 / (1 + x)),
        x->(-1 / (1 + x) ^ 2)
    ),
    (
        :(Base.log1p),
        log1p,
        x->(1 / (1 + x)),
        x->(-1 / (1 + x) ^ 2)
    ),
    (
        :log10,
        log10,
        x -> inv(x)/2.30258509299404590109361379290930926799774169921875,
        x -> (-abs2(inv(x)))/2.30258509299404590109361379290930926799774169921875
    ),
    (
        :(Base.log10),
        log10,
        x -> inv(x)/2.30258509299404590109361379290930926799774169921875,
        x -> (-abs2(inv(x)))/2.30258509299404590109361379290930926799774169921875
    ),
    (
        :sin,
        sin,
        cos,
        x -> -sin(x)
    ),
    (
        :(Base.sin),
        sin,
        cos,
        x -> -sin(x)
    ),
    (
        :cos,
        cos,
        x -> -sin(x),
        x -> -cos(x)
    ),
    (
        :(Base.cos),
        cos,
        x -> -sin(x),
        x -> -cos(x)
    ),
    (
        :tan,
        tan,
        x -> 1 + tan(x)^2,
        x -> 2*tan(x)*(1 + tan(x)^2)
    ),
    (
        :(Base.tan),
        tan,
        x -> 1 + tan(x)^2,
        x -> 2*tan(x)*(1 + tan(x)^2)
    ),
    (
        :asin,
        asin,
        x->(1 / sqrt(1 - x ^ 2)),
        x->(-(-(2x) * (0.5 / sqrt(1 - x ^ 2))) / sqrt(1 - x ^ 2) ^ 2)
    ),
    (
        :(Base.asin),
        asin,
        x->(1 / sqrt(1 - x ^ 2)),
        x->(-(-(2x) * (0.5 / sqrt(1 - x ^ 2))) / sqrt(1 - x ^ 2) ^ 2)
    ),
    (
        :acos,
        acos,
        x->(-1 / sqrt(1 - x ^ 2)),
        x->(-(-(-(2x) * (0.5 / sqrt(1 - x ^ 2)))) / sqrt(1 - x ^ 2) ^ 2)
    ),
    (
        :(Base.acos),
        acos,
        x->(-1 / sqrt(1 - x ^ 2)),
        x->(-(-(-(2x) * (0.5 / sqrt(1 - x ^ 2)))) / sqrt(1 - x ^ 2) ^ 2)
    ),
    (
        :csc,
        csc,
        x -> (-csc(x))*cot(x),
        x -> (-(-csc(x))*cot(x))*cot(x) + (-(1 + cot(x)^2))*(-csc(x))
    ),
    (
        :sec,
        sec,
        x -> sec(x)*tan(x),
        x -> sec(x)*tan(x)*tan(x) + (1 + tan(x)^2)*sec(x)
    ),
    (
        :cot,
        cot,
        x -> -(1 + cot(x)^2),
        x -> -2*cot(x)*(-(1 + cot(x)^2))
    ),
    (
        :atan,
        atan,
        x -> inv(1 + x^2),
        x -> (-abs2(inv(1 + x^2)))*2x
    ),
    (
        :acot,
        acot,
        x -> -inv(1 + x^2),
        x -> -(-abs2(inv(1 + x^2)))*2x
    ),
    (
        :sind,
        sind,
        x -> 0.0174532925199432954743716805978692718781530857086181640625*cosd(x),
        x -> 0.0174532925199432954743716805978692718781530857086181640625*-0.0174532925199432954743716805978692718781530857086181640625*sind(x)
    ),
    (
        :cosd,
        cosd,
        x -> -0.0174532925199432954743716805978692718781530857086181640625*sind(x),
        x -> -0.0174532925199432954743716805978692718781530857086181640625*0.0174532925199432954743716805978692718781530857086181640625*cosd(x)
    ),
    (
        :tand,
        tand,
        x -> 0.0174532925199432954743716805978692718781530857086181640625*(1 + tand(x)^2),
        x -> 0.0174532925199432954743716805978692718781530857086181640625*2*tand(x)*0.0174532925199432954743716805978692718781530857086181640625*(1 + tand(x)^2)
    ),
    (
        :cscd,
        cscd,
        x -> -0.0174532925199432954743716805978692718781530857086181640625*cscd(x)*cotd(x),
        x -> -0.0174532925199432954743716805978692718781530857086181640625*-0.0174532925199432954743716805978692718781530857086181640625*cscd(x)*cotd(x)*cotd(x) + -0.0174532925199432954743716805978692718781530857086181640625*(1 + cotd(x)^2)*-0.0174532925199432954743716805978692718781530857086181640625*cscd(x)
    ),
    (
        :secd,
        secd,
        x -> 0.0174532925199432954743716805978692718781530857086181640625*secd(x)*tand(x),
        x -> 0.0174532925199432954743716805978692718781530857086181640625*0.0174532925199432954743716805978692718781530857086181640625*secd(x)*tand(x)*tand(x) + 0.0174532925199432954743716805978692718781530857086181640625*(1 + tand(x)^2)*0.0174532925199432954743716805978692718781530857086181640625*secd(x)
    ),
    (
        :cotd,
        cotd,
        x -> -0.0174532925199432954743716805978692718781530857086181640625*(1 + cotd(x)^2),
        x -> -0.0174532925199432954743716805978692718781530857086181640625*2*cotd(x)*-0.0174532925199432954743716805978692718781530857086181640625*(1 + cotd(x)^2)
    ),
    (
        :atand,
        atand,
        x -> 57.29577951308232286464772187173366546630859375/(1 + x^2),
        x -> -57.29577951308232286464772187173366546630859375*2*x/(1 + x^2)^2
    ),
    (
        :acotd,
        acotd,
        x -> -57.29577951308232286464772187173366546630859375/(1 + x^2),
        x -> 57.29577951308232286464772187173366546630859375*2*x/(1 + x^2)^2
    ),
    (
        :sinh,
        sinh,
        cosh,
        sinh
    ),
    (
        :cosh,
        cosh,
        sinh,
        cosh
    ),
    (
        :tanh,
        tanh,
        x -> 1 - tanh(x)^2,
        x -> -2*tanh(x)*(1 - tanh(x)^2)
    ),
    (
        :csch,
        csch,
        x -> (-coth(x))*csch(x),
        x -> csch(x)^2*csch(x) + (-coth(x))*csch(x)*(-coth(x))
    ),
    (
        :sech,
        sech,
        x -> (-tanh(x))*sech(x),
        x -> (-(1 - tanh(x)^2))*sech(x) + (-tanh(x))*sech(x)*(-tanh(x))
    ),
    (
        :coth,
        coth,
        x -> -csch(x)^2,
        x -> -2*csch(x)*(-coth(x))*csch(x)
    ),
    (
        :atanh,
        atanh,
        x -> inv(1 - x^2),
        x -> (-abs2(inv(1 - x^2)))*(-2x)
    ),
    (
        :(Base.atanh),
        atanh,
        x -> inv(1 - x^2),
        x -> (-abs2(inv(1 - x^2)))*(-2x)
    ),
    (
        :acoth,
        acoth,
        x -> inv(1 - x^2),
        x -> (-abs2(inv(1 - x^2)))*(-2x)
    ),
    (
        :erfi,
        erfi,
        x -> 1.1283791670955125585606992899556644260883331298828125*exp(x^2),
        x -> 1.1283791670955125585606992899556644260883331298828125*exp(x^2)*2x
    ),
    (
        :erfcinv,
        erfcinv,
        x -> -0.8862269254527579409597137782839126884937286376953125*exp(erfcinv(x)^2),
        x -> -0.8862269254527579409597137782839126884937286376953125*exp(erfcinv(x)^2)*2*erfcinv(x)*-0.8862269254527579409597137782839126884937286376953125*exp(erfcinv(x)^2)
    ),
    (
        :erfcx,
        erfcx,
        x -> 2*x*erfcx(x) - 1.1283791670955125585606992899556644260883331298828125,
        x -> 2*erfcx(x) + (2*x*erfcx(x) - 1.1283791670955125585606992899556644260883331298828125)*2*x
    ),
    (
        :invdigamma,
        invdigamma,
        x -> inv(trigamma(invdigamma(x))),
        x -> (-abs2(inv(trigamma(invdigamma(x)))))*polygamma(2 , invdigamma(x))*inv(trigamma(invdigamma(x)))
    ),
    (
        :bessely1,
        bessely1,
        x -> (bessely0(x) - bessely(2 , x))/2,
        x -> (-bessely1(x) + -(bessely(1 , x) - bessely(3 , x))/2)/2
    ),
    (
        :besselj1,
        besselj1,
        x -> (besselj0(x) - besselj(2 , x))/2,
        x -> (-besselj1(x) + -(besselj(1 , x) - besselj(3 , x))/2)/2
    ),
    (
        :dawson,
        dawson,
        x -> 1 - 2*x*dawson(x),
        x -> -(2*dawson(x) + (1 - 2*x*dawson(x))*2*x)
    ),
    (
        :airyaiprime,
        airyaiprime,
        x -> x*airyai(x),
        x -> airyai(x) + airyaiprime(x)*x
    ),
    (
        :erf,
        erf,
        x -> 1.1283791670955125585606992899556644260883331298828125*exp(-x*x),
        x -> 1.1283791670955125585606992899556644260883331298828125*exp(-x*x)*(-2x)
    ),
    (
        :digamma,
        digamma,
        trigamma,
        x -> polygamma(2 , x)
    ),
    (
        :gamma,
        gamma,
        x -> digamma(x)*gamma(x),
        x -> trigamma(x)*gamma(x) + digamma(x)*gamma(x)*digamma(x)
    ),
    (
        :airyai,
        airyai,
        airyaiprime,
        x -> x*airyai(x)
    ),
    (
        :airybi,
        airybi,
        airybiprime,
        x -> x*airybi(x)
    ),
    (
        :erfinv,
        erfinv,
        x -> 0.8862269254527579409597137782839126884937286376953125*exp(erfinv(x)^2),
        x -> 0.8862269254527579409597137782839126884937286376953125*exp(erfinv(x)^2)*2*erfinv(x)*0.8862269254527579409597137782839126884937286376953125*exp(erfinv(x)^2)
    ),
    (
        :bessely0,
        bessely0,
        x -> -bessely1(x),
        x -> -(bessely0(x) - bessely(2 , x))/2
    ),
    (
        :erfc,
        erfc,
        x -> -1.1283791670955125585606992899556644260883331298828125*exp(-x*x),
        x -> -1.1283791670955125585606992899556644260883331298828125*exp(-x*x)*(-2x)
    ),
    (
        :trigamma,
        trigamma,
        x -> polygamma(2 , x),
        x -> polygamma(3 , x)
    ),
    (
        :airybiprime,
        airybiprime,
        x -> x*airybi(x),
        x -> airybi(x) + airybiprime(x)*x
    ),
    (
        :besselj0,
        besselj0,
        x -> -besselj1(x),
        x -> -(besselj0(x) - besselj(2 , x))/2
    )
]



const f_nargs_2 = [
    (
        :+,
        +,
        one,
        one,
        zero,
        zero,
        zero
    ),
    (
        :-,
        -,
        one,
        mone,
        zero,
        zero,
        zero
    ),
    (
        :*,
        *,
        x2,
        x1,
        zero,
        one,
        zero
    ),
    (
        :^,
        ^,
        ((x1,x2)  -> x2*x1^(x2 - 1)),
        ((x1,x2)  -> x1^x2*log(x1)),
        ((x1,x2)  -> x2*(x2 - 1)*x1^(x2 - 2)),
        ((x1,x2)  -> x2*x1^(x2  -  1)*log(x1) + x1^(x2 - 1)),
        ((x1,x2)  -> x1^x2*log(x1)*log(x1))),
    (
        :/,
        /,
        ((x1,x2)  -> 1/x2),
        ((x1,x2)  -> -x1/x2^2),
        zero,
        ((x1,x2)  -> -1/x2^2),
        ((x1,x2)  -> 2x1/x2^3)
    ),
    (
        :beta,
        beta,
        (x1,x2) -> beta(x1 , x2)*(digamma(x1) - digamma(x1 + x2)),
        (x1,x2) -> beta(x1 , x2)*(digamma(x2) - digamma(x1 + x2)),
        (x1,x2) -> beta(x1 , x2)*(digamma(x1) - digamma(x1 + x2))*(digamma(x1) - digamma(x1 + x2)) + (trigamma(x1) + -trigamma(x1 + x2))*beta(x1 , x2),
        (x1,x2) -> beta(x1 , x2)*(digamma(x2) - digamma(x1 + x2))*(digamma(x1) - digamma(x1 + x2)) + (-trigamma(x1 + x2))*beta(x1 , x2),
        (x1,x2) -> beta(x1 , x2)*(digamma(x2) - digamma(x1 + x2))*(digamma(x2) - digamma(x1 + x2)) + (trigamma(x2) + -trigamma(x1 + x2))*beta(x1 , x2)
    ),
    (
        :logbeta,
        logbeta,
        (x1,x2) -> digamma(x1) - digamma(x1 + x2),
        (x1,x2) -> digamma(x2) - digamma(x1 + x2),
        (x1,x2) -> trigamma(x1) + -trigamma(x1 + x2),
        (x1,x2) -> -trigamma(x1 + x2),
        (x1,x2) -> trigamma(x2) + -trigamma(x1 + x2)
    ),
    (
        :(<=),
        (<=),
        zero,
        zero,
        zero,
        zero,
        zero
    ),
    (
        :(>=),
        (>=),
        zero,
        zero,
        zero,
        zero,
        zero
    ),
    (
        :(==),
        (==),
        zero,
        zero,
        zero,
        zero,
        zero
    ),
    (
        :<,
        <,
        zero,
        zero,
        zero,
        zero,
        zero
    ),
    (
        :>,
        >,
        zero,
        zero,
        zero,
        zero,
        zero
    ),
    (
        :and,
        and,
        zero,
        zero,
        zero,
        zero,
        zero
    ),
    (
        :or,
        or,
        zero,
        zero,
        zero,
        zero,
        zero
    ),
]
