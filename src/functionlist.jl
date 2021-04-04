struct FOne1 <: Function end; @inline (::FOne1)(x)=1.; f_one1 = FOne1()
struct FZero1 <: Function end; @inline (::FZero1)(x)=0.;  f_zero1 = FZero1()
struct FMOne1 <: Function end; @inline (::FMOne1)(x)=-1.;  f_mone1 = FMOne1()
struct FOne2 <: Function end; @inline (::FOne2)(x1,x2)=1.; f_one2 = FOne2()
struct FZero2 <: Function end; @inline (::FZero2)(x1,x2)=0.;  f_zero2 = FZero2()
struct FMOne2 <: Function end; @inline (::FMOne2)(x1,x2)=-1.;  f_mone2 = FMOne2()
struct FX1 <: Function end; @inline (::FX1)(x1,x2)=x1;  f_x1 = FX1()
struct FX2 <: Function end; @inline (::FX2)(x1,x2)=x2;  f_x2 = FX2()

const f_nargs_1 = [
    (
        :+,
        f_one1,
        f_zero1
    ),
    (
        :-,
        f_mone1,
        f_zero1
    ),
    (
        :inv,
        x -> -abs2(inv(x)),
        x -> -(2*inv(x))*(-abs2(inv(x)))
    ),
    (
        :abs2,
        x -> 2x,
        x -> 2
    ),
    (
        :exp,
        exp,
        exp
    ),
    (
        :exp2,
        x -> exp2(x)*0.69314718055994528622676398299518041312694549560546875,
        x -> exp2(x)*0.69314718055994528622676398299518041312694549560546875*0.69314718055994528622676398299518041312694549560546875
    ),
    (
        :exp10,
        x -> exp10(x)*2.30258509299404590109361379290930926799774169921875,
        x -> exp10(x)*2.30258509299404590109361379290930926799774169921875*2.30258509299404590109361379290930926799774169921875
    ),
    (
        :log,
        inv,
        x -> -abs2(inv(x))
    ),
    (
        :log2,
        x -> inv(x)/0.69314718055994528622676398299518041312694549560546875,
        x -> (-abs2(inv(x)))/0.69314718055994528622676398299518041312694549560546875
    ),
    (
        :log10,
        x -> inv(x)/2.30258509299404590109361379290930926799774169921875,
        x -> (-abs2(inv(x)))/2.30258509299404590109361379290930926799774169921875
    ),
    (
        :sin,
        cos,
        x -> -sin(x)
    ),
    (
        :cos,
        x -> -sin(x),
        x -> -cos(x)
    ),
    (
        :tan,
        x -> 1 + tan(x)^2,
        x -> 2*tan(x)*(1 + tan(x)^2)
    ),
    (
        :csc,
        x -> (-csc(x))*cot(x),
        x -> (-(-csc(x))*cot(x))*cot(x) + (-(1 + cot(x)^2))*(-csc(x))
    ),
    (
        :sec,
        x -> sec(x)*tan(x),
        x -> sec(x)*tan(x)*tan(x) + (1 + tan(x)^2)*sec(x)
    ),
    (
        :cot,
        x -> -(1 + cot(x)^2),
        x -> -2*cot(x)*(-(1 + cot(x)^2))
    ),
    (
        :atan,
        x -> inv(1 + x^2),
        x -> (-abs2(inv(1 + x^2)))*2x
    ),
    (
        :acot,
        x -> -inv(1 + x^2),
        x -> -(-abs2(inv(1 + x^2)))*2x
    ),
    (
        :sind,
        x -> 0.0174532925199432954743716805978692718781530857086181640625*cosd(x),
        x -> 0.0174532925199432954743716805978692718781530857086181640625*-0.0174532925199432954743716805978692718781530857086181640625*sind(x)
    ),
    (
        :cosd,
        x -> -0.0174532925199432954743716805978692718781530857086181640625*sind(x),
        x -> -0.0174532925199432954743716805978692718781530857086181640625*0.0174532925199432954743716805978692718781530857086181640625*cosd(x)
    ),
    (
        :tand,
        x -> 0.0174532925199432954743716805978692718781530857086181640625*(1 + tand(x)^2),
        x -> 0.0174532925199432954743716805978692718781530857086181640625*2*tand(x)*0.0174532925199432954743716805978692718781530857086181640625*(1 + tand(x)^2)
    ),
    (
        :cscd,
        x -> -0.0174532925199432954743716805978692718781530857086181640625*cscd(x)*cotd(x),
        x -> -0.0174532925199432954743716805978692718781530857086181640625*-0.0174532925199432954743716805978692718781530857086181640625*cscd(x)*cotd(x)*cotd(x) + -0.0174532925199432954743716805978692718781530857086181640625*(1 + cotd(x)^2)*-0.0174532925199432954743716805978692718781530857086181640625*cscd(x)
    ),
    (
        :secd,
        x -> 0.0174532925199432954743716805978692718781530857086181640625*secd(x)*tand(x),
        x -> 0.0174532925199432954743716805978692718781530857086181640625*0.0174532925199432954743716805978692718781530857086181640625*secd(x)*tand(x)*tand(x) + 0.0174532925199432954743716805978692718781530857086181640625*(1 + tand(x)^2)*0.0174532925199432954743716805978692718781530857086181640625*secd(x)
    ),
    (
        :cotd,
        x -> -0.0174532925199432954743716805978692718781530857086181640625*(1 + cotd(x)^2),
        x -> -0.0174532925199432954743716805978692718781530857086181640625*2*cotd(x)*-0.0174532925199432954743716805978692718781530857086181640625*(1 + cotd(x)^2)
    ),
    (
        :atand,
        x -> 57.29577951308232286464772187173366546630859375/(1 + x^2),
        x -> -57.29577951308232286464772187173366546630859375*2*x/(1 + x^2)^2
    ),
    (
        :acotd,
        x -> -57.29577951308232286464772187173366546630859375/(1 + x^2),
        x -> 57.29577951308232286464772187173366546630859375*2*x/(1 + x^2)^2
    ),
    (
        :sinh,
        cosh,
        sinh
    ),
    (
        :cosh,
        sinh,
        cosh
    ),
    (
        :tanh,
        x -> 1 - tanh(x)^2,
        x -> -2*tanh(x)*(1 - tanh(x)^2)
    ),
    (
        :csch,
        x -> (-coth(x))*csch(x),
        x -> csch(x)^2*csch(x) + (-coth(x))*csch(x)*(-coth(x))
    ),
    (
        :sech,
        x -> (-tanh(x))*sech(x),
        x -> (-(1 - tanh(x)^2))*sech(x) + (-tanh(x))*sech(x)*(-tanh(x))
    ),
    (
        :coth,
        x -> -csch(x)^2,
        x -> -2*csch(x)*(-coth(x))*csch(x)
    ),
    (
        :atanh,
        x -> inv(1 - x^2),
        x -> (-abs2(inv(1 - x^2)))*(-2x)
    ),
    (
        :acoth,
        x -> inv(1 - x^2),
        x -> (-abs2(inv(1 - x^2)))*(-2x)
    ),
    (
        :erfi,
        x -> 1.1283791670955125585606992899556644260883331298828125*exp(x^2),
        x -> 1.1283791670955125585606992899556644260883331298828125*exp(x^2)*2x
    ),
    (
        :loggamma,
        digamma,
        trigamma
    ),
    (
        :erfcinv,
        x -> -0.8862269254527579409597137782839126884937286376953125*exp(erfcinv(x)^2),
        x -> -0.8862269254527579409597137782839126884937286376953125*exp(erfcinv(x)^2)*2*erfcinv(x)*-0.8862269254527579409597137782839126884937286376953125*exp(erfcinv(x)^2)
    ),
    (
        :erfcx,
        x -> 2*x*erfcx(x) - 1.1283791670955125585606992899556644260883331298828125,
        x -> 2*erfcx(x) + (2*x*erfcx(x) - 1.1283791670955125585606992899556644260883331298828125)*2*x
    ),
    (
        :invdigamma,
        x -> inv(trigamma(invdigamma(x))),
        x -> (-abs2(inv(trigamma(invdigamma(x)))))*polygamma(2 , invdigamma(x))*inv(trigamma(invdigamma(x)))
    ),
    (
        :bessely1,
        x -> (bessely0(x) - bessely(2 , x))/2,
        x -> (-bessely1(x) + -(bessely(1 , x) - bessely(3 , x))/2)/2
    ),
    (
        :besselj1,
        x -> (besselj0(x) - besselj(2 , x))/2,
        x -> (-besselj1(x) + -(besselj(1 , x) - besselj(3 , x))/2)/2
    ),
    (
        :dawson,
        x -> 1 - 2*x*dawson(x),
        x -> -(2*dawson(x) + (1 - 2*x*dawson(x))*2*x)
    ),
    (
        :airyaiprime,
        x -> x*airyai(x),
        x -> airyai(x) + airyaiprime(x)*x
    ),
    (
        :erf,
        x -> 1.1283791670955125585606992899556644260883331298828125*exp(-x*x),
        x -> 1.1283791670955125585606992899556644260883331298828125*exp(-x*x)*(-2x)
    ),
    (
        :digamma,
        trigamma,
        x -> polygamma(2 , x)
    ),
    (
        :gamma,
        x -> digamma(x)*gamma(x),
        x -> trigamma(x)*gamma(x) + digamma(x)*gamma(x)*digamma(x)
    ),
    (
        :airyai,
        airyaiprime,
        x -> x*airyai(x)
    ),
    (
        :airybi,
        airybiprime,
        x -> x*airybi(x)
    ),
    (
        :erfinv,
        x -> 0.8862269254527579409597137782839126884937286376953125*exp(erfinv(x)^2),
        x -> 0.8862269254527579409597137782839126884937286376953125*exp(erfinv(x)^2)*2*erfinv(x)*0.8862269254527579409597137782839126884937286376953125*exp(erfinv(x)^2)
    ),
    (
        :bessely0,
        x -> -bessely1(x),
        x -> -(bessely0(x) - bessely(2 , x))/2
    ),
    (
        :erfc,
        x -> -1.1283791670955125585606992899556644260883331298828125*exp(-x*x),
        x -> -1.1283791670955125585606992899556644260883331298828125*exp(-x*x)*(-2x)
    ),
    (
        :trigamma,
        x -> polygamma(2 , x),
        x -> polygamma(3 , x)
    ),
    (
        :airybiprime,
        x -> x*airybi(x),
        x -> airybi(x) + airybiprime(x)*x
    ),
    (
        :besselj0,
        x -> -besselj1(x),
        x -> -(besselj0(x) - besselj(2 , x))/2
    )
]

const f_nargs_2 = [
    (
        :+,
        f_one2,
        f_one2,
        f_zero2,
        f_zero2,
        f_zero2
    ),
    (
        :-,
        f_one2,
        f_mone2,
        f_zero2,
        f_zero2,
        f_zero2
    ),
    (
        :*,
        f_x2,
        f_x1,
        f_zero2,
        f_one2,
        f_zero2
    ),
    (
        :^,
        ((x1,x2)  -> x2*x1^(x2 - 1)),
        ((x1,x2)  -> x1^x2*log(x1)),
        ((x1,x2)  -> x2*(x2 - 1)*x1^(x2 - 2)),
        ((x1,x2)  -> x2*x1^(x2  -  1)*log(x1) + x1^(x2 - 1)),
        ((x1,x2)  -> x1^x2*log(x1)*log(x1))),
    (
        :/,
        ((x1,x2)  -> 1/x2),
        ((x1,x2)  -> -x1/x2^2),
        f_zero2,
        ((x1,x2)  -> -1/x2^2),
        ((x1,x2)  -> 2x1/x2^3)
    ),
    (
        :beta,
        (x1,x2) -> beta(x1 , x2)*(digamma(x1) - digamma(x1 + x2)),
        (x1,x2) -> beta(x1 , x2)*(digamma(x2) - digamma(x1 + x2)),
        (x1,x2) -> beta(x1 , x2)*(digamma(x1) - digamma(x1 + x2))*(digamma(x1) - digamma(x1 + x2)) + (trigamma(x1) + -trigamma(x1 + x2))*beta(x1 , x2),
        (x1,x2) -> beta(x1 , x2)*(digamma(x2) - digamma(x1 + x2))*(digamma(x1) - digamma(x1 + x2)) + (-trigamma(x1 + x2))*beta(x1 , x2),
        (x1,x2) -> beta(x1 , x2)*(digamma(x2) - digamma(x1 + x2))*(digamma(x2) - digamma(x1 + x2)) + (trigamma(x2) + -trigamma(x1 + x2))*beta(x1 , x2)
    ),
    (
        :logbeta,
        (x1,x2) -> digamma(x1) - digamma(x1 + x2),
        (x1,x2) -> digamma(x2) - digamma(x1 + x2),
        (x1,x2) -> trigamma(x1) + -trigamma(x1 + x2),
        (x1,x2) -> -trigamma(x1 + x2),
        (x1,x2) -> trigamma(x2) + -trigamma(x1 + x2)
    )
]
