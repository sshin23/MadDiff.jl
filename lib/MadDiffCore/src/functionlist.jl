@register_univariate(
    +,
    one,
    zero
)
@register_univariate(
    -,
    _mone,
    zero
)
@register_univariate(
    inv,
    x -> -abs2(inv(x)),
    x -> -(2*inv(x))*(-abs2(inv(x)))
)
@register_univariate(
    abs,
    x->(ifelse(x >= 0,
                    one(x), -one(x))),        zero
)
@register_univariate(
    sqrt,
    x->(0.5 / sqrt(x)),
    x->((0.5 * -(0.5 / sqrt(x))) / sqrt(x) ^ 2)
)
@register_univariate(
    cbrt,
    x->(0.3333333333333333 / cbrt(x) ^ 2),
    x->((0.3333333333333333 * -(2 * (0.3333333333333333 / cbrt(x) ^ 2) * cbrt(x))) / (cbrt(x) ^ 2) ^ 2)
)
@register_univariate(
    abs2,
    x -> 2x,
    x -> 2
)
@register_univariate(
    exp,
    exp,
    exp
)
@register_univariate(
    exp2,
    x -> exp2(x)*0.69314718055994528622676398299518041312694549560546875,
    x -> exp2(x)*0.69314718055994528622676398299518041312694549560546875*0.69314718055994528622676398299518041312694549560546875
)
@register_univariate(
    exp10,
    x -> exp10(x)*2.30258509299404590109361379290930926799774169921875,
    x -> exp10(x)*2.30258509299404590109361379290930926799774169921875*2.30258509299404590109361379290930926799774169921875
)
@register_univariate(
    log,
    inv,
    x -> -abs2(inv(x))
)
@register_univariate(
    log2,
    x -> inv(x)/0.69314718055994528622676398299518041312694549560546875,
    x -> (-abs2(inv(x)))/0.69314718055994528622676398299518041312694549560546875
)
@register_univariate(
    log1p,
    x->(1 / (1 + x)),
    x->(-1 / (1 + x) ^ 2)
)
@register_univariate(
    log10,
    x -> inv(x)/2.30258509299404590109361379290930926799774169921875,
    x -> (-abs2(inv(x)))/2.30258509299404590109361379290930926799774169921875
)
@register_univariate(
    sin,
    cos,
    x -> -sin(x)
)
@register_univariate(
    cos,
    x -> -sin(x),
    x -> -cos(x)
)
@register_univariate(
    tan,
    x -> 1 + tan(x)^2,
    x -> 2*tan(x)*(1 + tan(x)^2)
)
@register_univariate(
    asin,
    x->(1 / sqrt(1 - x ^ 2)),
    x->(-(-(2x) * (0.5 / sqrt(1 - x ^ 2))) / sqrt(1 - x ^ 2) ^ 2)
)
@register_univariate(
    acos,
    x->(-1 / sqrt(1 - x ^ 2)),
    x->(-(-(-(2x) * (0.5 / sqrt(1 - x ^ 2)))) / sqrt(1 - x ^ 2) ^ 2)
)
@register_univariate(
    csc,
    x -> (-csc(x))*cot(x),
    x -> (-(-csc(x))*cot(x))*cot(x) + (-(1 + cot(x)^2))*(-csc(x))
)
@register_univariate(
    sec,
    x -> sec(x)*tan(x),
    x -> sec(x)*tan(x)*tan(x) + (1 + tan(x)^2)*sec(x)
)
@register_univariate(
    cot,
    x -> -(1 + cot(x)^2),
    x -> -2*cot(x)*(-(1 + cot(x)^2))
)
@register_univariate(
    atan,
    x -> inv(1 + x^2),
    x -> (-abs2(inv(1 + x^2)))*2x
)
@register_univariate(
    acot,
    x -> -inv(1 + x^2),
    x -> -(-abs2(inv(1 + x^2)))*2x
)
@register_univariate(
    sind,
    x -> 0.0174532925199432954743716805978692718781530857086181640625*cosd(x),
    x -> 0.0174532925199432954743716805978692718781530857086181640625*-0.0174532925199432954743716805978692718781530857086181640625*sind(x)
)
@register_univariate(
    cosd,
    x -> -0.0174532925199432954743716805978692718781530857086181640625*sind(x),
    x -> -0.0174532925199432954743716805978692718781530857086181640625*0.0174532925199432954743716805978692718781530857086181640625*cosd(x)
)
@register_univariate(
    tand,
    x -> 0.0174532925199432954743716805978692718781530857086181640625*(1 + tand(x)^2),
    x -> 0.0174532925199432954743716805978692718781530857086181640625*2*tand(x)*0.0174532925199432954743716805978692718781530857086181640625*(1 + tand(x)^2)
)
@register_univariate(
    cscd,
    x -> -0.0174532925199432954743716805978692718781530857086181640625*cscd(x)*cotd(x),
    x -> -0.0174532925199432954743716805978692718781530857086181640625*-0.0174532925199432954743716805978692718781530857086181640625*cscd(x)*cotd(x)*cotd(x) + -0.0174532925199432954743716805978692718781530857086181640625*(1 + cotd(x)^2)*-0.0174532925199432954743716805978692718781530857086181640625*cscd(x)
)
@register_univariate(
    secd,
    x -> 0.0174532925199432954743716805978692718781530857086181640625*secd(x)*tand(x),
    x -> 0.0174532925199432954743716805978692718781530857086181640625*0.0174532925199432954743716805978692718781530857086181640625*secd(x)*tand(x)*tand(x) + 0.0174532925199432954743716805978692718781530857086181640625*(1 + tand(x)^2)*0.0174532925199432954743716805978692718781530857086181640625*secd(x)
)
@register_univariate(
    cotd,
    x -> -0.0174532925199432954743716805978692718781530857086181640625*(1 + cotd(x)^2),
    x -> -0.0174532925199432954743716805978692718781530857086181640625*2*cotd(x)*-0.0174532925199432954743716805978692718781530857086181640625*(1 + cotd(x)^2)
)
@register_univariate(
    atand,
    x -> 57.29577951308232286464772187173366546630859375/(1 + x^2),
    x -> -57.29577951308232286464772187173366546630859375*2*x/(1 + x^2)^2
)
@register_univariate(
    acotd,
    x -> -57.29577951308232286464772187173366546630859375/(1 + x^2),
    x -> 57.29577951308232286464772187173366546630859375*2*x/(1 + x^2)^2
)
@register_univariate(
    sinh,
    cosh,
    sinh
)
@register_univariate(
    cosh,
    sinh,
    cosh
)
@register_univariate(
    tanh,
    x -> 1 - tanh(x)^2,
    x -> -2*tanh(x)*(1 - tanh(x)^2)
)
@register_univariate(
    csch,
    x -> (-coth(x))*csch(x),
    x -> csch(x)^2*csch(x) + (-coth(x))*csch(x)*(-coth(x))
)
@register_univariate(
    sech,
    x -> (-tanh(x))*sech(x),
    x -> (-(1 - tanh(x)^2))*sech(x) + (-tanh(x))*sech(x)*(-tanh(x))
)
@register_univariate(
    coth,
    x -> -csch(x)^2,
    x -> -2*csch(x)*(-coth(x))*csch(x)
)
@register_univariate(
    atanh,
    x -> abs(x) > 1.0 ? NaN : inv(1 - x^2),
    x -> abs(x) > 1.0 ? NaN : (-abs2(inv(1 - x^2)))*(-2x)
)
@register_univariate(
    acoth,
    x -> abs(x) < 1.0 ? NaN : inv(1 - x^2),
    x -> abs(x) < 1.0 ? NaN : (-abs2(inv(1 - x^2)))*(-2x)
)

@register_bivariate(
    +,
    _one,
    _one,
    _zero,
    _zero,
    _zero
)
@register_bivariate(
    -,
    _one,
    _mone,
    _zero,
    _zero,
    _zero
)
@register_bivariate(
    *,
    _x2,
    _x1,
    _zero,
    _one,
    _zero
)
@register_bivariate(
    ^,
    ((x1,x2)  -> x2*x1^(x2 - 1)),
    ((x1,x2)  -> x1^x2*log(x1)),
    ((x1,x2)  -> x2*(x2 - 1)*x1^(x2 - 2)),
    ((x1,x2)  -> x2*x1^(x2  -  1)*log(x1) + x1^(x2 - 1)),
    ((x1,x2)  -> x1^x2*log(x1)*log(x1)))
@register_bivariate(
    /,
    ((x1,x2)  -> 1/x2),
    ((x1,x2)  -> -x1/x2^2),
    _zero,
    ((x1,x2)  -> -1/x2^2),
    ((x1,x2)  -> 2x1/x2^3),
)
@register_bivariate(
    (<=),
    _zero,
    _zero,
    _zero,
    _zero,
    _zero
)
@register_bivariate(
    (>=),
    _zero,
    _zero,
    _zero,
    _zero,
    _zero
)
@register_bivariate(
    (==),
    _zero,
    _zero,
    _zero,
    _zero,
    _zero
)
@register_bivariate(
    <,
    _zero,
    _zero,
    _zero,
    _zero,
    _zero
)
@register_bivariate(
    >,
    _zero,
    _zero,
    _zero,
    _zero,
    _zero
)
@register_bivariate(
    _and,
    _zero,
    _zero,
    _zero,
    _zero,
    _zero
)
@register_bivariate(
    _or,
    _zero,
    _zero,
    _zero,
    _zero,
    _zero
)


