module MadDiffSpecialFunctions

using MadDiffCore

import SpecialFunctions: SpecialFunctions, erfi, bessely, besselj, loggamma, erfcinv, hankelh2, hankelh1, erfcx, besselk, beta, invdigamma, bessely1, besselj1, dawson, airyaiprime, erf, digamma, gamma, airyai, airybi, erfinv, bessely0, erfc, trigamma, besseli, polygamma, logbeta, airybiprime, besselj0 

const _UNIVARIATE_FUNCTIONS = Function[]
const _BIVARIATE_FUNCTIONS = Function[]

@register_univariate(
    erfi,
    x -> 1.1283791670955125585606992899556644260883331298828125*exp(x^2),
    x -> 1.1283791670955125585606992899556644260883331298828125*exp(x^2)*2x
)
@register_univariate(
    erfcinv,
    x -> -0.8862269254527579409597137782839126884937286376953125*exp(erfcinv(x)^2),
    x -> -0.8862269254527579409597137782839126884937286376953125*exp(erfcinv(x)^2)*2*erfcinv(x)*-0.8862269254527579409597137782839126884937286376953125*exp(erfcinv(x)^2)
)
@register_univariate(
    erfcx,
    x -> 2*x*erfcx(x) - 1.1283791670955125585606992899556644260883331298828125,
    x -> 2*erfcx(x) + (2*x*erfcx(x) - 1.1283791670955125585606992899556644260883331298828125)*2*x
)
@register_univariate(
    invdigamma,
    x -> inv(trigamma(invdigamma(x))),
    x -> (-abs2(inv(trigamma(invdigamma(x)))))*polygamma(2 , invdigamma(x))*inv(trigamma(invdigamma(x)))
)
@register_univariate(
    bessely1,
    x -> (bessely0(x) - bessely(2 ,
                                x))/2,        x -> (-bessely1(x) + -(bessely(1 , x) - bessely(3 , x))/2)/2
)
@register_univariate(
    besselj1,
    x -> (besselj0(x) - besselj(2 ,
                                x))/2,        x -> (-besselj1(x) + -(besselj(1 , x) - besselj(3 , x))/2)/2
)
@register_univariate(
    dawson,
    x -> 1 - 2*x*dawson(x),
    x -> -(2*dawson(x) + (1 - 2*x*dawson(x))*2*x)
)
@register_univariate(
    airyaiprime,
    x -> x*airyai(x),
    x -> airyai(x) + airyaiprime(x)*x
)
@register_univariate(
    erf,
    x -> 1.1283791670955125585606992899556644260883331298828125*exp(-x*x),
    x -> 1.1283791670955125585606992899556644260883331298828125*exp(-x*x)*(-2x)
)
@register_univariate(
    digamma,
    trigamma,
    x -> polygamma(2 , x)
)
@register_univariate(
    gamma,
    x -> digamma(x)*gamma(x),
    x -> trigamma(x)*gamma(x) + digamma(x)*gamma(x)*digamma(x)
)
@register_univariate(
    airyai,
    airyaiprime,
    x -> x*airyai(x)
)
@register_univariate(
    airybi,
    airybiprime,
    x -> x*airybi(x)
)
@register_univariate(
    erfinv,
    x -> 0.8862269254527579409597137782839126884937286376953125*exp(erfinv(x)^2),
    x -> 0.8862269254527579409597137782839126884937286376953125*exp(erfinv(x)^2)*2*erfinv(x)*0.8862269254527579409597137782839126884937286376953125*exp(erfinv(x)^2)
)
@register_univariate(
    bessely0,
    x -> -bessely1(x),
    x -> -(bessely0(x) - bessely(2 , x))/2
)
@register_univariate(
    erfc,
    x -> -1.1283791670955125585606992899556644260883331298828125*exp(-x*x),
    x -> -1.1283791670955125585606992899556644260883331298828125*exp(-x*x)*(-2x)
)
@register_univariate(
    trigamma,
    x -> polygamma(2 ,
                   x),        x -> polygamma(3 , x)
)
@register_univariate(
    airybiprime,
    x -> x*airybi(x),
    x -> airybi(x) + airybiprime(x)*x
)
@register_univariate(
    besselj0,
    x -> -besselj1(x),
    x -> -(besselj0(x) - besselj(2 , x))/2
)

@register_bivariate(
    beta,
    (x1,x2) -> beta(x1 , x2)*(digamma(x1) - digamma(x1 + x2)),
    (x1,x2) -> beta(x1 , x2)*(digamma(x2) - digamma(x1 + x2)),
    (x1,x2) -> beta(x1 , x2)*(digamma(x1) - digamma(x1 + x2))*(digamma(x1) - digamma(x1 + x2)) + (trigamma(x1) + -trigamma(x1 + x2))*beta(x1 , x2),
    (x1,x2) -> beta(x1 , x2)*(digamma(x2) - digamma(x1 + x2))*(digamma(x1) - digamma(x1 + x2)) + (-trigamma(x1 + x2))*beta(x1 , x2),
    (x1,x2) -> beta(x1 , x2)*(digamma(x2) - digamma(x1 + x2))*(digamma(x2) - digamma(x1 + x2)) + (trigamma(x2) + -trigamma(x1 + x2))*beta(x1 , x2)
)
@register_bivariate(
    logbeta,
    (x1,x2) -> digamma(x1) - digamma(x1 + x2),
    (x1,x2) -> digamma(x2) - digamma(x1 + x2),
    (x1,x2) -> trigamma(x1) + -trigamma(x1 + x2),
    (x1,x2) -> -trigamma(x1 + x2),
    (x1,x2) -> trigamma(x2) + -trigamma(x1 + x2)
)

end # module
