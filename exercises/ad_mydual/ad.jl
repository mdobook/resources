struct Dual{TF}
    value::TF
    deriv::TF
end

Base.:sin(x::Dual) = Dual(sin(x.value), cos(x.value)*x.deriv)
Base.:cos(x::Dual) = Dual( , )
Base.:exp(x::Dual) = Dual( , )
Base.:sqrt(x::Dual) = Dual( , )
Base.:+(x::Dual, y::Dual) = Dual(x.value + y.value, x.deriv + y.deriv)
Base.:/(x::Dual, y::Dual) = Dual( , )
Base.:^(x::Dual, p::Number) = Dual( , )

func(x) = exp(x)/sqrt(sin(x)^3 + cos(x)^3)

x = Dual(2.0, 1.0)
y = func(x)
println(y)
# exact derivative = 18.880984872278777