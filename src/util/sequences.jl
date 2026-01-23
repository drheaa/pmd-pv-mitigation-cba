"""
seq_components(Aa, Ab, Ac)

Symmetrical components:
A0 is zero sequence, A1 is positive, A2 is negative sequence.
"""
function seq_components(Aa::Complex, Ab::Complex, Ac::Complex)
    a = cis(2pi / 3)
    A0 = (Aa + Ab + Ac) / 3
    A1 = (Aa + a * Ab + a^2 * Ac) / 3
    A2 = (Aa + a^2 * Ab + a * Ac) / 3
    return (A0=A0, A1=A1, A2=A2)
end

"""
seq_magnitudes(Aa, Ab, Ac)

Returns magnitudes and negative-to-positive ratio.
A2_over_A1 is commonly used for unbalance metrics.
"""
function seq_magnitudes(Aa::Complex, Ab::Complex, Ac::Complex)
    s = seq_components(Aa, Ab, Ac)
    A0 = abs(s.A0)
    A1 = abs(s.A1)
    A2 = abs(s.A2)
    ratio = A1 > 1e-12 ? (A2 / A1) : NaN
    return (A0=A0, A1=A1, A2=A2, A2_over_A1=ratio)
end
