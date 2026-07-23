#!/usr/bin/env python3
# 2026-07-22: xi-GENERIC coincidence-shape certificates for {734,735}.
# These re-prove the coincidence leaves WITHOUT the gauge shift (xi free),
# so each sub-case kill stands alone, gauge-free.
#  a2-coinc-xi: Ac=Bc=0, h=h2(r-xi)^2, K=k3(r-xi)^3, P=p1(r-xi)
#               => c = c0 + (p1+2*h2*xi)(r-xi)^2/2   (from P = 2h + c' - r h')
#  b-coinc-xi:  Bc=0, h=h2(r-xi)^2, Ac=a2(r-xi)^2, K=k3(r-xi)^3, c generic
import sympy as sp, subprocess

r, u0 = sp.symbols('r u0')
al = sp.Symbol('al')
N = sp.sympify(open('q8_layers.srepr').read())
S = sp.Symbol

def cert(tag, eqs, vars_, hyps, targets, cap=1200):
    coeffs = []
    for e in eqs:
        e = sp.expand(e)
        if e == 0:
            continue
        for x in sp.Poly(e, r).all_coeffs():
            if x == 0:
                continue
            num, _ = sp.fraction(sp.together(x))
            coeffs.append(sp.expand(num))
    svars = [f's{i}' for i in range(len(hyps))]
    L = [f'ring rr = 0, ({", ".join(vars_ + svars + ["tt"])}), dp;']
    gens = [str(x).replace('**', '^') for x in coeffs]
    gens += [f'{svars[i]}*({str(hp).replace("**", "^")})-1' for i, hp in enumerate(hyps)]
    L.append('ideal LD = ' + ','.join(gens) + ';')
    for name, expr in targets.items():
        es = str(sp.expand(expr)).replace('**', '^')
        L.append(f'ideal J{name} = LD, tt*({es})-1;')
        L.append(f'if (reduce(1, groebner(J{name})) == 0) {{ "UNIT {name}"; }} else {{ "NONUNIT {name}"; }}')
    L.append('exit;')
    try:
        p = subprocess.run(['/usr/bin/Singular', '-q'], input='\n'.join(L) + '\n',
                           capture_output=True, text=True, timeout=cap)
        out = ' | '.join(l for l in p.stdout.strip().splitlines() if l.strip())
    except subprocess.TimeoutExpired:
        out = f'TIMEOUT {cap}s (UNDECIDED)'
    print(f'CERT {tag}: {out}', flush=True)

xi, h2, k3, p1, a2, c0 = sp.symbols('xi h2 k3 p1 a2 c0')
z = r - xi

# --- a2-coinc-xi
sub = {S('ac0'): 0, S('ac1'): 0, S('ac2'): 0,
       S('bc0'): 0, S('bc1'): 0, S('bc2'): 0, S('bc3'): 0}
hpoly = sp.expand(h2*z**2)
Kpoly = sp.expand(k3*z**3)
cpoly = sp.expand(c0 + (p1 + 2*h2*xi)*z**2/2)
for j in range(3):
    sub[S(f'bh{j}')] = hpoly.coeff(r, j)
    sub[S(f'ch{j}')] = cpoly.coeff(r, j)
for j in range(4):
    sub[S(f'cc{j}')] = Kpoly.coeff(r, j)
# sanity: P = 2h + c' - r h' == p1*(r-xi)
P = sp.expand(2*hpoly + sp.diff(cpoly, r) - r*sp.diff(hpoly, r))
print('CHECK a2-coinc-xi P == p1*(r-xi):',
      'PASS' if sp.expand(P - p1*z) == 0 else 'FAIL')
eqs = [N[m].subs(sub) for m in range(2, 9)]
cert('a2-coinc-xi', eqs, ['h2', 'k3', 'p1', 'c0', 'xi', 'al'],
     [al, h2, k3], {'p1': p1})

# --- b-coinc-xi
sub = {S('bc0'): 0, S('bc1'): 0, S('bc2'): 0, S('bc3'): 0}
Apoly = sp.expand(a2*z**2)
for j in range(3):
    sub[S(f'bh{j}')] = hpoly.coeff(r, j)
    sub[S(f'ac{j}')] = Apoly.coeff(r, j)
for j in range(4):
    sub[S(f'cc{j}')] = Kpoly.coeff(r, j)
eqs = [N[m].subs(sub) for m in range(2, 9)]
cert('b-coinc-xi', eqs, ['h2', 'k3', 'a2', 'ch0', 'ch1', 'ch2', 'xi', 'al'],
     [al, h2, k3], {'a2': a2})
