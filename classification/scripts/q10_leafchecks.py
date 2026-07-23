#!/usr/bin/env python3
# 2026-07-22: {734,735} hand-tree leaf verification battery.
# Part 1: verify the structural master forms for branches a2, b against
#         the srepr layers (exact polynomial equality).
# Part 2: per-leaf Singular Rabinowitsch unit certificates on the SPECIALIZED
#         (branch + shape substituted) coefficient systems. Tiny rings.
import sympy as sp, subprocess

r, u0 = sp.symbols('r u0')
al = sp.Symbol('al')
N = sp.sympify(open('q8_layers.srepr').read())

def poly(name, d):
    return sum(sp.Symbol(f'{name}{j}')*r**j for j in range(d+1))
h, c = poly('bh', 2), poly('ch', 2)
K, Q = poly('cc', 3), poly('bc', 3)
Ac = poly('ac', 2)
def dr(f): return sp.diff(f, r)
W = lambda f, g: sp.expand(f*dr(g) - dr(f)*g)
P = sp.expand(2*h + dr(c) - r*dr(h))
T = sp.expand(r*dr(h) - dr(c))

zac = {sp.Symbol(f'ac{j}'): 0 for j in range(3)}
zbc = {sp.Symbol(f'bc{j}'): 0 for j in range(4)}
zcc = {sp.Symbol(f'cc{j}'): 0 for j in range(4)}
zbh = {sp.Symbol(f'bh{j}'): 0 for j in range(3)}

def chk(tag, lhs, rhs):
    print(f'CHECK {tag}:', 'PASS' if sp.expand(lhs - rhs) == 0 else 'FAIL')

# ---- Part 1a: a2 master forms (Ac=0, Bc=0), N* := 2N/al
a2 = {m: sp.expand(N[m].subs(zac).subs(zbc)*2/al) for m in range(2, 6)}
chk('a2 N2*', a2[2], sp.expand(-3*P*h + 2*h*T - 2*T**2 + 2*dr(K) + 2*W(h, c) - (c - r*h)*dr(T)))
chk('a2 N3*', a2[3], sp.expand(-3*P*h*(h - T) + 2*T*(dr(K) + W(h, c)) + 2*W(h, K)
                               + (c - r*h)*dr(P*h) - 2*K*dr(T)))
chk('a2 N4*', a2[4], sp.expand(-3*P*h*(dr(K) + W(h, c)) + 2*T*W(h, K) + 2*K*dr(P*h) - h*K*dr(T)))
chk('a2 N5*', a2[5], sp.expand(h*(K*dr(P*h) - 3*P*W(h, K))))

# ---- Part 1b: b master forms (Bc=0)
b2 = sp.expand(N[2].subs(zbc)*2)
chk('b N2x2', b2, sp.expand(6*Ac - 3*al*P*h + al*(2*(h - T)*T + 2*dr(K) + 2*W(h, c) - (c - r*h)*dr(T))))
b3 = sp.expand(N[3].subs(zbc)*2)
chk('b N3x2', b3, sp.expand(8*Ac*h + 6*(h - T)*Ac - 3*al*(h - T)*P*h + 2*al*T*(dr(K) + W(h, c))
                            + 2*al*W(h, K) - 2*(c - r*h)*dr(Ac) + al*(c - r*h)*dr(P*h) - 2*al*K*dr(T)))
b4 = sp.expand(N[4].subs(zbc))
chk('b N4', b4, sp.expand(4*(h - T)*Ac*h + 3*(dr(K) + W(h, c))*(Ac - al*P*h/2) + al*T*W(h, K)
                          - (c - r*h)*dr(Ac*h) - 2*K*(dr(Ac) - al*dr(P*h)/2) - al*h*K*dr(T)/2))
b5 = sp.expand(N[5].subs(zbc))
chk('b N5', b5, sp.expand(4*(dr(K) + W(h, c))*Ac*h + 3*W(h, K)*(Ac - al*P*h/2)
                          - 2*K*dr(Ac*h) - h*K*(dr(Ac) - al*dr(P*h)/2)))
b6 = sp.expand(N[6].subs(zbc))
chk('b N6', b6, sp.expand(h*(4*Ac*W(h, K) - K*dr(Ac*h))))

# ---- Part 2: Singular unit certificates per leaf
def cert(tag, eqs, vars_, hyps, targets, cap=600):
    """eqs: list of sympy exprs (=0). hyps: nonzero exprs (Rabinowitsch s_i).
    targets: dict name -> expr; each must be UNIT when inverted (t)."""
    coeffs = []
    for e in eqs:
        if sp.expand(e) == 0:
            continue
        for x in sp.Poly(sp.expand(e), r).all_coeffs():
            if x == 0:
                continue
            num, _den = sp.fraction(sp.together(x))  # den is a rational constant
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

albl = [al]
# L-a2-h0 (a2a): Ac=Bc=Bh=0; targets: ch1, ch2 (P = c') must be forced 0
eqs = [N[m].subs(zac).subs(zbc).subs(zbh) for m in range(2, 9)]
cert('a2-h0', eqs, ['ch0', 'ch1', 'ch2', 'cc0', 'cc1', 'cc2', 'cc3', 'al'],
     albl, {'ch1': sp.Symbol('ch1'), 'ch2': sp.Symbol('ch2')})

# L-a2-K0: Ac=Bc=Cc=0; targets: P0 = 2bh0+ch1, P1 = bh1+2ch2 forced 0
eqs = [N[m].subs(zac).subs(zbc).subs(zcc) for m in range(2, 9)]
cert('a2-K0', eqs, ['bh0', 'bh1', 'bh2', 'ch0', 'ch1', 'ch2', 'al'],
     albl, {'P0': sp.Symbol('bh0')*2 + sp.Symbol('ch1'),
            'P1': sp.Symbol('bh1') + 2*sp.Symbol('ch2')})

# L-a2-full (monolithic attempt): Ac=Bc=0; targets P0, P1
eqs = [N[m].subs(zac).subs(zbc) for m in range(2, 9)]
cert('a2-FULL', eqs, ['bh0', 'bh1', 'bh2', 'ch0', 'ch1', 'ch2',
                      'cc0', 'cc1', 'cc2', 'cc3', 'al'],
     albl, {'P0': sp.Symbol('bh0')*2 + sp.Symbol('ch1'),
            'P1': sp.Symbol('bh1') + 2*sp.Symbol('ch2')}, cap=900)
