#!/usr/bin/env python3
# 2026-07-22: branch-specialized normalized layers for the {734,735}
# hand tree. Loads q8_layers.srepr (N2..N8, Ah = -al*P/2 already substituted).
# Branches from N8 = 5*Ac*Bc*W(Bc,Cc) and N7|_{Ac=0} = -2*al*P*Bc*W(Bc,Cc):
#   a2: Ac=0, Bc=0            (P !== 0 assumed in the hand argument)
#   a3: Ac=0, Cc=0            (after gauge; Bc !== 0, P !== 0)
#   b : Bc=0                  (Ac !== 0)
#   g : Cc=0                  (Ac !== 0, Bc !== 0, after gauge)
# Print each nonzero specialized layer FACTORED.
import sympy as sp

r, u0 = sp.symbols('r u0')
al = sp.Symbol('al')
N = sp.sympify(open('q8_layers.srepr').read())

def sub_zero(expr, names):
    d = {}
    for nm, deg in names:
        for j in range(deg+1):
            d[sp.Symbol(f'{nm}{j}')] = 0
    return sp.expand(expr.subs(d))

BR = {
    'a2 (Ac=0,Bc=0)': [('ac', 2), ('bc', 3)],
    'a3 (Ac=0,Cc=0)': [('ac', 2), ('cc', 3)],
    'b  (Bc=0)':      [('bc', 3)],
    'g  (Cc=0)':      [('cc', 3)],
}
for tag, zs in BR.items():
    print(f'===== BRANCH {tag} =====')
    for m in range(2, 9):
        e = sub_zero(N[m], zs)
        if e == 0:
            print(f'N{m} == 0')
            continue
        f = sp.factor(e)
        s = sp.sstr(f).replace('**', '^')
        print(f'N{m} [deg {sp.degree(e, r)}]: {s if len(s) <= 1500 else s[:1500] + " ...TRUNC"}')
    print()
