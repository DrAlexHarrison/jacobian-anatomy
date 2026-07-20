/* ==========================================================================
   NEGATIVE CONTROLS (PARI/GP). Each feeds a DELIBERATELY BROKEN input to the
   SAME machinery used for a real claim and asserts the property FAILS.
   "PASS" = the break was correctly detected (the property did NOT hold).
   ========================================================================== */
default(realprecision, 60);
report(tag, cond) = if(cond, print(tag": PASS (break detected)"), print(tag": FAIL (suite blind)"));

/* ---- NC1 : det class. Replace (4+3xy) by (4+2xy): det must be nonconstant. ---- */
G1 = (1+x*y)^3*z + y^2*(1+x*y)*(4+2*x*y);
G2 = y + 3*x*(1+x*y)^2*z + 3*x*y^2*(4+2*x*y);
G3 = 2*x - 3*x^2*y - x^3*z;
Jb = [deriv(G1,x),deriv(G1,y),deriv(G1,z); deriv(G2,x),deriv(G2,y),deriv(G2,z); deriv(G3,x),deriv(G3,y),deriv(G3,z)];
dJb = matdet(Jb);
print("NC1 perturbed det = ", dJb);
report("NC1 det_nonconstant_when_broken", (deriv(dJb,x)!=0) || (deriv(dJb,y)!=0) || (deriv(dJb,z)!=0));

/* ---- NC2 : fiber-eval class. Corrupt witness z 13/2 -> 13/3. ---- */
F1 = (1+x*y)^3*z + y^2*(1+x*y)*(4+3*x*y);
F2 = y + 3*x*(1+x*y)^2*z + 3*x*y^2*(4+3*x*y);
F3 = 2*x - 3*x^2*y - x^3*z;
ev(P,X,Y,Z) = subst(subst(subst(P,x,X),y,Y),z,Z);
wbad = [ev(F1,1,-3/2,13/3), ev(F2,1,-3/2,13/3), ev(F3,1,-3/2,13/3)];
print("NC2 F(1,-3/2,13/3) = ", wbad, "  (target -1/4,0,0)");
report("NC2 broken_witness_misses_target", wbad != [-1/4,0,0]);

/* ---- NC4 : eliminant class. Perturbed map's resultant not divisible by claimed cubic. ---- */
a0=3/7; b0=-2/5; c0=11/13;
R13 = polresultant(G1-a0, G3-c0, z);
R23 = polresultant(G2-b0, G3-c0, z);
Gxbad = polresultant(R13, R23, y);
D0 = 27*a0^2*c0^2 - 18*a0*b0*c0 + 16*a0 + b0^3*c0 - b0^2;
Xc = D0*x^3 + (4-3*b0*c0)*x - 2*c0;
print("NC4 (perturbed-map resultant) mod (claimed cubic) = ", Gxbad % Xc);
report("NC4 broken_eliminant_not_divisible", (Gxbad % Xc) != 0);

/* ---- NC5 : discriminant square-factor class. Perturb D -> D+1. ---- */
D = 27*a^2*c^2 - 18*a*b*c + 16*a + b^3*c - b^2;
Cbad = (D+1)*x^3 + (4-3*b*c)*x - 2*c;
discbad = poldisc(Cbad, x);
pattern = -4 * (D+1) * (27*a*c^2 - 9*b*c + 8)^2;
print("NC5 disc(perturbed) - (-4)(D+1)(27ac^2-9bc+8)^2 = ", discbad - pattern);
report("NC5 broken_discriminant_no_square", (discbad - pattern) != 0);

/* ---- NC6 : empty-vs-generic. A generic target's x-eliminant is degree 3, not a constant. ---- */
Dg = 27*(3/7)^2*(1/3)^2 - 18*(3/7)*(-2/5)*(1/3) + 16*(3/7) + (-2/5)^3*(1/3) - (-2/5)^2;
Cg = Dg*x^3 + (4-3*(-2/5)*(1/3))*x - 2*(1/3);
print("NC6 generic-target x-eliminant degree = ", poldegree(Cg), " (want 3, not a constant)");
report("NC6 generic_eliminant_not_constant", poldegree(Cg) == 3);

/* ---- NC7 : equivariance class. Wrong z-weight z -> z*t. ---- */
s1bad = subst(subst(subst(F1,x,t*x),y,y/t),z,z*t);
print("NC7 t^2*F1(tx,y/t,z*t) - F1 = ", t^2*s1bad - F1, " (want nonzero)");
report("NC7 wrong_weight_breaks_equivariance", (t^2*s1bad - F1) != 0);

/* ---- NC8 : punt-identity class. Coefficient 2 -> 3 on F1. ---- */
r = x/(1+x*y);
puntbad = 3*F1*r^3 - F2*r^2 + 2*r - F3;
print("NC8 broken punt (3*F1) simplifies to = ", simplify(puntbad), " (want nonzero)");
report("NC8 broken_punt_identity_nonzero", simplify(puntbad) != 0);

quit;
