/* ==========================================================================
   PARI/GP verification of Alpoge's Jacobian-Conjecture counterexample.
   INDEPENDENT of both sympy and Singular (different CAS, different algorithms:
   determinants by cofactor, elimination by RESULTANTS not Groebner bases,
   discriminant by poldisc, fibers by explicit numeric root-solving).
     F1 = (1+xy)^3 z + y^2 (1+xy)(4+3xy)
     F2 = y + 3x(1+xy)^2 z + 3 x y^2 (4+3xy)
     F3 = 2x - 3 x^2 y - x^3 z
   ========================================================================== */
default(realprecision, 60);

F1 = (1+x*y)^3*z + y^2*(1+x*y)*(4+3*x*y);
F2 = y + 3*x*(1+x*y)^2*z + 3*x*y^2*(4+3*x*y);
F3 = 2*x - 3*x^2*y - x^3*z;

report(tag, cond) = if(cond, print(tag": PASS"), print(tag": FAIL"));

/* ---- CLAIM 1: Jacobian determinant == -2 (cofactor determinant) ---- */
J = [deriv(F1,x),deriv(F1,y),deriv(F1,z); deriv(F2,x),deriv(F2,y),deriv(F2,z); deriv(F3,x),deriv(F3,y),deriv(F3,z)];
dJ = matdet(J);
print("CLAIM1 det(Jacobian) = ", dJ);
report("CLAIM1 det_is_-2", simplify(dJ) == -2);

/* ---- CLAIM 2: fiber over (-1/4,0,0) ---- */
ev(P, X, Y, Z) = subst(subst(subst(P, x, X), y, Y), z, Z);
w1 = [ ev(F1,0,0,-1/4),     ev(F2,0,0,-1/4),     ev(F3,0,0,-1/4) ];
w2 = [ ev(F1,1,-3/2,13/2),  ev(F2,1,-3/2,13/2),  ev(F3,1,-3/2,13/2) ];
w3 = [ ev(F1,-1,3/2,13/2),  ev(F2,-1,3/2,13/2),  ev(F3,-1,3/2,13/2) ];
print("CLAIM2 F(0,0,-1/4)=", w1, " F(1,-3/2,13/2)=", w2, " F(-1,3/2,13/2)=", w3);
tgt = [-1/4, 0, 0];
report("CLAIM2a three_witnesses_map_to_target", (w1==tgt) && (w2==tgt) && (w3==tgt));
/* x-eliminant at target (-1/4,0,0) is -4x^3+4x = -4 x (x-1)(x+1): roots {0,1,-1}
   are EXACTLY the witness x-coordinates. */
Dt = 27*(-1/4)^2*0^2 - 18*(-1/4)*0*0 + 16*(-1/4) + 0^3*0 - 0^2;  /* = -4 */
Ct = Dt*x^3 + (4 - 3*0*0)*x - 2*0;
report("CLAIM2b x-cubic factors as -4 x (x-1)(x+1) [roots {0,1,-1}]", (Ct - (-4*x*(x-1)*(x+1))) == 0);

/* ---- CLAIM 3: generic fiber has exactly 3 points (explicit solve) ----
   For a numeric target, x-coords are the roots of the x-cubic. F3=c is linear
   in z, giving z=(2x-3x^2 y - c)/x^3; substitute into F2=b to get a polynomial
   in y; keep (x,y) whose z also satisfies F1=a. Count distinct genuine triples. */
/* polroots after trimming numerically-zero leading coefficients (floating
   roots can leave a ~1e-50 leading term that only triggers a cosmetic warning) */
cleanroots(P) =
{
  my(v = Vec(P));
  while(#v > 1 && abs(v[1]) < 1e-30, v = v[2..#v]);
  if(#v <= 1, return([]));
  return(polroots(Pol(v)));
}
fibersolve(a0,b0,c0) =
{
  my(D0, Cx, y0, z0, P2y, xr, yr, sols, key);
  D0 = 27*a0^2*c0^2 - 18*a0*b0*c0 + 16*a0 + b0^3*c0 - b0^2;
  Cx = D0*x^3 + (4-3*b0*c0)*x - 2*c0;
  xr = cleanroots(Cx);
  sols = List();
  for(i=1, #xr,
    my(X = xr[i]);
    if(abs(X) < 1e-12, next);              /* guard x=0 (not a root here) */
    /* z as a function of y from F3=c0 */
    my(zofy = (2*X - 3*X^2*y - c0)/X^3);
    P2y = subst(subst(F2, x, X), z, zofy) - b0;   /* polynomial in y */
    yr = cleanroots(P2y);
    for(j=1, #yr,
      my(Y = yr[j], Z);
      Z = (2*X - 3*X^2*Y - c0)/X^3;
      /* verify ALL three equations to weed out spurious roots */
      if( abs(ev(F1,X,Y,Z)-a0) < 1e-9 && abs(ev(F2,X,Y,Z)-b0) < 1e-9
          && abs(ev(F3,X,Y,Z)-c0) < 1e-9,
        key = Str(round(real(X)*10^6), ",", round(imag(X)*10^6), ",",
                  round(real(Y)*10^6), ",", round(imag(Y)*10^6));
        listput(sols, key);
      );
    );
  );
  return(#Set(Vec(sols)));
}
n1 = fibersolve(3/7, -2/5, 1/3);
n2 = fibersolve(11/4, 5/9, -7/2);
n3 = fibersolve(-8/3, 13/11, 1/6);
print("CLAIM3 fiber sizes at 3 generic targets: ", n1, " ", n2, " ", n3);
report("CLAIM3 generic_fiber_count_3", (n1==3) && (n2==3) && (n3==3));

/* ---- CLAIM 4: eliminants via RESULTANTS (independent of Groebner) ----
   Specialize (a,b,c) to random rationals, eliminate z then y by resultants,
   and confirm the claimed x-cubic divides the resulting x-polynomial. */
checkelim(a0,b0,c0) =
{
  my(R13, R23, Gx, D0, Xc);
  R13 = polresultant(F1-a0, F3-c0, z);   /* eliminate z */
  R23 = polresultant(F2-b0, F3-c0, z);
  Gx  = polresultant(R13, R23, y);       /* eliminate y -> poly in x */
  D0  = 27*a0^2*c0^2 - 18*a0*b0*c0 + 16*a0 + b0^3*c0 - b0^2;
  Xc  = D0*x^3 + (4-3*b0*c0)*x - 2*c0;
  return( Gx != 0 && (Gx % Xc) == 0 );   /* claimed cubic divides resultant */
}
e1 = checkelim(3/7, -2/5, 11/13);
e2 = checkelim(5/3, 7/4, -2/9);
print("CLAIM4 resultant-eliminant divisible by claimed x-cubic: ", e1, " ", e2);
report("CLAIM4 x_eliminant_matches_resultant", e1 && e2);

/* ---- CLAIM 5: discriminant of the x-cubic (poldisc) ----
   disc = -4 (27ac^2-9bc+8)^2 D  exactly, over Q(a,b,c). PRIMARY tool = PARI. */
D  = 27*a^2*c^2 - 18*a*b*c + 16*a + b^3*c - b^2;
C  = D*x^3 + (4-3*b*c)*x - 2*c;
disc = poldisc(C, x);
target5 = -4 * (27*a*c^2 - 9*b*c + 8)^2 * D;
print("CLAIM5 poldisc(x-cubic) factored:");
print(factor(disc));
report("CLAIM5 discriminant_square_factor", (disc - target5) == 0);

/* ---- CLAIM 6: special fibers on Gamma via eliminant degeneration ----
   Empty fibers: the x-eliminant specializes to a NONZERO CONSTANT (no x-root).
   Single point: it specializes to a LINEAR polynomial (one x-root).
   (Singular's std(I)==1 / vdim=1 is authoritative; this is the PARI cross-view.) */
xcubic(a0,b0,c0) =
{
  my(D0 = 27*a0^2*c0^2 - 18*a0*b0*c0 + 16*a0 + b0^3*c0 - b0^2);
  return(D0*x^3 + (4-3*b0*c0)*x - 2*c0);
}
g1 = xcubic(4/27,  4/3,  1);
g2 = xcubic(1/27,  2/3,  2);
g3 = xcubic(4/243, -4/9, -3);
print("CLAIM6 empty-target x-eliminants: ", g1, " | ", g2, " | ", g3);
report("CLAIM6a empty_fibers_eliminant_is_nonzero_constant", (poldegree(g1)==0 && g1!=0) && (poldegree(g2)==0 && g2!=0) && (poldegree(g3)==0 && g3!=0));
gs = xcubic(-16/27, 0, 1);
print("CLAIM6 single-point target x-eliminant: ", gs, "  (root x = ", polroots(gs)[1], ")");
report("CLAIM6b single_point_eliminant_is_linear_one_root", poldegree(gs)==1);

/* ---- CLAIM 7: G_m equivariance, weights (-2,-1,+1) ---- */
s1 = subst(subst(subst(F1, x, t*x), y, y/t), z, z/t^2);
s2 = subst(subst(subst(F2, x, t*x), y, y/t), z, z/t^2);
s3 = subst(subst(subst(F3, x, t*x), y, y/t), z, z/t^2);
report("CLAIM7 equivariance_weights_-2_-1_+1", (t^2*s1 - F1 == 0) && (t*s2 - F2 == 0) && (s3 - t*F3 == 0));

/* ---- CLAIM 8: Punt covering identity, r = x/(1+xy) ---- */
r = x/(1+x*y);
punt = 2*F1*r^3 - F2*r^2 + 2*r - F3;
print("CLAIM8 2*F1*r^3 - F2*r^2 + 2*r - F3 simplifies to: ", simplify(punt));
report("CLAIM8 punt_covering_identity", simplify(punt) == 0);

quit;
