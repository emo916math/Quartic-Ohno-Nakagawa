test = True;

import re
from sage.interfaces.latte import count
from sage.geometry.polyhedron.cdd_file_format import cdd_Hrepresentation
from sage.interfaces.maxima_lib import maxima_lib, max_to_sr
attach("Fraction.sage")
attach("SubstitutionInGeneratingFunctions.sage")

#########################################################################



#########################################################################


# A function for computing the multivariate generating function for lattice
# points in a region (using LattE) and outputting the result as a Sage
# rational function.
# Output if add==True: ([[num, [dens]]], ring):
# a numerator and a list of denominators to be multiplied, together with
# the ring they are defined over
# Output if add==False: ([[num1, [dens1]], [num2, [dens2]], ...], ring):
# a list of fractions to be added.
def lat_pt_gen_func(cdd_type='rational', ieqs=[], eqns=[],
  variables=None, ring=None, NAMES=None,
  add=True, verb_conversion=False, suppress_warnings=False,
  timing=False, **kwargs):
  
  verb = verb_conversion;
  
  if variables is None:
    variables = ieqs[0].parent().gens();
    
  # The case of zero variables is anomalous.
  if len(variables) == 0:
    for ieq in ieqs:
      assert ieq in QQ, ieq;
      if ieq < 0:
        return Fraction();
    for eqn in eqns:
      assert eqn in QQ, eqn;
      if eqn != 0:
        return Fraction();
    return Fraction(ZZ, 1);
  
  if timing:
    t0 = time.time()
  ieqs = [linear_to_coeffs(rel,variables) for rel in ieqs];
  eqns = [linear_to_coeffs(rel,variables) for rel in eqns];
  numvars = len(variables) if variables else len(ieqs[0]) - 1;
  
  # if verb:
  #   print(ieqs, eqns, "with variables", variables);
  
  if NAMES is None:
    NAMES = [str(v).upper() for v in variables]
  RRR = PolynomialRing(ZZ, names=NAMES)
  QQQ = RRR.fraction_field()
    
  # First detect empty and 1-point regions,
  # which crash some versions of LattE.
  pdn = Polyhedron(ieqs=ieqs, eqns=eqns)
  if pdn.dim() == -1:
    if verb:
      print("Empty region.")
    return Fraction(RRR, 0, [])
  elif pdn.dim() == 0:
    if verb:
      print("1-point region.")
    v = pdn.vertices();
    assert len(v) == 1;
    if verb and variables:
      vtx = v[0].homogeneous_vector();
      assert len(vtx) == numvars + 1;
      print({variables[i] : vtx[i] for i in range(numvars)})
    if v[0].is_integral():
      vtx = v[0].homogeneous_vector();
      assert len(vtx) == numvars + 1;
      assert vtx[numvars] == 1;
      fraction = Fraction(RRR, 
        prod((RRR.gen(i))^vtx[i] for i in range(numvars)), []);
      if verb:
        print(fraction);
      return fraction;
    else:
      return Fraction(RRR, 0, [])
  else:
    if verb:
      print("Using polyhedron", pdn);
    # Use LattE.
    cddin = pdn.cdd_Hrepresentation();
    if verb:
      print("Running LattE...")
    found = False;
    for iter in range(100):
      try:
        fn = count(cddin, cdd=True,
          multivariate_generating_function=True, raw_output=True, **kwargs)
        found = True;
      except IOError:
        # Occurs only when polyhedron is in fact empty.
        # if not(suppress_warnings):
        #   print("Warning: apparently empty polyhedron found with vertices",
        #    pdn.vertices());
        return Fraction(RRR, 0, [])
      except RuntimeError as e:
        print(repr(e));
        print("Shuffling the inequalities and equations");
        ieqs = list(pdn.inequality_generator());
        shuffle(ieqs);
        eqns = list(pdn.equation_generator());
        shuffle(eqns);
        cddin = cdd_Hrepresentation(cdd_type="rational",
          ieqs=ieqs, eqns=eqns);
      if found:
        break;
    if not(found):
      raise RuntimeError("Unable to get an answer for\n" + cddin);
    
    for i in range(len(NAMES)):
      fn = fn.replace("x[" + str(i) + "]", NAMES[i]);
    fn = re.sub(r'\((\-[0-9]+)\)', lambda m:m.group(1), fn)
    
    if verb:
      print("Output length:", len(fn))
      if len(fn) < 300:
        print(fn);
      else:
        print(fn[0:100]);
    
    term_strs = fn.split("\n + ")
    if verb:
      print(len(term_strs), "terms");
    fractions = [];
    if timing:
      td0 = time.time();
    for s in term_strs: # Look to speed this loop.
      nd = s.split("/");
      if len(nd) == 1:
        num = QQQ(nd[0]);
        fractions.append(Fraction(RRR, num, []));
      else:
        num, den = nd;
        num = QQQ(num);
        dens = [QQQ(f.replace("(","").replace(")","")) for f in den.split(")*(")];
        frac = Fraction(RRR, num, dens);
        frac.decomplicate();
        fractions.append(frac);
    if timing:
      print("  Placing LattE output in ring took", time.time() - td0);
    
    if verb:
      print(str(fractions)[0:200])
      print(len(fractions[0].denoms), "denominator binomials per term")
    all_den_factors = {fac for f in fractions for fac in f.denoms}
    if verb:
      print(len(all_den_factors), "distinct denominator binomials")
  
  if add:
    return Fraction.add_all(fractions, verb=verb);
  else:
    return fractions;


# Test for lat_pt_gen_func.
if test:
  # Old syntax, no longer reliably supported.
  
  # frac1 = lat_pt_gen_func(ieqs=[[0,1]], NAMES='T');
  # R1 = frac1.base_ring;
  # num1 = frac1.num;
  # dens1 = frac1.denoms;
  # assert num1 in R1;
  # assert all(den in R1 for den in dens1)
  # assert frac1.evaluate() == 1/(1 - R1.0)
  # 
  # frac2 = lat_pt_gen_func(ieqs=[[2,1,0],[-3,0,1],[1,-7,-5]], NAMES=['A','B'])
  # R2 = frac2.base_ring;
  # A = R2.0; B = R2.1;
  # assert frac2.evaluate() == A^-2 * B^3;
  # 
  # frac3 = lat_pt_gen_func(ieqs=[[3,1],[-4,-1]], NAMES='T')
  # assert frac3.num == 0;
  # 
  # frac4 = lat_pt_gen_func(ieqs=[[0,1,0],[0,0,1]], NAMES=['A','B'])
  # frac5 = lat_pt_gen_func(ieqs=[[0,1,0],[0,0,1],[1,0,0],[2,3,4],[0,5,6],[0,8,7]],
  #   NAMES=['A','B'], redundancy_check="full-cddlib");
  # assert repr(frac4) == repr(frac5);
  
  oldring.<a, b> = PolynomialRing(ZZ);
  frac4 = lat_pt_gen_func(ieqs=[a, b], NAMES=['A','B'])
  frac5 = lat_pt_gen_func(ieqs=[a, b, 1, 2 + 3*a + 4*b, 5*a + 6*b],
    NAMES=['A','B'], redundancy_check="full-cddlib");
  assert repr(frac4) == repr(frac5);
  
  oldring.<a, b, c> = PolynomialRing(QQ);
  frac10 = lat_pt_gen_func(eqns = [a - 4, b - 7, c + 1], variables = [a,b,c],
    add=False)
  R = frac10.base_ring;
  A = R.0; B = R.1; C = R.2;
  assert frac10.evaluate() == A^4*B^7*C^-1
  
  # Test fractional zones
  frac11 = lat_pt_gen_func(ieqs = [a + 1/2, 1/2*b, 1/2*c+1/6, 1/5 - a - b - c])
  assert frac11.evaluate() == 1;
  
  # Test LattE with no variables
  ans1 = lat_pt_gen_func(ieqs = [2 + 0*a, 1/2 + 0*a], eqns = [0*a,0*a], variables = [])
  assert ans1.evaluate() == 1, ans1;
  ans0 = lat_pt_gen_func(ieqs = [2 + 0*a, -1/2 + 0*a], eqns = [0*a,0*a], variables = [])
  assert ans0 == 0, ans0;