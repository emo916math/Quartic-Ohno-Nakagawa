test = True;
import itertools;
from collections import OrderedDict;

# A collection of utilities useful when summing exponential functions
# over lattice points in polyhedral regions. Includes the classes
# ExponentDictionary and RingFactorData.

#########################################################################

# Adds (or performs the indicated operation on) values for the same key.
# Similar to the addition method for Counter, but without
# the positivity checking.
# The given element zero should be at least a right identity for the
# operation (so subtraction is OK).
def dict_add(d1, d2, operation=lambda x,y: x+y, zero=0, in_place=False):
  if in_place:
    ret = d1;
  else:
    ret = copy(d1);
  for (key,val) in d2.items():
    if ret.get(key) is not None:
      newval = operation(ret.get(key), val);
      if newval == zero:
        try:
          ret.pop(key);
        except KeyError:
          pass;
      else:
        ret.update({key : newval});
    else:
      ret.update({key : operation(zero,val)});
  if not(in_place):
    return ret;

# Sums a list (or other iterable) of dictionaries according to the 
# preceding method.
def dict_sum(dicts, operation=lambda x,y: x+y, zero=0):
  ret = {};
  for d in dicts:
    dict_add(ret, d, operation, zero, in_place=True);
  return ret;
  

# Test for these functions
if test:
  d1 = {1 : 1, 2 : 2}
  d2 = {1 : 1, 3 : 3};
  assert dict_add(d1, d2) == {1 : 2, 2 : 2, 3 : 3};
  assert dict_add(d1, d2) == dict_add(d2, d1);
  assert dict_add(d1, d2, lambda x,y:x-y) == {2 : 2, 3 : -3}
  
  assert dict_sum([]) == {};
  assert dict_sum([d1]) == d1;
  assert dict_sum([d1, d2, d1, d1]) == {1 : 4, 2 : 6, 3 : 3};
  
  # Custom operations and zeros.
  assert dict_add(d1, d2, operation=lambda x,y:x*y, zero=1) == {2 : 2, 3 : 3}
  dd1 = {1 : [1], 2 : [2,3]};
  dd2 = {1 : [5], 6 : [4,7]};
  assert dict_add(dd1, dd2, zero = []) == {1 : [1, 5], 2 : [2,3], 6 : [4,7]}
  
  # in_place.
  dict_add(dd1, dd2, zero = [], in_place = True);
  assert dd1 == {1 : [1, 5], 2 : [2,3], 6 : [4,7]}

# A function to rewrite relations given as the nonnegativity of Sage symbolic
# polynomials as lists of coefficients, for input to LattE.
def linear_to_coeffs(rel, variables=None, ring=None):
  if isinstance(rel, list):
    return rel;
  if variables is None:
    if ring is None:
      ring = rel.parent();
    variables = ring.gens()
  else:
    if ring is None:
      ring = coercion_model.common_parent(*([rel] + list(variables)))
     
  ret = [ring(rel).monomial_coefficient(mon) for mon in (ring(1),) + tuple(variables)]
  rel2 = ret[0] + sum([ret[i+1] * variables[i] for i in range(len(variables))]);
  if rel2 != rel:
    raise ValueError(rel-rel2,rel,variables)
  return ret;

# Test for this function
if test:
  Atest.<a> = PolynomialRing(ZZ);
  assert linear_to_coeffs(2*a) == [0,2]
  Atest.<a,b> = PolynomialRing(ZZ);
  assert linear_to_coeffs(3*a+2*b - 1) == [-1,3,2];
  assert linear_to_coeffs(3*b, (b,)) == [0,3];
  try:
    linear_to_coeffs(a*b, (b,));
    assert False;
  except ValueError:
    pass;
  try:
    linear_to_coeffs(a*b, [Atest(2), a, b]);
    assert False;
  except ValueError:
    pass;

class ExponentDictionary:
  # expression const;
  # dict dic;
  def __init__(self, dic={}, const=1):
    self.dic = dic;
    self.const = const;
    
  def get(self, key):
    val = self.dic.get(key);
    return val if val is not None else 0;
  
  def __mul__(self, e):
    if isinstance(e, ExponentDictionary):
      return ExponentDictionary(dict_add(self.dic, e.dic), self.const * e.const);
    else:
      raise TypeError("Trying to multiply " + repr(self) + " by " + repr(e))
      return ExponentDictionary(self.dic.copy(), self.const * e);
  
  def __truediv__(self, e):
    if isinstance(e, ExponentDictionary):
      return ExponentDictionary(dict_add(self.dic, e.dic, lambda x,y: x-y), 
        self.const / e.const);
    return ExponentDictionary(self.dic.copy(), self.const / e);
    
  def __pow__(self, n):
    return ExponentDictionary({key : val*n for (key, val) in self.dic.items()}, self.const^n);
    
  def __str__(self):
    return str(self.dic) + " * (" + str(self.const) + ")";
  def __repr__(self):
    return "ExponentDictionary(" + repr(self.dic) + ", " + repr(self.const) + ")"
  def __eq__(self, other):
    return self.const == other.const and self.dic == other.dic; # limited: fooled by zero exponents
  def __ne__(self, other):
    return not(self == other);
  
  def copy(self):
    return ExponentDictionary(copy(self.dic), self.const);
  
  def normalize(self):
    # print("Normalizing", self);
    for (key, val) in self.dic.items():
      try: 
        cc = val.constant_coefficient();
      except AttributeError:
        if val in QQ:
          cc = val;
        else:
          raise
      if cc != 0:
        self.const *= key^cc;
        self.dic.update({key : val - cc});
  
if test:
  Rtest.<x, y> = PolynomialRing(QQ);
  d1 = ExponentDictionary({x : Rtest(2)}, 4) # 4*x^2
  d2 = ExponentDictionary({y : x}) # y^x
  
  assert d1.get(x) == 2;
  assert d2.get(x) == 0;
  
  assert d1 == eval(repr(d1));
  assert (d1 == d1.copy());
  assert (d1 != d1.copy()) == False;
  assert (d1 != d2) == True;
  
  assert d1 * d2 == ExponentDictionary({x : 2, y : x}, 4);
  assert (d1 * d2) * d2 == d1 * (d2 * d2);
  assert (d1 * d2) / d2 / d1 == ExponentDictionary();
  assert d1^2 == d1*d1;
  
  d1.normalize();
  assert d1 == ExponentDictionary({x : 0}, 4*x^2)
  d1.normalize();
  assert d1 == ExponentDictionary({x : 0}, 4*x^2)

def deep_subs(obj, subs, multiplier=None, verb=False):
  if isinstance(obj, list):
    return [deep_subs(e, subs, multiplier) for e in obj]
  elif isinstance(obj, set):
    return {deep_subs(e, subs, multiplier) for e in obj}
  elif isinstance(obj, tuple):
    return tuple(deep_subs(e, subs, multiplier) for e in obj)
  elif isinstance(obj, Counter):
    return Counter({deep_subs(key, subs, multiplier) : val
      for (key, val) in obj.items()});
  elif isinstance(obj, dict):
    return {key : deep_subs(val, subs, multiplier)
      for (key, val) in obj.items()};
  elif isinstance(obj, ExponentDictionary):
    return ExponentDictionary(
      const = deep_subs(obj.const, subs, multiplier),
      dic = deep_subs(dic, subs)
    );
  else:
    if verb:
      print("Trying to substitute", subs, "into", obj,
        "which lies in", obj.parent());
    try:
      # First try directly. This is useful if obj is an object, such as
      # a Fraction object, which has a "substitute" method
      subs = obj.substitute({key : val for (key,val) in subs.items()});
    except (TypeError, AttributeError):
      # Coerce the input and all the keys to one parent.
      parent1 = reduce(coercion_model.common_parent, subs.keys(), obj.parent());
      # print(parent);
      subs = parent1(obj).substitute(
        {parent1(key) : val for (key,val) in subs.items()});
    if multiplier is None:
      return subs;
    return multiplier * subs;

# A function to substitute in generating functions. Each of NEWVARS appears
# whenever any of the oldvars appears in its exponent and the corresponding
# OLDVAR appears in obj.
def substitute(obj, oldvars, OLDVARS, NEWVARS, exponents, verb=False):
  assert len(oldvars) == len(OLDVARS);
  assert len(NEWVARS) == len(exponents);
  
  allcoeffs = [linear_to_coeffs(exp, oldvars) for exp in exponents];
  multiplier = prod([NEWVARS[i]^allcoeffs[i][0] for i in range(len(NEWVARS))]);
  subs = {OLDVARS[j] : prod([NEWVARS[i]^allcoeffs[i][j+1] for i in range(len(NEWVARS))])
    for j in range(len(oldvars))};
  return deep_subs(obj, subs, multiplier, verb=verb);

def cast(fn,key):
  try:
    return fn.parent()(key);
  except AttributeError:
    return key;
if test:
  r.<a,b> = PolynomialRing(ZZ);
  R0.<A,B> = PolynomialRing(ZZ);
  R = R0.fraction_field();
  fn = A + B;
  assert substitute(fn, [a,b], [A,B], [A,B], [2*a+3*b,a]) == A^3 + A^2*B;
  assert substitute(fn, [a,b], [A,B], [A,B], [2*a+3*b+1,a]) == A*(A^3 + A^2*B);
  fn = 5*A^3;
  assert substitute(fn, [a], [A], [A,B], [2*a, -3*a + 1]) == 5*A^6*B^-8;
  fn = [A, [[A, A], [], {A}]];
  assert substitute(fn, [a], [A], [A,B], [a, 0*a]) == fn;
  
  assert deep_subs(fn, {}) == fn;
  assert deep_subs(1, {}) == 1;
  assert deep_subs(1, {a : 4}) == 1;

# A function to substitute in generating functions. Each of NEWVARS appears whenever any of the oldvars
# appears in its exponent and the corresponding OLDVAR appears in obj.
def expdict_to_subs_multiplier(oldvars, OLDVARS, expdict):
  if len(oldvars) != len(OLDVARS):
    raise ValueError("Length mismatch when keying " + str(oldvars) +
      " by " + str(OLDVARS))
  
  dict_of_coeffs = {NEWVAR : linear_to_coeffs(expdict.get(NEWVAR), oldvars)
    for NEWVAR in expdict.dic};
  multiplier = expdict.const * prod([NEWVAR^dict_of_coeffs.get(NEWVAR)[0] for NEWVAR in 
    dict_of_coeffs]);
  subs = {OLDVARS[j] : prod([NEWVAR^dict_of_coeffs.get(NEWVAR)[j+1] for NEWVAR in dict_of_coeffs])
    for j in range(len(oldvars))};
  return subs, multiplier;

def subs_by_expdict(obj, oldvars, OLDVARS, expdict):
  subs, multiplier = expdict_to_subs_multiplier(oldvars, OLDVARS, expdict);
  return deep_subs(obj, subs, multiplier);
  

if test:
  old.<e> = PolynomialRing(ZZ);
  OLD.<E> = PolynomialRing(ZZ);
  NEW.<E1, N> = PolynomialRing(ZZ);
  d = ExponentDictionary({E1 : e, N : e + 1}, N^100);
  a = E^3;
  assert subs_by_expdict(a, old.gens(), OLD.gens(), d) == E1^3 * N^104;
  d = ExponentDictionary({E1 : e, N : e + 7}, 1);
  assert subs_by_expdict(a, old.gens(), OLD.gens(), d) == E1^3 * N^10;
  d = ExponentDictionary({E1 : e, N : old(4)}, 1);
  assert subs_by_expdict(a, old.gens(), OLD.gens(), d) == E1^3 * N^4;
  
def deep_convert(stuff, ring):
  if isinstance(stuff, list):
    return [deep_convert(e, ring) for e in stuff]
  elif isinstance(stuff, set):
    return {deep_convert(e, ring) for e in stuff}
  elif isinstance(stuff, tuple):
    return tuple(deep_convert(e, ring) for e in stuff)
  elif isinstance(stuff, Counter):
    return Counter({deep_convert(key, ring) : val for (key, val) in stuff.items()});
  elif isinstance(stuff, dict):
    return {key : deep_convert(val, ring) for (key,val) in stuff.items()};
  else:
    return ring(stuff);

if test:
  R1 = QQ; R2 = GF(5);
  assert deep_convert(1/2, R2) == 3;
  assert deep_convert([R1(0),{R1(1),1/2}], R2) == [0, {1,3}];

class RingFactorData:
  def __init__(self, name="",ieqs=[], eqns=[], variables=[], temps={},
      factor=ExponentDictionary(), Ftypes = [], subcases=[]):
    self.name = name;
    self.ieqs = ieqs;
    self.eqns = eqns;
    self.variables = variables;
    self.temps = temps;
    if isinstance(factor, ExponentDictionary):
      self.factor = factor;
    elif isinstance(factor, dict):
      self.factor = ExponentDictionary(dic = factor);
    elif factor.parent().is_ring():
      self.factor = ExponentDictionary(const = factor);
    else:
      raise ValueError(factor)
    self.Ftypes = Ftypes;
    self.subcases = [sw if isinstance(sw, SubcaseSwitch) else SubcaseSwitch(sw)
      for sw in subcases]
  
  def copy(self):
    return RingFactorData(
      name     = self.name     ,
      ieqs     = copy(self.ieqs     ),
      eqns     = copy(self.eqns     ),
      variables= copy(self.variables),
      temps    = copy(self.temps    ),
      factor   = copy(self.factor   ),
      Ftypes   = copy(self.Ftypes   ),
      subcases = deepcopy(self.subcases ),
      );
    
  def dename(self): # A version without the name.
    return RingFactorData(
      name     = "",
      ieqs     = self.ieqs     ,
      eqns     = self.eqns     ,
      variables= self.variables,
      temps    = self.temps    ,
      factor   = self.factor   ,
      Ftypes   = self.Ftypes   ,
      subcases = self.subcases ,
      );
    
  def defactor(self): # A version without the factor.
    return RingFactorData(
      name     = self.name     ,
      ieqs     = self.ieqs     ,
      eqns     = self.eqns     ,
      variables= self.variables,
      temps    = self.temps    ,
      factor   = ExponentDictionary(),
      Ftypes   = self.Ftypes   ,
      subcases = self.subcases ,
      );
  
  def __add__(self, other):
    try:
      other.name;
    except:
      print(other);
      raise;
    if self.name and other.name:
      new_name = self.name + "," + other.name;
    else:
      new_name = self.name + other.name;
    if self.Ftypes and other.Ftypes:
      new_Ftypes = set(self.Ftypes).intersection(set(other.Ftypes))
    else:
      new_Ftypes = self.Ftypes + other.Ftypes;
    
    conflict = [];
    new_temps = self.temps.copy();
    new_temps.update(other.temps);
    for key in set(self.temps).intersection(set(other.temps)):
      eqn = self.temps.get(key) - other.temps.get(key);
      if eqn != 0:
        conflict += [eqn];
        if False:
          print("Consolidating", {key : self.temps.get(key)}, "and",
            {key : other.temps.get(key)}, "into equation 0 = ", eqn) 
            
    
    return RingFactorData(new_name, self.ieqs + other.ieqs,
      conflict + self.eqns + other.eqns,
      deduplicate(self.variables + other.variables),
      new_temps,
      self.factor * other.factor,
      new_Ftypes,
      other.subcases + self.subcases);
    # reverse order on subcases: appears to speed up traversal.
  
  def __repr__(self):
    return ("RingFactorData(" + repr(self.name)
      + ", ieqs=" + str(self.ieqs)
      + ", eqns=" + str(self.eqns)
      + ", variables=" + str(self.variables)
      + ", temps=" + str(self.temps)
      + ", factor=" + repr(self.factor)
      + ", Ftypes=" + repr(self.Ftypes)
      + ", subcases:" + repr(self.subcases) 
      + ")")
    
  def __str__(self):
    return repr(self);
    
  def __eq__(self, other):
    return self.name == other.name and\
        self.ieqs == other.ieqs and\
        self.eqns == other.eqns and\
        set(self.variables) == set(other.variables) and\
        self.temps == other.temps and\
        self.factor == other.factor and\
        set(self.Ftypes) == set(other.Ftypes) and\
        self.subcases == other.subcases;
    
  # Replaces variables by values in ieqs and eqns, also in factor if the
  # optional parameter factor_too is set to True. Does not touch name, Ftypes, subcases.
  def substitute(self, dic, factor_too = False):
    newvars = copy(self.variables);
    for v in dic.keys():
      if v in newvars:
        newvars.pop(newvars.index(v));
    return RingFactorData(
      name = self.name,
      ieqs = deep_subs(self.ieqs, dic),
      eqns = deep_subs(self.eqns, dic),
      variables = newvars,
      temps = deep_subs(self.temps, dic),
      factor = ExponentDictionary(const = self.factor.const,
        dic = deep_subs(self.factor.dic, dic)) if factor_too else self.factor,
      Ftypes = self.Ftypes,
      subcases = self.subcases
    );
  
  # Returns true if self is trivially empty because it includes a condition
  # of one of the types 1 = 0, 1 <= 0, or 2x = 1.
  def is_trivially_empty(self):
    for eqn in self.eqns:
      if eqn in QQ and eqn != 0: # 1 = 0
        return True;
      if eqn.content() != (eqn - eqn.constant_coefficient()).content(): # 2x = 1
        return True;
    for ieq in self.ieqs:
      if ieq in QQ and ieq < 0: # 1 <= 0
        return True;
    return False;
  
  # Returns true if self is an empty zone, that is,
  # if no values of the variables can satisfy all the inequalities and equations,
  # even ignoring the condition that they be integral.
  def is_empty(self):
    if self.is_trivially_empty():
      return True;
    ieqs_coeffs = [];
    for rel in self.ieqs:
      try:
        ieqs_coeffs += [linear_to_coeffs(rel, self.variables)];
      except ValueError:
        pass; # At intermediate nodes, ignore ieqs with foreign variables.
    eqns_coeffs = [];
    for rel in self.eqns:
      try:
        eqns_coeffs += [linear_to_coeffs(rel, self.variables)];
      except ValueError:
        pass; # At intermediate nodes, ignore eqns with foreign variables.
    pdn = Polyhedron(ieqs=ieqs_coeffs, eqns=eqns_coeffs);
    # return pdn #$ Debugging purposes only
    return pdn.dim() == -1;
    
  # Eliminates a variable from self.
  # If self does not use the variable var, this is a no-op.
  # Otherwise, self must have either
  #   - an equation that uniquely determines the value of var, or 
  #   - a value specified in its factor for the variable repl_key, which is then used.
  def elim(self, var, repl_key=None):
    if not(var in self.variables):
      return;
    value = None;
    # Find an equation containing var alone:                                          
    for eqn in self.eqns:
      a = eqn.coefficient(var);
      if a != 0 and a*var - eqn in ZZ:
        value = ZZ((a*var - eqn)/a);
        break;
    if value is None:
      value = self.factor.dic.get(repl_key);
    if value is None:
      raise ValueError("No specification for {} in {}".format(var, self));
    
    repl = {var : value};
    self.ieqs = deep_subs(self.ieqs, repl);
    self.eqns = deep_subs(self.eqns, repl);
    self.variables.remove(var);
    self.factor.dic = deep_subs(self.factor.dic, repl);
    
    return value;
  
  def elim_temps(self):
    if not(self.temps):
      return self;
    ret = copy(self);
    ret.temps = {};
    for i in range(len(self.temps) + 1): # The number of variables is a bound
      ret2 = ret.substitute(self.temps, factor_too=True);
      if ret2 == ret:
        return ret;
      ret = ret2;
    print("Warning: circular substitution in", self);
    return ret;
  
  def has_subcases(self):
    return len(self.subcases) > 0;
  
  def successors(self):
    keep = copy(self);
    # Remove the highest-priority subcase switch, the first one in the list in
    # the case of a tie.
    keep.subcases = copy(keep.subcases);
    sw = keep.subcases.pop(keep.subcases.index(
      max(keep.subcases, key = lambda sc:sc.priority)))
    return [keep + case for case in sw.switch];
  
  # Find a subzone with a specified name.
  def fish(self, name=""):
    if self.name == name:
      return self;
    if self.has_subcases():
      for case in self.successors():
        if name.startswith(case.name):
          f = case.fish(name);
          if f is not None:
            return f;
  
  def fish_all(self, name=""):
    ret = [];
    if self.name == name:
      ret += [self];
    if self.has_subcases():
      for case in self.successors():
        if name.startswith(case.name):
          f = case.fish_all(name);
          ret += f;
    return ret;

def deduplicate(list_):
  return list(OrderedDict.fromkeys(list_))

class SubcaseSwitch:
  def __init__(self, switch, priority = 0):
    assert priority in RR, priority;
    self.priority = priority;
    assert all(isinstance(sc, RingFactorData) for sc in switch)
    self.switch = switch;
    
  def __repr__(self):
    return repr({self.priority : [sc.name for sc in self.switch]})

if test:
  tester.<a,b,c> = PolynomialRing(ZZ);
  TESTER.<A,B,C> = PolynomialRing(ZZ);
  Z1 = RingFactorData(
    "angle",
    ieqs = [a, b-a],
    variables = [a,b],
    factor = ExponentDictionary({B : b})
  );
  Z2 = RingFactorData(
    "line",
    eqns = [b - 2],
    factor = ExponentDictionary(const = 7)
  );
  
  Z12 = Z1 + Z2;
  assert Z12.name == "angle,line";
  assert Z12.ieqs == Z1.ieqs;
  assert Z12.eqns == Z2.eqns;
  assert Z12.variables == [a,b];
  assert Z12.factor == ExponentDictionary({B : b}, const = 7);
  
  assert not(Z12.is_trivially_empty());
  assert not(Z12.is_empty());
  
  assert Z12.elim(b) == 2;
  assert Z12 == RingFactorData('angle,line',
    ieqs=[a, -a + 2],
    eqns=[0],
    variables=[a],
    factor=ExponentDictionary(const=7, dic={B:2}),
    Ftypes=[]
  )
  
  assert Z12.elim(a,B) == 2;
  assert Z12 == RingFactorData('angle,line',
    ieqs=[2, 0],
    eqns=[0],
    variables=[],
    factor=ExponentDictionary(const=7, dic={B:2}),
    Ftypes=[]
  );
  
  Z1temp = RingFactorData(
    "angle",
    ieqs = [a + c - 3, b-a],
    variables = [a,b,c],
    temps = {c : 3},
    factor = ExponentDictionary({B : b})
  );