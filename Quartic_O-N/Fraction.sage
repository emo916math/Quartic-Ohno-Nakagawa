from collections import Counter
import time
attach("SubstitutionInGeneratingFunctions.sage")

test = True;

class Fraction:
  # Instance variables:
  # base_ring
  # num
  # denoms
  
  # Parameters:
  #   base_ring - a base ring
  #   num - the numerator
  #   denoms - either a list of denominators, or a counter or dictionary,
  #   each entry d:e meaning d^e.
  def __init__(self, base_ring=ZZ, num=0, denoms=[]):
    if not isinstance(base_ring, Ring):
      raise ValueError(base_ring);  
    self.base_ring = base_ring;
    
    self.num = self.normalize(num);
    denoms0 = Counter(denoms);
    self.denoms = Counter();
    for (den,exp) in denoms0.items():
      den = self.normalize(den);
      if den == 0:
        raise ZeroDivisionError();
      self.denoms.update({den : exp});
  
  def normalize(self, value):
    if value in self.base_ring:
      return self.base_ring(value);
    if value in self.base_ring.fraction_field():
      return self.base_ring.fraction_field()(value);
    raise ValueError("Value " + str(value) + " lies in " +
          str(value.parent()) + " and does not lie in the fraction field of " +
          str(self.base_ring))
  
  def __str__(self):
    try:
      abbreviate = self.num.hamming_weight() > 30;
      if abbreviate:
        ret = "(" + str(self.num.hamming_weight()) + " terms)"
      else:
        ret = "(" + str(self.num) + ")";
    except AttributeError:
      numstr = str(self.num);
      abbreviate = len(numstr) > 200;
      if abbreviate:
        ret = "(" + numstr[:10] + "..." + numstr[-10:] + ")"
      else:
        ret = "(" + numstr + ")"
    ret += " / ";
    ret += "(" + str(Factorization(self.denoms.items(),
      sort=False, simplify=False)) + ")";
    return ret;
  
  def __repr__(self):
    # return "Fraction(" + repr(self.base_ring) + "," +\
    #   repr(self.num) + "," + repr(self.denoms) + ")";
    return str(self);
    
  def __eq__(self, other):
    if other == 0: # can be one level recursive
      return self.num == 0;
    try:
      return self.num == other.num and self.denoms == other.denoms;
    except:
      return super().__eq__(other)
    
  def evaluate(self):
    base_field = self.base_ring.fraction_field();
    return base_field(self.num) / prod(
      [base_field(den)^exp for (den,exp) in self.denoms.items()]);

  def decomplicate(self, factorize=False):
    base_field = self.base_ring.fraction_field();
    
    # Deal with the numerator.
    dens_new = Counter();
    if self.num in self.base_ring:
      num_new = self.base_ring(self.num);
      assert num_new.parent() is self.base_ring;
    else:
      num_new = self.base_ring(self.num.numerator());
      d = self.base_ring(self.num.denominator());
      if factorize:
        fzn = factor(d);
        num_new //= fzn.unit();
        dens_new += Counter(dict(fzn));
      else:
        unit = ZZ(sgn(d.lc()));
        num_new //= unit;
        dens_new += Counter([d*unit]);
        
    # Deal with the denominators.
    for (den, exp) in self.denoms.items():
      if den in self.base_ring:
        den_new = self.base_ring(den);
      else:
        mult_by = self.base_ring(den.denominator()/sgn(den.numerator().lc()));
        num_new *= self.base_ring(mult_by^exp);
        den_new = self.base_ring(den*mult_by);
      
      if factorize:
        try:
          fzn = factor(den_new);
          num_new *= self.base_ring(fzn.unit()^(-exp));
          factorized_successfully = True;
        except(NotImplementedError, TypeError):
          factorized_successfully = False;
      else:
        factorized_successfully = False;
       
      if factorized_successfully:
        dens_new += Counter(dict(fzn^exp))
      else:
        unit = ZZ(sgn(den_new.lc()));
        num_new //= unit^exp;
        dens_new += Counter({den_new*unit : exp})
   
    # modify the fraction in place
    self.num = num_new;
    assert num_new.parent() is self.base_ring, self;
    self.denoms = dens_new;
    assert all(den.parent() is self.base_ring for den in self.denoms), self
    
  # Used to check if a polynomial is divisible by a binomial in a quick,
  # probabilistically accurate way. Monomial divisors are ignored.
  @staticmethod
  def divisibility_tester(base_ring,
      divisor, p=next_prime(10^5)):
    if divisor not in base_ring:
      raise ValueError("Domain mismatch: cannot divide by " + str(divisor))
    if divisor.hamming_weight() != 2:
      raise ValueError("Non-binomial divisor: " + str(divisor))
    exp1, exp2 = divisor.exponents();
    term1 = base_ring({exp1 : divisor.dict().get(exp1)});
    term2 = base_ring({exp2 : divisor.dict().get(exp2)});
    assert term1 + term2 == divisor;
    try:
      exps = list(exp1.esub(exp2));
    except AttributeError:
      assert base_ring.ngens() == 1;
      # The output format is different for 1-variable poly rings.
      exps = [exp1 - exp2];
    g, adjust_exps = extended_gcd(exps);
    if g != 1:
      raise UnsuitableDenominatorError(
        "Nontrivial power in divisor: " + str(divisor))
    
    F = GF(p);
    success = False;
    while not(success):
      try:
        subs = [F.random_element() for g in base_ring.gens()];
        ratio = -term1(subs) / term2(subs);
        for i in range(base_ring.ngens()):
          subs[i] *= ratio^(-adjust_exps[i]);
        assert term1(subs) == -term2(subs);
        assert divisor(subs) == 0;
        success = True;
      except ZeroDivisionError:
        pass;
    return subs;
  
  @staticmethod
  def quick_is_divisible_by_binomial(base_ring,
      dividend, divisor, p=next_prime(10^5)):
    subs = Fraction.divisibility_tester(base_ring, divisor, p);
    return dividend(subs) == 0;
  
  @staticmethod
  def is_divisible_by_binomial(base_ring, dividend, divisor):
    try:
      if not(Fraction.quick_is_divisible_by_binomial(
        base_ring, dividend, divisor, 1009)):
          return false;
      if not(Fraction.quick_is_divisible_by_binomial(
        base_ring, dividend, divisor, next_prime(10^7))):
          return false;
    except ValueError:
      pass;
    return dividend % divisor == 0;
  
  # Used to cancel a denominator that divides the numerator.
  # If the optional parameter candidates is not provided, all denom factors are
  # tried.
  def simplify(self, candidates=None, factorize=False, verb="auto"):
    self.decomplicate(factorize=factorize);
    
    if verb == "auto":
      try:
        verb = self.num.hamming_weight() > 1000;
      except AttributeError:
        verb = False;
    if verb:
      print(self)
    num_new = self.num; dens_new = self.denoms;
    if candidates is None:
      candidates = list(reversed(list(self.denoms.keys())));
      # reverse order to take into account the
      # strange pattern that the last denom is more likely to give
      # a cancellation
    i = 0;
    for den in candidates:
      den = self.base_ring(den);
      i += 1;
      if verb:
        print("Trying denominator", i, "of", len(candidates), end=' ') 
        sys.stdout.flush();
      while dens_new.get(den) is not None:
        if verb:
          print("%", end=' ');
          sys.stdout.flush();
        try:
          if verb and den.hamming_weight() == 2:
            divisible = Fraction.is_divisible_by_binomial(
              self.base_ring, self.base_ring(num_new), self.base_ring(den));
          else:
            if verb and den.hamming_weight() >= 3:
              print("<multinomial>", end=' ');
              sys.stdout.flush();
            divisible = (num_new % den == 0);
        except(ArithmeticError, AttributeError, UnsuitableDenominatorError):
          divisible = num_new / den in self.base_ring;
        if divisible:
          if verb:
            print("//", end=' ');
            sys.stdout.flush();
          divide = num_new // den; # truncated division
          # assert den * divide == num_new;
          num_new = self.base_ring(divide);
          if verb:
            print(num_new.hamming_weight(), "terms now");
          dens_new -= Counter([den]);
        else:
          if verb:
            print("/\\");
          break;
    
    # Modify the fraction in place.
    self.num = num_new;
    assert num_new.parent() is self.base_ring;
    self.denoms = dens_new;
    assert all(den.parent() is self.base_ring for den in self.denoms)
    
  def add(self, other, verb=True):
    if verb:
      t0 = time.time()
      
    new_ring = coercion_model.common_parent(self.base_ring, other.base_ring)
    num1 = self.num;
    num2 = other.num;
    dens1ctr = self.denoms;
    dens2ctr = other.denoms;
    dens_new = dens1ctr | dens2ctr;
    num_new = prod((den^cnt for (den, cnt) in
      (dens2ctr - dens1ctr).items()), num1) +\
      prod((den^cnt for (den, cnt) in
      (dens1ctr - dens2ctr).items()), num2)
    if verb:
      t1 = time.time()
      print("Addition time:", t1 - t0, end = " ");
      sys.stdout.flush()
    ret = Fraction(new_ring, num_new, dens_new);
    ret.simplify(dens1ctr & dens2ctr, verb = "auto" if verb else False)
    if verb:
      print("Simplification time:", time.time() - t1)
    return ret;
  
  def __add__(self, other):
    return self.add(other, verb=False);
    
  def adding_complexity(self, other):
    dens1ctr = self.denoms;
    dens2ctr = other.denoms;
    dens_new = dens1ctr | dens2ctr;
    return len(dens_new)
    # return sum(den.degree() * exp for (den,exp) in dens_new.items())
    
  def scalar_mult(self, factor):
    return Fraction(self.base_ring, self.num * factor, self.denoms);
  
  def __sub__(self, other):
    return self + scalar_mult(other, -1);
    
  # A method for adding a long list of fractions,
  # searching to avoid overly large intermediate numerators
  # and denominators.
  
  # Note: Adds all fractions in place.
  @staticmethod
  def add_all(fractions, verb=False, min_complexity = 10, max_complexity = 20):
    if verb:
      print("Adding all of", fractions);
      t0 = time.time();
    if len(fractions) == 0:
      return Fraction();
    complexity_record = 0;
    for bound in [min_complexity..max_complexity]:
      if verb:
        print("Trying complexity", bound);
      # Run once over all the fraction pairs.
      i = 0;
      while i < len(fractions) - 1:
        j = i+1;
        while j < len(fractions):
          c = fractions[i].adding_complexity(fractions[j])
          if c <= bound:
            if verb:
              print("Adding", i, "and", j, "with complexity", c)
            complexity_record = max(complexity_record, c);
            newfrac = fractions[i].add(fractions[j], verb);
            fractions[i] = newfrac;
            fractions.pop(j);
            if verb:
              print(len(fractions), "left")
          else:
            j += 1;
          if len(fractions) <= 16:
            bound = max_complexity + 1; # i.e. break
        i += 1;
          
    while len(fractions) > 1:
      ibest = 0; jbest = 1; cbest = Infinity;
      for i in [0..len(fractions) - 2]:
        for j in [i+1..len(fractions) - 1]:
          c = fractions[i].adding_complexity(fractions[j])
          if c < cbest:
            ibest = i; jbest = j; cbest = c;
      i = ibest; j = jbest; c = cbest;
      if verb:
        print("Adding", i, "and", j, "with complexity", c)
      complexity_record = max(complexity_record, c);
      newfrac = fractions[i].add(fractions[j], verb);
      fractions[i] = newfrac;
      fractions.pop(j);
      if verb:
        print(len(fractions), "left")
      
    assert len(fractions) == 1;
    if verb:
      print("Total time:", time.time() - t0)
    if verb:
      print("Maximum complexity was", complexity_record)
    return fractions[0];
  
  # The crossing space of a fraction measures where it takes its highest values
  def crossing_space(self):
    # No funny business.
    self.decomplicate();
    
    base_ring = self.base_ring;
    if base_ring.ngens() <= 1:
      raise NotImplementedError("crossing space of a 1-variable fraction")
    crossing_vecs = [];
    for divisor in self.denoms:
      if divisor.hamming_weight() != 2: # Hamming weights >= 3 arise from fzn
        continue;
      # if divisor.hamming_weight() > 2:
      #   raise ValueError("Non-binomial divisor: " + str(divisor))
      exp1, exp2 = divisor.exponents();
      term1 = base_ring({exp1 : divisor.dict().get(exp1)});
      term2 = base_ring({exp2 : divisor.dict().get(exp2)});
      assert term1 + term2 == divisor;
      exps = list(exp1.esub(exp2));
      crossing_vecs += [exps];
    V = QQ^(base_ring.ngens());
    ret = V.subspace(crossing_vecs)
    return ret;
    
  @staticmethod
  def add_all_new(fractions, verb=False, min_complexity=10, max_complexity=20):
    
    if verb:
      print(len(fractions), "fractions")
    if len(fractions) == 0:
      return Fraction()
    base_ring = fractions[0].base_ring;
    
    # Hash the fractions, first by dimension of the crossing space, then by the
    # crossing space itself.
    fracs_by_crossing_space = [{} for dim in range(base_ring.ngens() + 1)];
    unmatched = [];
    for frac in fractions:
      try:
        sp = frac.crossing_space();
      except:
        # Crossing space not implemented for some reason.
        # Turn off verbosity, because likely fracs are very simple.
        return Fraction.add_all(fractions, False, min_complexity, max_complexity)
      i = sp.dimension();
      fracs_by_crossing_space[i] = dict_add(
        fracs_by_crossing_space[i], {sp : [frac]}, zero=[])
    if verb:
      for i in reversed(range(base_ring.ngens() + 1)):
        if len(fracs_by_crossing_space[i]) > 0:
          print("In dimension", i, "there are",
            len(fracs_by_crossing_space[i]), "crossing spaces, with",
            [len(val) for (sp,val) in fracs_by_crossing_space[i].items()],
            "terms")
    def check_sanity():
      assert len(fracs_by_crossing_space) == base_ring.ngens() + 1,\
        (len(fracs_by_crossing_space), base_ring.ngens() + 1);
      for i in range(len(fracs_by_crossing_space)):
        for (sp, fracs) in fracs_by_crossing_space[i].items():
          for frac in fracs:
            assert frac.crossing_space() == sp, (i, sp, frac);
    def combine(fracs_sp, i, j):
      assert i < j;
      if verb:
        cplxty = fracs_sp[i].adding_complexity(fracs_sp[j])
        print("Adding", i, "and", j, "with complexity", cplxty);
      newfrac = fracs_sp[i].add(fracs_sp[j], verb);
      sp_new = newfrac.crossing_space();
      if sp_new == sp:
        fracs_sp[i] = newfrac;
        fracs_sp.pop(j);
        ret = 0;
      else:
        oldfrac1 = fracs_sp.pop(j);
        oldfrac2 = fracs_sp.pop(i);
        d_new = sp_new.dimension();
        assert d_new < d, (oldfrac1, oldfrac2, newfrac, sp, sp_new);
        fracs_by_crossing_space[d_new] = dict_add(
          fracs_by_crossing_space[d_new], {sp_new : [newfrac]}, zero=[])
        ret = 1;
        if verb:
          print("%%%%%%%%% Dropped to dimension", d_new)
      if verb:
        print(len(fracs_sp), "left")
      return ret;
    
    complexity_record = 0;
    for d in reversed(range(base_ring.ngens() + 1)):
      check_sanity();
      if verb:
        print("############ Summing crossing spaces of dimension", d)
      while len(fracs_by_crossing_space[d]) > 0:
        sp, fracs_sp = fracs_by_crossing_space[d].popitem();
        if verb:
          print("########## Crossing space:\n", sp)
        while len(fracs_sp) > 0:
          for bound in [min_complexity..max_complexity]:
            if verb:
              print("Trying complexity", bound);
              print(len(fracs_sp), "left")
            # Run once over all the fraction pairs for each complexity.
            i = 0;
            while i < len(fracs_sp) - 1:
              j = i+1;
              while j < len(fracs_sp):
                c = fracs_sp[i].adding_complexity(fracs_sp[j])
                if c <= bound:
                  complexity_record = max(complexity_record, c);
                  dropped = combine(fracs_sp, i, j);
                  if dropped:
                    j = i+1;
                else:
                  j += 1;
                if len(fracs_sp) <= 16:
                  break;
              i += 1;
              if len(fracs_sp) <= 16:
                break;
            if len(fracs_sp) <= 16:
                break;
              
          while len(fracs_sp) > 1:
            ibest = 0; jbest = 1; cbest = Infinity;
            for i in [0..len(fracs_sp) - 2]:
              for j in [i+1..len(fracs_sp) - 1]:
                c = fracs_sp[i].adding_complexity(fracs_sp[j])
                if c < cbest:
                  ibest = i; jbest = j; cbest = c;
            i = ibest; j = jbest; c = cbest;
            complexity_record = max(complexity_record, c);
            dropped = combine(fracs_sp, i, j);
            
          if len(fracs_sp) == 1:
            # Check if fraction is in fact hollow.
            frac = fracs_sp.pop();
            try:
              hollow = frac.is_hollow();
            except NotImplementedError:
              hollow = False;
            if hollow:
              if verb:
                print("Hollow")
              dec = frac.hollow_decomp();
              if all(term.crossing_space().dimension() < d for term in dec):
                if verb:
                  print("Decomposing into dimensions", end=" ")
                for term in dec:
                  spt = term.crossing_space();
                  dt = spt.dimension();
                  if verb:
                    print(dt, end=" ");
                  fracs_by_crossing_space[dt] = dict_add(
                    fracs_by_crossing_space[dt], {spt : [term]}, zero=[])
                if verb:
                  print();
              else:
                unmatched.append(frac);
            else:
              unmatched.append(frac);
    if verb:
      print("There were", len(unmatched), "unmatched fractions");
    final_sum = Fraction.add_all(unmatched, verb=verb);
    if verb:
      print("Maximum complexity in phase I was", complexity_record)
    return final_sum;
  
  def substitute(self, vals, target_ring):
    self.decomplicate();
    if isinstance(vals, dict):
      return Fraction(target_ring, 
        self.num.substitute(vals),
        dict_sum({den.substitute(vals) : exp}
          for (den, exp) in self.denoms.items()))
    else:
      return Fraction(target_ring, 
        self.num(vals),
        dict_sum({den(vals) : exp}
          for (den, exp) in self.denoms.items()))
  
  def is_hollow(self):
    dens_ideal = self.base_ring.ideal([den for den in self.denoms if den not in QQ])
    Rnew = coercion_model.common_parent(self.base_ring, QQ)
    return self.num in Rnew * dens_ideal;
  
  def hollow_decomp(self, verb = False):
    dens_ideal = self.base_ring.ideal([den for den in self.denoms if den not in QQ])
    num_Q = self.num / 1;
    try:
      if verb:
        print("(", end='', flush=True)
      coefs = num_Q.lift(dens_ideal);
      if verb:
        print(")", end='', flush=True)
    except ValueError:
      return None; # fraction is not hollow
    ret = [];
    for (coef, denom) in zip(coefs, self.denoms):
      try:
        term = Fraction(self.base_ring, coef,
          dict_add(self.denoms, {denom : -1}));
      except ValueError:
        term = Fraction(self.base_ring, coef.numerator(),
          dict_sum([self.denoms, {denom : -1}, {coef.denominator() : 1}]));
      term.simplify(verb=verb);
      ret += [term];
    return ret;
    
  def full_hollow_decomp(self, verb=False):
    if self.is_hollow():
      return sum((term.full_hollow_decomp(verb=verb)
        for term in self.hollow_decomp(verb=verb)), []);
    else:
      return [self];
  
# Helper method:
# Given a list [a1,...,ar], returns their gcd, g, and a list of integers
# x1,...,xr such that a1*x1 + ... + ar*xr = g.
def extended_gcd(nums):
  if len(nums) == 0:
    return 0, [];
  if len(nums) == 1:
    num, = nums;
    return abs(num), [sgn(num)];
  g1, x1 = extended_gcd(nums[0:len(nums) - 1]);
  g, a, b = xgcd(nums[len(nums) - 1], g1);
  return g, [b * x for x in x1] + [a]

def add_type(totals, Ftype, shuffle_=False, verb=True):
  tot = totals.get(Ftype);
  if shuffle_:
    shuffle(tot)
  return Fraction.add_all_new(tot, verb=verb);
  
def add_all_types(totals, shuffle_=False, verb=True):
  t0 = time.time();
  answers = {};
  for Ftype in totals.keys():
    print("Beginning type", Ftype);
    answers.update({Ftype : add_type(totals, Ftype, shuffle_=shuffle_,
        verb=verb)});
  print("Total addition time:", time.time() - t0)
  return answers;

# Testing for the foregoing class.
if test:
# For evaluate() and decomplicate()
  Btest.<x,y> = PolynomialRing(ZZ);
  Qtest = Btest.fraction_field();
  exs = [
    (1, [x]),
    (1, [1/x]),
    (y/(x+1),
    Counter({-x/y : 3})),
    (Qtest(1), [x]),
    (1, [Qtest(x)]),
    (Qtest(1), [Qtest(x)]),
    (x, {x:3}),
    (x, [x^3])
  ]
  for data in exs:
    num, denoms = data;
    frac = Fraction(Btest, num, denoms);
    ev1 = frac.evaluate();
    frac.decomplicate();
    ev2 = frac.evaluate();
    frac.decomplicate(factorize=True);
    ev3 = frac.evaluate();
    assert ev1 == ev2 and ev2 == ev3;
    
# For extended_gcd()
  for nums in [[], [-3], [4, -7], [1, 3, 0], [0, 0, 1], [6, 10, 15]]:
    g, x = extended_gcd(nums);
    assert g in ZZ;
    assert g == gcd(nums)
    assert len(x) == len(nums);
    assert (a in ZZ for a in x);
    assert sum(a*b for (a,b) in zip(nums,x)) == g;

# For quick_is_divisible_by_binomial();
  TR.<x,y,z> = PolynomialRing(ZZ);
  assert Fraction.quick_is_divisible_by_binomial(TR, x^2, x - 1) == False;
  assert Fraction.quick_is_divisible_by_binomial(TR, x^2 - 1, x - 1) == True;
  divisor = 5*x^4 - 7*y^2*z^3;
  tp = TR.random_element();
  assert Fraction.quick_is_divisible_by_binomial(TR, tp*divisor, divisor) == True;
  assert Fraction.quick_is_divisible_by_binomial(TR, tp*divisor + x^2, divisor) == False;

# For __add__ and simplify().
  Atest.<z> = PolynomialRing(ZZ);
  Btest.<x,y> = PolynomialRing(ZZ);
  for data in [
    [Atest, (1, [z,z-1]), (1, [z-1,z-2])], 
    [Atest, (1, [z,z-1]), (3, [z-1,z-2])],
    [Btest, (x, {x-y}), (y, {x-y})],
    [Btest, (2*x, {x-y}), (-2*y, {x-y})],
    [Btest, (2*x, {x-y}), (-2*y, {x-y})],
    [Btest, (1, [x,x]), (-1,[x,x])],
    [Btest, (1, [x,x]), (x-1,[x,x])]
  ]:
    base_ring, fr1, fr2 = data;
    num1, dens1 = fr1;
    num2, dens2 = fr2;
    frac1 = Fraction(base_ring, num1, dens1);
    frac2 = Fraction(base_ring, num2, dens2);
    frac_sum = frac1 + frac2;
    assert frac1.evaluate() + frac2.evaluate() == frac_sum.evaluate();
    for den in frac_sum.denoms:
      assert not(den.divides(frac_sum.num));

# For __eq__.
  assert Fraction() == Fraction(QQ, 0, [])
  assert Fraction(ZZ, 1, [3]) == Fraction(ZZ, 1, [3]);
  assert Fraction(ZZ, 1, [3]) != Fraction(ZZ, 1, [4]);
  assert Fraction() == 0
  assert 0 == Fraction()

# For add_all().
  fracs = [];
  Btest.<x,y> = PolynomialRing(ZZ);
  assert(Fraction.add_all(
    [Fraction(Btest, 1, [x^2 - x, x^2 - y]),
    Fraction(Btest, -1, [x - y, x^2 - x]),
    Fraction(Btest, 1, [x - y, x^2 - y])])).evaluate() == 0
  for i in range(20):
    num = Btest.random_element();
    den = Btest.random_element();
    try:
      fracs += [Fraction(Btest, num, [den])];
      fracs += [Fraction(Btest, -num, [den])];
    except ZeroDivisionError:
      pass;
  shuffle(fracs);
  total = Fraction.add_all(fracs, min_complexity = 1, max_complexity = 1);
  assert total.num == 0;
  total = Fraction.add_all(fracs, min_complexity = 4, max_complexity = 7);
  assert total.num == 0;
# For add_all_new().
  Btest.<x,y> = PolynomialRing(ZZ);
  assert Fraction.add_all_new(
    [Fraction(Btest, 1, [x^2 - x, x^2 - y]),
    Fraction(Btest, -1, [x - y, x^2 - x]),
    Fraction(Btest, 1, [x - y, x^2 - y])]).evaluate() == 0
  fracs = [Fraction(Btest, x, [y])];
  for i in range(5):
    num = 1;
    den = x^i - y; 
    try:
      fracs += [Fraction(Btest, num, [den])];
      fracs += [Fraction(Btest, -num, [den])];
    except ZeroDivisionError:
      pass;
  shuffle(fracs);
  total = Fraction.add_all_new(fracs, min_complexity = 1, max_complexity = 1);
  total == Fraction(Btest, x, [y]);
  total = Fraction.add_all_new(fracs, min_complexity = 4, max_complexity = 7);
  total == Fraction(Btest, x, [y])


# For substitute()
  Btest.<x,y> = PolynomialRing(ZZ);
  frac = Fraction(Btest, x, [y]);
  target_ring = GF(7)
  vals = {x : 3, y : 5}
  new_frac = frac.substitute(vals, target_ring);
  assert new_frac.evaluate() == target_ring(3/5);
  frac = Fraction(Btest, 1, [x + 1, y - 1]);
  new_frac = frac.substitute(vals, target_ring);
  assert new_frac.evaluate() == target_ring(1/(4*4)); # fix longstanding bug
  
# For crossing_space()
  Ctest.<x,y,z> = PolynomialRing(ZZ);
  frac = Fraction(Ctest, 1, [x^2 - y, x + 3*z^2])
  assert frac.crossing_space() == (QQ^3).subspace([[2,-1,0],[1,0,-2]])
  frac = Fraction(Ctest, x + y + z, [])
  assert frac.crossing_space() == (QQ^3).subspace([])
  frac = Fraction(Ctest, x + y + z, [x])
  assert frac.crossing_space() == (QQ^3).subspace([])
  
# For hollow_decomp() and full_hollow_decomp
  frac = Fraction(Ctest, 1, [x, y]);
  assert frac.hollow_decomp() is None;
  frac = Fraction(Ctest, x + 2*y, [x, y]);
  assert frac.hollow_decomp() == [Fraction(Ctest, 1, [y]), Fraction(Ctest, 2, [x])]
  for frac in [
    Fraction(Ctest, x, [x + z, x - z]),
    Fraction(Ctest, x - 1, [x^5*y^6*z^5 - 1, x^3*y^3*z^4 - 1, x^2*y^2*z^3 - 1])
  ]:
    dec = frac.hollow_decomp();
    assert frac.evaluate() == sum(f.evaluate() for f in dec)
    dec = frac.full_hollow_decomp();
    assert frac.evaluate() == sum(f.evaluate() for f in dec)
    