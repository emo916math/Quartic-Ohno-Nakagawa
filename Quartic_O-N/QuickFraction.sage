from collections import Counter
import time
attach("Fraction.sage")

test = True;
p=next_prime(10^9);

class QuickFraction(Fraction):
  # Instance variables:
  # base_ring -- a ring to work over
  # denoms
  # value -- a list of fractions
  # complexity -- the maximum complexity of a simplification that went
  # into this QuickFraction
  
  
  def __init__(self, *args):
    if len(args) == 0:
      # QuickFraction representing an empty sum.
      self.base_ring = ZZ; # coercible into any ring
      self.denoms = Counter();
      self.value = [];
    elif len(args) == 1:
      # Make a QuickFraction out of a Fraction
      fraction, = args;
      fraction.simplify();
      self.base_ring = fraction.base_ring;
      self.denoms = fraction.denoms;
      self.value = [fraction];
    elif len(args) == 3:
      # Make a QuickFraction out of base ring, denominators, and value
      base_ring, denoms, value = args;
      self.base_ring = base_ring;
      assert isinstance(denoms, Counter);
      for den in denoms:
        assert den in base_ring;
      self.denoms = denoms;
      for frac in value:
        # assert isinstance(frac, Fraction);
        assert frac.base_ring is base_ring;
      self.value = value;
    self.complexity = len(self.denoms);
    
  def __str__(self):
    ret = "QuickFraction(" + str(len(self.value)) + " fractions)"
    ret += " / (" + str(Factorization(self.denoms.items(),
      sort=False, simplify=False)) + "). Complexity: " + str(self.complexity);
    return ret;
    
  def __repr__(self):
    return str(self);
    
  def decomplicate(self, factorize=False):
    for den in denoms:
      assert den in base_ring; # do not attempt to clear them.
    
  def simplify(self, candidates=None, factorize=False, verb="auto"):
    if verb == "auto":
      try:
        verb = len(self.value) > 200 or len(self.denoms) > 10;
      except AttributeError:
        verb = False;
    if verb:
      print("Simplifying", self)
    dens_new = self.denoms;
    
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
          if den.hamming_weight() == 2:
            unknown = True;
            while unknown:
              try:
                subs = Fraction.divisibility_tester(self.base_ring, den, p);
              except UnsuitableDenominatorError:
                divisible = False;
                unknown = False;
                break;
              def test_sum(frac):
                if frac.denoms.get(den) is None:
                  return 0;
                frac_copy = Fraction(frac.base_ring, frac.num, copy(frac.denoms));
                if frac_copy.denoms.pop(den) != 1:
                  raise UnsuitableDenominatorError(); 
                return frac_copy.substitute(subs, GF(p)).evaluate();
              try:
                divisible = sum(test_sum(frac) for frac in self.value) == 0;
              except ZeroDivisionError:
                # Get a new divisibility tester
                continue;
              except UnsuitableDenominatorError:
                divisible = False;
              unknown = False; 
          else:
            if verb and den.hamming_weight() >= 3:
              print("<multinomial>", end=' ');
              sys.stdout.flush();
            divisible = False; # we'll do better
        except(ArithmeticError, AttributeError):
          raise;
          # divisible = num_new / den in self.base_ring;
        if divisible:
          if verb:
            print("//");
          dens_new -= Counter([den]);
        else:
          if verb:
            print("/\\");
          break;
    
    # Modify the fraction in place.
    self.denoms = dens_new;
    assert all(den.parent() is self.base_ring for den in self.denoms)
    
  def add(self, other, verb=True):
    new_ring = coercion_model.common_parent(self.base_ring, other.base_ring)
    dens1ctr = self.denoms;
    dens2ctr = other.denoms;
    dens_new = dens1ctr | dens2ctr;
    newfrac = QuickFraction(
      new_ring, dens_new, self.value + other.value
    );
    newfrac.complexity = max(newfrac.complexity,
      self.complexity, other.complexity)
    newfrac.simplify(dens1ctr & dens2ctr, verb = "auto" if verb else False);
    return newfrac;
  
  def evaluate(self):
    return self.value;
  
  def __eq__(self, other):
    try:
      return other == self.value;
    except:
      try:
        return other.value == self.value;
      except:
        return super().__eq__(other)

  def is_hollow(self):
    raise NotImplementedError("Not easy to tell whether a QuickFraction is hollow")

if test:
# For add_all_new().
  fracs = [];
  Btest.<x,y> = PolynomialRing(ZZ);
  assert len(Fraction.add_all(
    [QuickFraction(Fraction(Btest, 1, [x^2 - x, x^2 - y])),
    QuickFraction(Fraction(Btest, -1, [x - y, x^2 - x])),
    QuickFraction(Fraction(Btest, 1, [x - y, x^2 - y]))]).denoms) == 0
  for i in range(5):
    num = 1;
    den = x^i - y; 
    try:
      fracs += [QuickFraction(Fraction(Btest, num, [den]))];
      fracs += [QuickFraction(Fraction(Btest, -num, [den]))];
    except ZeroDivisionError:
      pass;
  shuffle(fracs);
  total = Fraction.add_all_new(fracs, min_complexity = 1, max_complexity = 1);
  assert len(total.denoms) == 0
  
  assert Fraction() != QuickFraction()
  assert QuickFraction() != Fraction()
  
class UnsuitableDenominatorError(Exception):
  pass