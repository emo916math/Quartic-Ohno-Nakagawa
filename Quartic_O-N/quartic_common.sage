import time
from collections import Counter

old_signs = False;

attach("LatticePointGeneratingFunctions.sage")
attach("quartic_examine.sage")
attach("QuickFraction.sage")
attach("run_in_pieces.sage")

# I. Hashing

p = next_prime(100000);
Fp = GF(p);

# V. Running and displaying.

PartialHash.<L0h,L2h> = PolynomialRing(Fp);

def display(tot_dic, hashing, add):
  if not(add):
    for (key, val) in tot_dic.items():
      print(key, ":", len(val), "terms")
  elif hashing in ZZ:
    print();
    for (key, val) in tot_dic.items():
      n_value = PartialHash(val);
      print(n_value);
      n = Matrix([[1 if n_value.coefficient(L0h^l0 * L2h^l2)
         .constant_coefficient() else 0 for l0 in [0..hashing]]
         for l2 in [0..hashing]])
      if n:
        print(n.str(zero=''), key);
  elif hashing == "full":
    print(tot_dic, end=' ');
  elif hashing == "none":
    for (key, val) in tot_dic.items():
      print(key, ":", end=' ');
      length = len(str(val))
      if length > 10000:
        print("<" + str(length) + " characters>")
      else:
        try:
          fzn = val.factor();
          fzn.sort(key = lambda x:x);
          print(fzn);
        except:
          print(val);
  elif hashing == "thought":
    for (key, val) in tot_dic.items():
      # Change to SHORT only _after_ all fraction operations.
      print(key, ":", shorten(val));
  else:
    raise ValueError(hashing);
  sys.stdout.flush()

def elim_l0l2(ring_factor_data):
  try:
    dic = {g0 : ring_factor_data.factor.get(G0),
      g1 : ring_factor_data.factor.get(G1)}
  except NameError: # different variables in new version of 1^2 1.
    dic = {l0 : ring_factor_data.factor.get(L0),
      l2 : ring_factor_data.factor.get(L2)}
  return ring_factor_data + RingFactorData(temps = dic);
# To be followed by elim_temps.

def has_no_l0l2(poly):
  return poly.degree(l0) == 0 and poly.degree(l2) == 0;
def ignore_l0l2(ring_factor_data):
  new_name = ring_factor_data.name;
  new_ieqs = [ieq for ieq in ring_factor_data.ieqs if has_no_l0l2(ieq)];
  new_eqns = [eqn for eqn in ring_factor_data.eqns if has_no_l0l2(eqn)];
  new_variables = [var for var in
    ring_factor_data.variables if has_no_l0l2(var)];
  new_factor = ring_factor_data.factor;
  new_Ftypes = ring_factor_data.Ftypes;
  return RingFactorData(new_name, new_ieqs, new_eqns, new_variables,
    new_factor, new_Ftypes);
  
# Parameter add can have 2 values:
# False: leave it as a large sum of fractions
# True: add everything as we go

# 0-5 verbosity levels:
# (0) Print nothing
# (1) Print salient input and output values
# (2) Print totals of each nonempty subzone
# (3) Print zones traversed, incl. empty ones
# (4) Print fraction addition
# (5) Dump info on everything

def run(
    zone,
    probe=RingFactorData(), hashing = "none", add=False, quick=False,
    primal=True, dual=False, asymmetric=False,
    zone_string = "", zone_list=None, excluded_zones=[], Ftype_list=None,
    top_level = True,
    verbosity = 1):
  
  if top_level and not(asymmetric):
    t0 = time.time();
    if verbosity >= 1:
      if probe.ieqs or probe.eqns or probe.temps:
        print("Probing inequalities", probe.ieqs,
          "equations", probe.eqns,
          "temps", probe.temps);
      if zone_string != "":
        print("Running zones that contain", repr(zone_string));
      if Ftype_list:
        print("Running types", Ftype_list);
        
  if add:
    if quick:
      zero = QuickFraction();
    else:
      zero = Fraction();
  else:
    zero = [];
  
  if asymmetric:
    if not(primal and dual):
      raise ValueError();
    Ftype = Ftype_list[0] if Ftype_list and len(Ftype_list) == 1 else None;
    probe2 = RingFactorData(
      ieqs=probe_dual(probe.ieqs, Ftype),
      eqns=probe_dual(probe.eqns, Ftype),
      temps=probe_dual(probe.temps, Ftype));
    Ftype_list2 = [dual_formulas.get(ty)[0]
      for ty in Ftype_list] if Ftype_list else None;
      
    run1, list1 = run(zone=zone, probe=probe, hashing=hashing, add=add,
      quick=quick, primal=True, dual=False,
      zone_string=zone_string, zone_list=zone_list, Ftype_list=Ftype_list,
      verbosity=verbosity);
    if verbosity > 0:
      print("============================");
    run2, list2 = run(zone=zone, probe=probe2, hashing=hashing, add=add,
      quick=quick, primal=False, dual=True,
      zone_string=zone_string, zone_list=zone_list, Ftype_list=Ftype_list2,
      verbosity=verbosity);
    
    if verbosity > 0:
      print("Discrepancy:", end=' ');
    
    assert not(old_signs);
    discrep = dict_add(run1, run2, zero=zero);
    if hashing == "partial":
      hashing = ZZ(e - probe.eqns[0]);
    if verbosity > 0:
      display(discrep, hashing, add);
      print()
    return discrep, list1, list2;
   
  if verbosity >= 3:
    print("Zone", zone.name);
  
  subzone = zone + probe;
  if False:
    print("Zone is now", subzone)
    print("Without temps:", subzone.elim_temps())
  if Ftype_list and subzone.Ftypes:
    # Check that it includes at least one relevant type.
    Ftype_relevant = False;
    if primal:
      if len(set(Ftype_list).intersection(subzone.Ftypes)) > 0:
        Ftype_relevant = True;
    if dual:
      Ftype_list2 = [dual_formulas.get(ty)[0] for ty in Ftype_list];
      if len(set(Ftype_list2).intersection(subzone.Ftypes)) > 0:
        Ftype_relevant = True;
  else:
    Ftype_relevant = True;
  if not(Ftype_relevant) or\
      (zone_list and all(not(entry.startswith(zone.name)) for entry in zone_list))\
      or subzone.is_empty() or subzone.elim_temps().is_empty():
    # Empty or excluded zone.
    if verbosity >= 6:
      print("Empty or excluded zone.")
    if (primal and dual) or not(top_level):
      return {}, {}, [];
    else:
      return {}, [];
  
  total_hash_p = {}; total_hash_d = {}; list_of_nonempty_zones = [];
  
  if zone.has_subcases():
    # Run recursively each subcase.
    for new_zone in zone.successors():
      hash_p, hash_d, nonempty_subzones = run(
        new_zone, probe=probe, hashing=hashing, add=add, quick=quick,
        primal=primal, dual=dual,
        zone_string=zone_string, zone_list=zone_list, excluded_zones=excluded_zones,
        Ftype_list=Ftype_list, top_level=False, verbosity=verbosity);
      if verbosity >= 9/2:
        print("Adding up zone", zone.name);
      dict_add(total_hash_p, hash_p, zero=zero, in_place=True);
      dict_add(total_hash_d, hash_d, zero=zero, in_place=True);
      list_of_nonempty_zones += nonempty_subzones;
      
  else:
    subzone = zone + probe;
    if subzone.name in excluded_zones\
        or zone_string not in subzone.name\
        or (zone_list and subzone.name not in zone_list):
      return {}, {}, [];
    if verbosity >= 6:
      print();
      print("Zone is:", subzone);
    subzone = elim_l0l2(subzone);
    if verbosity >= 6:
      print();
      print("After elim_l0l2:", subzone);
    subzone = subzone.elim_temps();
    if verbosity >= 6:
      print();
      print("After elim_temps:", subzone);
    
    
    if subzone.is_trivially_empty():
      # Excluded zone.
      if verbosity >= 6:
        print("Trivially empty zone.")
      return {}, {}, [];
    
    lat = lat_pt_gen_func(
      ieqs = subzone.ieqs,
      eqns = subzone.eqns,
      variables = subzone.variables,
      redundancy_check="full-cddlib",
      add=False,
      verbose = verbosity >= 5, # pass to LattE
      verb_conversion = verbosity >= 5,
      suppress_warnings=False,
      timing = verbosity >= 4
    );
    if verbosity >= 5:
      print(str(lat)[:50]);
    if lat == 0:
      if verbosity >= 5:
        print("Empty zone.")
      return {}, {}, [];
    
    try:
      subzone.factor.normalize();
    except:
      raise RuntimeError("Cannot normalize zone " + str(subzone));
    if verbosity >= 5:
      print();
      print("Zone is", subzone);
    
    if isinstance(lat, Fraction) and lat.num == 0:
      # Zone turned out to be empty.
      return {}, {}, [];
    
    ed_subs_only = ExponentDictionary(subzone.factor.dic)
    
    # Get hashing.
    if hashing == "none":
      hash_vals = {}; hash_ring = RINGS_RING;
    elif hashing == "thought":
      hash_vals = thought; hash_ring = THOUGHT;
    elif hashing == "partial" or hashing in ZZ:
      hash_vals = partial_hash; hash_ring = PartialHash;
      if hashing == "partial":
        try:
          hashing = ZZ(e - probe.eqns[0]);
        except(IndexError, TypeError):
          raise ValueError("Partial hash invalid: no value given for e")
    elif hashing == "full":
      hash_vals = full_hash; hash_ring = Fp;
    else:
      raise ValueError(hashing);
    
    # Tests whether zero-division occurs when hashing the zone answer.
    def is_good(frac):
      dens = frac.denoms;
      for den in dens:
        den_rings = subs_by_expdict(
          den, subzone.variables, RRR.gens(), ed_subs_only) ;
        den_hashed = den_rings.substitute(hash_vals);
        if den_hashed == 0:
          return False;
      return True; # TODO: Take into account dual also, to avoid complications 
                   # from dependencies in the hash values.
    
    good = []; bad = [];
    if isinstance(lat, Fraction):
      good = [lat];
      RRR = lat.base_ring;
    else:
      RRR = lat[0].base_ring;
      for frac in lat:
        if is_good(frac):
          good += [frac];
        else:
          bad += [frac];
    while len(bad) > 0:
      if verbosity >= 3:
        print(len(bad), "bad fractions")
      if len(bad) == 1:
        print("Lone fraction", bad[0])
        raise ValueError("Infinite zone " + subzone.name + " for probe " +\
          str(probe) + " in hashing " + repr(hashing)) 
      f1 = bad.pop();
      f2 = bad.pop();
      frac = f1 + f2;
      if is_good(frac):
        good += [frac];
      else:
        bad += [frac];
    if verbosity == 2:
      print("Zone", subzone.name)
    if verbosity >= 3:
      print(len(good), "good fractions")
    
    ring_mult = subzone.factor.const;
    if len(subzone.variables) > 0:
      ring_subs, mu = expdict_to_subs_multiplier(subzone.variables, RRR.gens(),
        ed_subs_only);
      assert mu == 1;
    else:
      ring_subs = {};
    ring_tot = [frac.substitute(ring_subs, RINGS_RING)
        .scalar_mult(ring_mult) for frac in good];
    
    zone_added = False;
    if primal:
      hash_tot_p = [frac.substitute(hash_vals, hash_ring)
        for frac in ring_tot]
      for frac in hash_tot_p:
        frac.simplify(factorize=True);
      if quick:
        hash_tot_p = [QuickFraction(f) for f in hash_tot_p]
      if add:
        hash_tot_p = Fraction.add_all_new(hash_tot_p,
          verb = verbosity >= 4);
      if hash_tot_p != 0:
        list_of_nonempty_zones += [subzone.name];
        zone_added = True;
        for Ftype in subzone.Ftypes:
          if Ftype_list and Ftype not in Ftype_list:
            continue;
          new_entry = {Ftype : hash_tot_p}
          if verbosity >= 2:
            display(new_entry, hashing, add)
          total_hash_p = dict_add(total_hash_p, new_entry, zero=zero);
    if verbosity >= 2 and primal and dual:
      print ("| ", end = "");
    if dual:
      for Ftype in subzone.Ftypes:
        if not(dual_formulas.get(Ftype)):
          raise KeyError(Ftype)
        new_type, ans_duals = dual_formulas.get(Ftype);
        if Ftype_list and new_type not in Ftype_list:
          continue;
        ans_dual, _, _ = ans_duals; # shed thought-dual and probe-dual
        
        if old_signs:
          ring_tot_d = [frac.substitute(ans_dual, RINGS_RING).scalar_mult(1)
          for frac in ring_tot];
        else:
          ring_tot_d = [frac.substitute(ans_dual, RINGS_RING).scalar_mult(-1)
            for frac in ring_tot];
        hash_tot_d = [frac.substitute(hash_vals, hash_ring)
          for frac in ring_tot_d]
        for frac in hash_tot_d:
          frac.simplify(factorize=True);
        if quick:
          hash_tot_d = [QuickFraction(f) for f in hash_tot_d]
        if add:
          hash_tot_d = Fraction.add_all_new(hash_tot_d,
            verb = verbosity >= 4)
        if hash_tot_d == 0:
          continue;
        if not(zone_added):
          list_of_nonempty_zones += [subzone.name]; 
        new_entry_d = {new_type : hash_tot_d};
        if verbosity >= 2:
            display(new_entry_d, hashing, add)
        total_hash_d = dict_add(total_hash_d, new_entry_d, zero=zero);
    if verbosity >= 2:
      print();
      
  if primal and dual and top_level:
    if old_signs:
        discrep = dict_add(total_hash_p, total_hash_d,
          operation=(lambda x,y: x-y),
          zero=zero);
    else:
      discrep = dict_add(total_hash_p, total_hash_d,
          operation=(lambda x,y: x+y),
          zero=zero);
      
  if verbosity >= 1 and top_level:
    print(len(list_of_nonempty_zones), "nonempty zones found:");
    for zn in list_of_nonempty_zones:
      print(zn);
    if primal:
      print()
      print("Total:", end=' ');
      display(total_hash_p, hashing, add);
    if dual:
      print()
      print("Dual total:", end=' ');
      display(total_hash_d, hashing, add);
    if primal and dual and top_level:
      print()
      print("Discrepancy:", end=' ');
      display(discrep, hashing, add);
    print("");
    print("Total time:", time.time() - t0)
    sys.stdout.flush();
  if not(top_level):
    # Return a standardized format for adding.
    return total_hash_p, total_hash_d, list_of_nonempty_zones;
  if primal and dual:
    return discrep, total_hash_p, total_hash_d, list_of_nonempty_zones;
  if primal:
    return total_hash_p, list_of_nonempty_zones;
  if dual:
    return total_hash_d, list_of_nonempty_zones;
    # list_of_nonempty_zones is used for examination


# To run it for a value of the reduced discriminant only.
def quick_check(level, hashing = "full", **kwargs):
  return run(hashing, "", probe = RingFactorData(eqns = [red_disc - level]),
    dual = True, **kwargs)

def test_case(ieqs = [], eqns = [], expected = 0, verbosity = 1/2):
  run_ans, _ = run(zone = main_zone, hashing = "none", add = True,
    probe = RingFactorData(ieqs = ieqs, eqns = eqns),
    verbosity = verbosity);
  run_ans = {Ft : fr.evaluate() for (Ft, fr) in run_ans.items()};
  ans = deep_subs(run_ans, thought);
  if ans != expected:
    if verbosity > 0:
      print("Error for inequalities", ieqs, "equations", eqns, " Expected:")
      print(expected)
      print("Got:")
      print(ans);
      print("Discrepancy:")
      print(dict_add(ans, expected, operation=lambda x,y:x-y, zero=0));
    return 1;
  else:
    if verbosity > 0:
      print("--------- Passed ---------")
    return 0;

# "cases" is defined per spl.t.
def test_all_cases(case_list = None, verbosity = 1/2):
  errors = 0;
  for (name, eqns, expected) in cases:
    if case_list is None or name in case_list:
      print(name);
      errors += test_case(eqns=eqns, expected=expected, verbosity=verbosity)
  print(errors, "errors")
  
def test_random_case(verbosity = 1/2):
  name, eqns, expected = cases[randint(0,len(cases) - 1)]
  print("Randomly chose", name);
  test_case(eqns=eqns, expected=expected, verbosity=verbosity)
  
def probe_dual(info, Ftype=None):
  if Ftype:
    return deep_subs(info, dual_formulas.get(Ftype)[1][2]);
  if not deep_subs(info, {l0 : 0, l2 : 0}) == info:
    raise ValueError("Cannot dualize " + str(info) +
      " without knowledge of Ftype")
  ret = deep_subs(info, {t : e - t});
  if isinstance(ret, dict) and t in ret.keys():
    ret.update({t : ret.get(e) - ret.get(t)});
  return ret;

THOUGHT.<qTH,xTH,yTH,zTH> = PolynomialRing(ZZ, order='invlex')
# THOUGHTL.<qTH,xTH,yTH,zTH> = LaurentPolynomialRing(ZZ, order='invlex')
def thought_dual(tot_dict):
  ret = {};
  for (Ftype, total) in tot_dict.items():
    try:
      total = THOUGHT(total);
    except TypeError:
      total = THOUGHT.fraction_field()(total);
    new_total = deep_subs(total,
      dual_formulas.get(Ftype)[1][1]);
    ret.update({dual_formulas.get(Ftype)[0]: new_total});
  return ret;
  
q_ring = PolynomialRing(ZZ, names="q")
SHORT = PolynomialRing(q_ring, names="x,y,z", order='lex')
short = {qTH : q_ring("q"), 
         xTH : SHORT("x"),
         yTH : SHORT("y"),
         zTH : SHORT("z")}
def shorten(total):
  if isinstance(total, Fraction):
    total = total.evaluate()
  return deep_subs(total,short)