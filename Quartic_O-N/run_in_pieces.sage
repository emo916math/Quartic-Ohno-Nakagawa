class Piece:
  # RingFactorData zone;
  # boolean symm;
  def __init__(self, zone, symm=False):
    self.zone = zone;
    self.symm = symm;
    
  def __str__(self):
    return "Piece(" + str(self.zone) + ", symm=" + str(self.symm) + ")"
    
  def __repr__(self):
    return str(self);
    
  def __add__(self,other):
    return Piece(self.zone + other.zone, self.symm and other.symm)

def run_piece(piece, hashing="none", quick=True, verbosity=0):
  # Allow naming pieces by name.
  if isinstance(piece, str):
    for p in pieces:
      if p.zone.name == piece:
        run_piece(p, hashing=hashing, quick=quick, verbosity=verbosity) # recursion
        return;
  
  summary = "Summary:\n";
  def print_important(*stuff):
    print(*stuff);
    return " ".join(str(item) for item in stuff) + "\n";

  print("Running piece", piece.zone.name);
  errors = 0;
  for Ftype in Ftypes_mod_symm if piece.symm else Ftypes:
    print("Function type", Ftype)
    t0 = time.time();
    ans = run(
      zone = main_zone,
      probe = piece.zone,
      Ftype_list = [Ftype],
      hashing = hashing,
      dual = True,
      add = True,
      asymmetric = not(piece.symm),
      quick = quick,
      verbosity=verbosity);
    discrep = ans[0].get(Ftype);
    if quick:
      if discrep is not None:
        cplx = discrep.complexity;
        summary += print_important("Complexity for ", Ftype, ":", cplx)
    else:
      if discrep is None or discrep == 0 or discrep == Fraction() or discrep == QuickFraction:
        summary += print_important("========== Piece passed in", time.time() - t0, "sec =========")
        
      else:
        summary += print_important("########## Piece failed in", time.time() - t0, "sec #########");
        summary += print_important("Discrepancy:")
        summary += print_important(discrep)
        errors += 1;
  if not(quick):
    summary += print_important(errors, "errors");
  
  print(summary);
  return errors;
  
def prove_provables(verbosity = 3):
  t0 = time.time(); errors = 0;
  for piece in pieces_easy:
    errors += run_piece(piece, quick=False, hashing="none", verbosity=verbosity)
  print("Total for provables:", errors, "errors")
  print("Total time for provables:", time.time() - t0)

def verify(verbosity = 3):
  t0 = time.time(); errors = 0;
  for piece in pieces_hard:
    errors += run_piece(piece, quick=False, hashing="full", verbosity=verbosity)
  print("Total for verifications:", errors, "errors")
  print("Total time for verifications:", time.time() - t0)