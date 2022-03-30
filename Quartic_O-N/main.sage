# Commands to prove the case of unramified resolvent (including the tame case)
reset()
load("quartic_ON_ur.sage");
verify(verbosity=3);
prove_provables(verbosity=2);

# Commands to prove the case of totally ramified resolvent
reset()
load("quartic_ON_1e3.sage");
verify(verbosity=2);
prove_provables(verbosity=2);

# Commands to prove the case of partially ramified resolvent
reset()
load("quartic_ON_1e21.sage");
verify(verbosity=2);
prove_provables(verbosity=2);