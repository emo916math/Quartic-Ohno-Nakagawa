# set_random_seed(0) # For consistency. Will be removed.

attach("quartic_common.sage")


zone_ring.<e, t, dpf, dpffl, o_d, b1f, b2f, spf, o_s, l_d, l_m, letter, h_eta, sbarfl, a1pf, o1, a2pf, o2, o3, a2pbf, o2bar, a2bar_type, o1_strong, o1_strong_parity, k, l, g0, g1, vK_xi1f, vK_xi1, vK_xi1_strong, vK_xi2f, k12, l_d_term, m11fl, k11ceil, n11fl, ncfl, ncflz, ntil, m22fl> = PolynomialRing(QQ);

RINGS_RING.<E_,T_,DPFL_,B1F_,B2F_,SPFL_,L_M_,q, G0,G1> = PolynomialRing(ZZ, order='lex');


# Basic indices
d0 = 2*dpf + o_d;
d0pr = d0/2;
spr = spf + o_s/2;
sbar = spr + d0pr;

a1p = a1pf + o1/4;
a2p = a2pf + o2/4;
a1pbar = a1p;
# a2pbar = a2pbar_times_4/4;
a2pbar = a2pbf + o2bar/4;

ob1 = o_s;
ob2 = 1 - ob1;
b1 = b1f + ob1/2;
b2 = b2f + ob2/2;
b1bar = b1;
b2bar = b2 + (d0 - 1)/2;
red_disc = b1f + b2f;

a1 = a1p - 2*b1;
a1f = a1 - o1/4;
a2 = a2p - 2*b1;
a2f = a2 - o2/4;
a1bar = a1;
a2bar = a2pbar - 2*b1;
a2bf = a2bar - o2bar/4;

a3 = (b1 + b2) + 2*sbar + 4*t - 4*e - a1 - a2;
a3p = a3 + 2*b1;
a3bar = (b1bar + b2bar) + 2*sbar + 4*t - 4*e - a1bar - a2bar;
a3pf = a3p - o3/4;

# vK_xi1 = vK_xi1f + o1/4;
vK_xi2 = vK_xi2f + o2/4;

# q = rq^2;

# For examination
os = [o_d, ob1, ob2, o_s, o1, o2, h_eta];
o_denom = 4; # For the two exam o's
flav002f = RingFactorData("002", temps = {o1 : 0, o2 : 0, o3 : 2, vK_xi1: vK_xi1f + o1/4}, factor = ({q : (-0)}))
flav020f = RingFactorData("020", temps = {o1 : 0, o2 : 2, o3 : 0, vK_xi1: vK_xi1f + o1/4}, factor = ({q : (-1)}))
flav200f = RingFactorData("200", temps = {o1 : 2, o2 : 0, o3 : 0, vK_xi1: vK_xi1f + o1/4}, factor = ({q : (-2)}))
flav123f = RingFactorData("123", temps = {o1 : 1, o2 : 2, o3 : 3, vK_xi1: vK_xi1f + o1/4}, factor = ({q : (-0)}))
flav132f = RingFactorData("132", temps = {o1 : 1, o2 : 3, o3 : 2, vK_xi1: vK_xi1f + o1/4}, factor = ({q : (-1)}))
flav231f = RingFactorData("231", temps = {o1 : 2, o2 : 3, o3 : 1, vK_xi1f: 0, vK_xi1: 0 }, factor = ({q : (-2)}))
flav213f = RingFactorData("213", temps = {o1 : 2, o2 : 1, o3 : 3, vK_xi1f: 0, vK_xi1: 0 }, factor = ({q : (-1)}))
flav312f = RingFactorData("312", temps = {o1 : 3, o2 : 1, o3 : 2, vK_xi1: vK_xi1f + o1/4}, factor = ({q : (-2)}))
flav321f = RingFactorData("321", temps = {o1 : 3, o2 : 2, o3 : 1, vK_xi1: vK_xi1f + o1/4}, factor = ({q : (-3)}))
a_flavors_0 = [
  flav002f, flav020f, flav200f,
  flav123f, flav132f, flav231f, flav213f, flav312f, flav321f
];
flavor_type = (o1 + o2 + o3 - 2)/4;

# For examination only
o2bar_0 = RingFactorData("o2bar=0", temps = {o2bar : 0})
o2bar_1 = RingFactorData("o2bar=1", temps = {o2bar : 1})
o2bar_2 = RingFactorData("o2bar=2", temps = {o2bar : 2})
o2bar_3 = RingFactorData("o2bar=3", temps = {o2bar : 3})
o2bar_flavors = [o2bar_0, o2bar_1, o2bar_2, o2bar_3]

a_flavors = [r + s for r in a_flavors_0 for s in o2bar_flavors]

x_floor = a1pf; x_o = o1; x_exam = a1p;
y_floor = a2pf; y_o = o2; y_exam = a2p;

z_floor = a2pbf; z_o = o2bar; z_exam = a2pbar;

body_ring.<x_floor_plot,y_floor_plot, z_floor_plot> = PolynomialRing(QQ)
slice_ring.<y_slice, z_slice> = PolynomialRing(QQ)

# Hashing
PartialHash.<G0h,G1h> = PolynomialRing(Fp);
partial_hash = {
  E_    : Fp.random_element(),
  T_    : Fp.random_element(),
  DPFL_ : Fp.random_element(),
  B1F_  : Fp.random_element(),
  B2F_  : Fp.random_element(),
  SPFL_ : Fp.random_element(),
  L_M_  : Fp.random_element(),
  q     : Fp.random_element(),
  G0 : G0h, G1 : G1h
};

full_hash = partial_hash.copy();
full_hash.update({
  G0 : Fp.random_element(),
  G1 : Fp.random_element(),
});

thought = {
  q : qTH,
  E_ : 1/zTH,
  T_ : 1,
  DPFL_ : 1,
  B1F_ : 1,
  B2F_ : 1,
  SPFL_ : 1,
  L_M_ : 1,
  G0 : xTH,
  G1 : yTH
}
short = {qTH : q_ring("q"), 
         xTH : SHORT("x"),
         yTH : SHORT("y"),
         zTH : SHORT(1)} # Hide z in short answers.

# Overrides the corresponding function in quartic_examine.sage.
def pgn_and_pts(row = None):
  assert x_exam == x_floor + x_o / o_denom;
  assert y_exam == y_floor + y_o / o_denom;
  assert z_exam == z_floor + z_o / o_denom;
  
  if row is None:
    plot_ring = body_ring;
    plot_subs = dict_add(resolvent_data,
      {x_floor : x_floor_plot, y_floor : y_floor_plot, z_floor : z_floor_plot,
      x_o : 0, y_o : 0, z_o : 0})
  else:
    plot_ring = slice_ring;
    plot_subs = dict_add(resolvent_data,
      {x_floor : row, y_floor : y_slice, z_floor : z_slice,
      x_o : 0, y_o : 0, z_o : 0})
  
  if False:
    print(plot_subs);
    for ieq in main_zone.ieqs:
      print("Step 0:", ieq)
      print("Step 1/2:", [term.subs(plot_subs) for term in ieq.monomials()])
      print("Step 1:", ieq.subs(plot_subs))
      print("Step 2:", plot_ring(ieq.subs(plot_subs)));
      print("Step 3:", linear_to_coeffs(plot_ring(ieq.subs(plot_subs))));
  plot_ieqs = [linear_to_coeffs(plot_ring(ieq.subs(plot_subs)))
    for ieq in main_zone.ieqs]
  plot_eqns = [linear_to_coeffs(plot_ring(eqn.subs(plot_subs)))
    for eqn in main_zone.eqns]
  pgn_0 = Polyhedron(ieqs = plot_ieqs, eqns = plot_eqns);
  
  pts = [];
  for a_flav in a_flavors:
    if row is None:
      vec = (QQ^3)((a_flav.temps.get(x_o)/o_denom,
                    a_flav.temps.get(y_o)/o_denom,
                    a_flav.temps.get(z_o)/o_denom));
    else:
      if a_flav.temps.get(x_o) / o_denom - row not in ZZ:
        # print("Skipping flavor", a_flav.name);
        continue;
      vec = (QQ^2)((a_flav.temps.get(y_o)/o_denom,
                    a_flav.temps.get(z_o)/o_denom));
      # print("         Flavor", a_flav.name, "vec =", vec)
    pgn = pgn_0.translation(-vec)
    pts += [pt + vec for pt in pgn.integral_points()
      if
      plot_ring(m12bar.subs(plot_subs))(*(pt + vec)) <= 0 or
      plot_ring(m22bar.subs(plot_subs))(*(pt + vec)) <= 0
    ]
  pts.sort();
  return pgn_0, pts;

a = [None, a1, a2, a3];
abar = [None, a1bar, a2bar, a3bar];

def mijbar(i,j):
  return b2bar - abar[i] - abar[j] - (2 if i == j else 3) * (e - t) + sbar;
def mij(i,j):
  return b2bar - a[i] - a[j] - (2 if i == j else 3) * (e - t) + sbar;
def nijbar(i,j):
  return b1 - abar[i] - abar[j] - (2 if i == j else 3) * (e - t) + 2*sbar;

m11 = mij(1,1); m12 = mij(1,2); m13 = mij(1,3); m22 = mij(2,2)
m12bar = mijbar(1,2);
m13bar = mijbar(1,3);
m22bar = mijbar(2,2);
m33bar = mijbar(3,3);
n11 = nijbar(1,1); n12bar = nijbar(1,2); n22bar = nijbar(2,2);

n_c = n11 - sbar; # the colorful normalization
m_c = m11 - d0pr;

l0 = (g0 - g1 - d0pr - h_eta)/2;
l2 = (2*e - g0 - g1 - d0pr - h_eta)/2;
assert g0 == e + l0 - l2;
assert g1 == e - d0pr - l0 - l2 - h_eta
# l1 = e - d0pr - l0 - l2;

# Possibilities for d0
d0_even = RingFactorData("d0_even",
  ieqs = [dpf, e - dpf],
  temps = {o_d : 0}
)
d0_odd = RingFactorData("d0_odd",
  temps = {o_d : 1, dpf : e}
);
d0_options = SubcaseSwitch([d0_even, d0_odd], priority = 1)

d0pr_even = RingFactorData("d0'_even",
  temps = {dpf : 2*dpffl}, variables = [dpffl])
d0pr_odd = RingFactorData("d0'_odd",
  temps = {dpf : 2*dpffl + 1}, variables = [dpffl])
d0_mod_4_options = SubcaseSwitch([d0pr_even, d0pr_odd], priority = -2)
dpfceil = dpf - dpffl;
dpfparity = dpfceil - dpffl;

# For the gray zone
# m11fl = floor(m11/2)
m11_even = RingFactorData("m11_even", eqns = [m11 - 2*m11fl],
  variables = [m11fl])
m11_odd = RingFactorData("m11_odd", eqns = [m11 - 2*m11fl - 1],
  variables = [m11fl])
m11_options = SubcaseSwitch([m11_even, m11_odd], priority = -2);
m11ceil = m11 - m11fl;
m11_h_fl = m11ceil*h_eta + m11fl*(1 - h_eta) # floor((m_c + h_eta)/2)

# For the ivory zone: k11 is near m11/2
k11 = m11/2 - o1/4;
k11_options = [
  RingFactorData(
    "k11_"+str(i),
    eqns = [k11ceil - k11 - i/4],
    variables = [k11ceil]
    )
  for i in [0,1,2,3]
]
k11_options_even = SubcaseSwitch([k11_options[0], k11_options[2]], priority = -2);
k11_options_odd  = SubcaseSwitch([k11_options[1], k11_options[3]], priority = -2);

m22_adj = m22 - o2/2;
m22_even = RingFactorData("m22_even", eqns = [m22_adj - 2*m22fl],
  variables = [m22fl]);
m22_odd = RingFactorData("m22_odd", eqns = [m22_adj - 2*m22fl - 1],
  variables = [m22fl]);
m22_options = SubcaseSwitch([m22_even, m22_odd], priority = -2)
m22ceil = m22_adj - m22fl;

# n11fl = floor(n11/2)
n11_even = RingFactorData("n11_even", eqns = [n11 - 2*n11fl],
  variables = [n11fl])
n11_odd = RingFactorData("n11_odd", eqns = [n11 - 2*n11fl - 1],
  variables = [n11fl])
n11_options = SubcaseSwitch([n11_even, n11_odd], priority = -2);
n11ceil = n11 - n11fl;
n11_parity = n11 - 2*n11fl;
ncflz_CD = n11fl - (sbar + h_eta)/2; # floor((n_c - h_eta)/2)
ncflz_E = n11ceil - (sbar + 1)/2;

# ncfl = floor(n_c/2)
nc_even = RingFactorData("nc_even", eqns = [n_c - 2*ncfl], variables = [ncfl])
nc_odd = RingFactorData("nc_odd", eqns = [n_c - 2*ncfl - 1], variables = [ncfl])
nc_options = SubcaseSwitch([nc_even, nc_odd], priority = -2);
ncceil = n_c - ncfl;
nc_parity = n_c - 2*ncfl;

# n11fl_brown = floor((n11 - d0 - o1/2)/2).
n11_brown = n11 - d0 - o1/2 + o_d/2;

n11_even_brown = RingFactorData("n11_brown_even", eqns = [n11_brown - 2*n11fl],
  variables = [n11fl]);
n11_odd_brown = RingFactorData("n11_brown_odd", eqns = [n11_brown - 2*n11fl - 1],
  variables = [n11fl]);
n11_brown_options = SubcaseSwitch([n11_even_brown, n11_odd_brown], priority = -2);
n11ceil_brown = n11_brown - n11fl;

# In types C and D, ntil = \floor{\frac{2e - s' - d_0' - n_c + 2 - 2h_\eta}{4}}.
# Used with n11_options.
ntil_near = RingFactorData("ntil_near", eqns = [e - n11ceil - h_eta - 2*ntil],     variables = [ntil])
ntil_far  = RingFactorData("ntil_far",  eqns = [e - n11ceil - h_eta - 2*ntil + 1], variables = [ntil]);
ntil_options = SubcaseSwitch([ntil_near, ntil_far], priority = -2)
ntil_minus = e - n11ceil - h_eta - ntil;

# Letter types
type_A = RingFactorData("A",
    temps = {letter : 1, o_s : 1, spf : -1, l_d : -1/2, vK_xi1_strong : (d0 - 1)/4,
       o1_strong_parity : 1 - o_d},
    subcases = [d0_options],
    )
type_B = RingFactorData("B",
    ieqs = [spf, (d0/2 - 1 - 1/2) - spf],
    temps = {letter : 2, o_s : 0, vK_xi1_strong : (2*spf + d0)/4, l_d : spf,
       o1_strong_parity : o_d},
    subcases = [d0_options],
    )
type_C = RingFactorData("C",
    ieqs = [e + 1 - l_d, l_m],
    temps = {letter : 3, o_s : 0, spf : dpf - 1, h_eta : 1, l_d : dpf - 1 + 2*l_m,
      ncflz : ncflz_CD, sbarfl : dpf - 1,
      vK_xi1_strong : (d0 - 1)/2 + (2 - o1)/4, # big enough
      o_d : 0, o1_strong : 2, o1_strong_parity : 0},
    variables = [l_m],
    factor = ({L_M_: l_m})
    )
type_D_d0_even = \
  RingFactorData("D_d0_even",
    ieqs = [spr - dpf, e + 1 - l_d, l_m],
    temps = {letter : 4, o_s : 0, spf : 2*sbarfl - dpf, h_eta : 0, l_d : dpf + 2*l_m,
      ncflz : ncflz_CD, vK_xi1_strong : d0/2, # big enough
      o_d : 0, o1_strong : 0, o1_strong_parity : 0},
    variables = [sbarfl, l_m],
    factor = ({L_M_: l_m})
    )
type_D_d0_odd = \
  RingFactorData("D_d0_odd",
    ieqs = [spr - e],
    temps = {letter : 4, o_s : 0, h_eta : 0, l_d : e, vK_xi1_strong : e + o1/4,
      o_d : 1, o1_strong_parity : 1},
    variables = [h_eta],
    ) + d0_odd
type_E = RingFactorData("E",
    ieqs = [spr - dpf - 1],
    temps = {letter : 5, o_s : 0, spf : 2*sbarfl - dpf + 1, h_eta : 0, l_d : dpf - 1,
      o1_strong : 2, ncflz : ncflz_E,
      vK_xi1_strong : (d0 - 1)/2 + (2 - o1)/4, # big enough
      o_d : 0, o1_strong : 2, o1_strong_parity : 0},
    variables = [sbarfl, h_eta, l_d],
    )

sbarceil = sbar - sbarfl;

letter_types = SubcaseSwitch(
  [type_A, type_B, type_C, type_D_d0_even, type_D_d0_odd, type_E], priority = 1)

# %F is either F, %balF, or %sideF: see p. XIII.31
E_options = SubcaseSwitch([
  RingFactorData("F", ieqs = [l_d - (g0 - g1), e - g0], Ftypes = ["F"]),
  RingFactorData("Ebal", ieqs = [l_d - (2*e - g0 - g1), g0 - e - 1], Ftypes = ["Ebal"]),
  RingFactorData("Eside", ieqs = [(g0 - g1) - l_d - 1/2, (2*e - g0 - g1) - l_d - 1/2,
      spr + g0, # xF_A
      g0 + (4 - letter)*(g1 + d0pr + 1), # xF_E
    ],
    Ftypes = ["Eside"]),
  RingFactorData("xF_A", eqns = [g0 - 0, g1 - 0, spr + 1/2], Ftypes = ["xF"]),
  RingFactorData("xF_E", eqns = [l0, letter - 5], Ftypes = ["xF"])
], priority = -1)

distinctnesses = [
  # If all i's are distinct. 
  RingFactorData(
    "distinct",
    ieqs = [4*(a[2] - a[1]) - 1, 4*(a[3] - a[2]) - 1],
    factor = ExponentDictionary({q : (2*a3pf - 2*a1pf)})
    ),
  
  # If a1 == a2 < a3.
  RingFactorData(
    "a1==a2",
    ieqs = [4*(a[3] - a[2]) - 1],
    eqns = [a[2] - a[1]],
    factor = ExponentDictionary({q : (2*a3pf - 2*a1pf)}, 1/(1 + 1/q))
    ),
  # If a1 < a2 == a3.
  RingFactorData(
    "a2==a3",
    ieqs = [4*(a[2] - a[1]) - 1],
    eqns = [a[3] - a[2]],
    factor = ExponentDictionary({q : (2*a3pf - 2*a1pf)}, 1/(1 + 1/q))
    ),
]

a2bar_min = a2bar - a2;
a2bar_d0 = (d0 - 1)/2 - (a2bar - a2);
a2bar_a3 = a3 - a2bar;
a2bar_xi1 = vK_xi1 - (a2bar - a2);
a2bar_xi2 = vK_xi2 - (a2bar - a2);

# k12 = ceil(a2bar - a2 - o2/4). There are simplified formulas for k12 in
# some a2bar types.
a2bar__min = RingFactorData(
    "a2bar:min",
    ieqs = [o1 - o2 + 1, o2 - o1 + 1],
    # easy flavors: 002, 123, 213, 231, 321 where one of xi1, xi2 is always K-led.
    eqns = [a2bar - a2],
    temps = {a2bar_type : 0, k12 : 0, o2bar : o2},
    )
a2bar__d0 = RingFactorData(
    "a2bar:d0",
    eqns = [
      a2bar_d0,
      o1 + o2 - o3 - 2, # flavors 020, 200, 132, 312 only
    ],
    temps = {a2bar_type : 3},
    subcases = [SubcaseSwitch(
      [
        RingFactorData("#0#", temps = {o2 : 0, o2bar : 2 - 2*o_d, k12 : dpf}),
        RingFactorData("#1#", temps = {o2 : 1, o2bar : 3 - 2*o_d, k12 : dpf}),
        RingFactorData("#2#", temps = {o2 : 2, o2bar : 0 + 2*o_d, k12 : dpf - 1 + o_d}),
        RingFactorData("#3#", temps = {o2 : 3, o2bar : 1 + 2*o_d, k12 : dpf - 1 + o_d})
      ]
    , priority = 0)]
    )
a2bar__a3 = RingFactorData(
    "a2bar:a3",
    eqns = [
      a2bar_a3,
      o1 + o2 - o3 - 2, # flavors 020, 200, 132, 312 only
    ],
    ieqs = [a2bar_d0 - 1/4],
    temps = {a2bar_type : 4, k12 : a3 - a2 - o2/4, # See XII.346
      o2bar : o3},
    )
a2bar__xi1 = RingFactorData(
    "a2bar:xi1",
    eqns = [a2bar_xi1,
      o1 + o2 - o3 - 2, # flavors 020, 200, 132, 312 only
    ],
    ieqs = [a2bar_d0 - 1/4, a2bar_a3 - 1/4],
    temps = {
      a2bar_type : 2, 
      vK_xi1f : a2bar - a2 - o1/4,
      k12 : vK_xi1f + (o1 - o2 + 2)/4, # see XIII.24
      o2bar : 2 - o3,
    },
    )
a2bar__xi2 = RingFactorData(
    "a2bar:xi2",
    eqns = [a2bar_xi2,
      o1 + o2 - o3 - 2, # flavors 020, 200, 132, 312 only
    ],
    ieqs = [
      a2bar_d0 - 1/4, a2bar_a3 - 1/4
    ],
    temps = {
      a2bar_type : 1, 
      vK_xi2f : a2bar - a2 - o2/4,
      k12 : vK_xi2f,
      o2bar : o3 },
      factor = 1 - 1/q, # in all xi2_options? TODO
    )

a2bar_options = [a2bar__min, a2bar__d0, a2bar__a3, a2bar__xi1, a2bar__xi2]

# o1_repl = 2*sbar - 4*sbarfl;
o1_repl = o1; # for compatibility
XX.<xx> = PolynomialRing(QQ);
mod2 = (-1/3) * (xx - 4) * xx * (xx - 2)^2 # returns right answer on inputs 0,1,2,3,4

normalize_g1 = SubcaseSwitch([
  RingFactorData(">",
    ieqs = [g1]
    ),
  RingFactorData("^",
    eqns = [g1 + 1],
    factor = G1, # keeps the written g1-value the same. See XIII.69.
    )
])
black = RingFactorData(
  "black",
  ieqs = [n_c - 2*e - 1/2], # plum
  temps = {
    l_d_term : l_d
  },
  factor = ({G0 : 2*e, G1 : zone_ring(0),
     q : (2*e - m11 - n11 + sbar/2 + d0/2 + o1_repl/4)
  }),
  Ftypes = ["F", "Fx"],
  )
plum = RingFactorData(
  "plum",
  ieqs = [
    2*e - n_c, # black
    n_c - (2*e - 2*d0pr + 2), # purple
    n_c - (2*e - 2*spr), # blue, green
  ],
  eqns = [o_s], # not letter A
  temps = {
    l_d_term : l_d
  },
  factor = ({G0 : e + ncfl, G1 : 0,
    q : (e - m11 + (d0 + o1_repl - 2*spr)/4 - ncceil)
  }),
  Ftypes = ["F"],
  subcases = [nc_options],
  )
purple = RingFactorData(
  "purple",
  ieqs = [
    2*e - n_c, # black
    (2*e - 2*d0pr + 1) - n_c, # plum
    n_c - (2*e - 2*sbarceil + 1 + h_eta), # blue/green
    m11 - 2*e - 1, # gray
    letter - 4, # types D(even), E only
  ],
  temps = {
    l_d_term : l_d,
    o_d : 0
  },
  factor = ({G0 : e + ncfl, G1: e - d0pr - ncfl,
      # L0 : ncfl, L2 : zone_ring(0),
    q : (e - m11 + (d0 + o1_repl - 2*spr)/4 - ncceil)}),
  Ftypes = ["F", "Fx"],
  subcases = [nc_options]
  # If problems with fractional exponents, try using sbarfl
  )
blue = RingFactorData(
  "blue",
  ieqs = [
    (2*e - d0pr - spr) - n_c, # purple
    n_c - (2*e + d0pr - spr - 2*l_d - 1), # green
    n_c - (2*e - spr - d0pr)/3, # red
    m_c - (e + (n_c + spr - d0pr + 1)/2) - 1/2, # gray
    letter - 3, 4 - letter, # types C, D(even) only
  ],
  temps = {
    l_d_term : l_d,
    o_d : 0
  },
  factor = ({G0 : e + ncflz_CD - ntil, G1 : e - d0pr - ncflz_CD - ntil - h_eta, 
      # L0: ncflz_CD, L2 : ntil,
    q : (-m11 + d0pr + (2*h_eta + o1)/4 + ntil_minus)}),
  Ftypes = ["F"],
  subcases = [n11_options, ntil_options, normalize_g1]
  )
greenAB = RingFactorData(
  "greenAB",
  ieqs = [
    (2*e - 2*spr - 1) - n_c, # plum
    2 - letter, # types A and B only
  ],
  temps = {
    l_d_term : spr
  },
  factor = ({q : (-m11 + (d0 + o1_repl + 2*spr)/4),
    G0 : n_c + spr, G1 : 0}),
  subcases = [
    [
      RingFactorData("F", Ftypes = ["F"]),
      RingFactorData("eF",
        ieqs = [m11 - (2*e + d0pr - spr)], # m11-smear: see p. XII.345
        Ftypes = ["eF"]),
  ]]
  )
green = RingFactorData(
  "green",
  ieqs = [
    (2*e - 2*sbarceil + h_eta) - n_c, # purple
    (2*e + d0pr - spr - 2*l_d - 2) - n_c, # blue
    n_c - (l_d - d0pr + 2), # red
    m_c - (n_c + spr + l_d - d0pr + 1) - 1/2, # gray
    letter - 3, # types C,D(even),E only
  ],
  temps = {
    l_d_term : l_d,
    o_d : 0,
  },
  factor = ({q : (-m11 + (d0 + o1_repl + 2*l_d)/4),
       G0 : n_c + spr/2 + l_d/2,
       G1 : - d0pr - 2*ncflz + n_c + spr/2 + l_d/2 - h_eta
    # L0 : ncflz, L2 : e - n_c + ncflz - spr/2 - l_d/2
  }),
  subcases = [
    n11_options, d0_mod_4_options,
    [
      RingFactorData("F", Ftypes = ["F"]),
      RingFactorData("eF",
        ieqs = [m11 - (2*e + d0pr - l_d)], # m11-smear
        Ftypes = ["eF"]),
    ],
    normalize_g1]
  )
red = RingFactorData(
  "red",
  ieqs = [
    (2*e - spr - d0pr - 1)/3 - n_c, # blue
    (l_d - d0pr + 1) - n_c, # green
    m_c - (2*n_c + spr) - 1, # gray
    letter - 3, 4 - letter, # types C, D(even) only
    k - ncflz_CD, e - n_c + (-d0pr - spr - h_eta)/2 - k, # unsmeared sum bounds (+)
    m11 - (2*e - 2*l2 - h_eta), # m11-smear -- see p. 322 (+)
  ],
  temps = {
    l_d_term : 2*k + h_eta + d0pr,
    o_d : 0,
  },
  variables = [k],
  factor = ExponentDictionary(
    {
      G0 : 2*k + n11ceil + h_eta, G1 : -d0pr + n11ceil,
      # L0: k, L2: e - k - n11ceil - h_eta,
      q : (-m11 + d0pr + h_eta + k)
    },
  ),
  subcases = [
    n11_options, ntil_options, d0_mod_4_options, [
    RingFactorData("+"),
    RingFactorData(
      "-",
      ieqs = [
        k - ncflz_CD, e - n_c + (-d0pr - spr - h_eta)/2 - 1 - k, # unsmeared sum bounds (-)
        m11 - (2*e - 2*l2 - h_eta + 2), # m11-smear -- see p. 322 (-)
      ],
      factor = ExponentDictionary({G0 : 1, G1 : -1 }, -1),
      subcases = [normalize_g1]
      )
  ], E_options]
  )
gray = RingFactorData(
  "gray",
  ieqs = [
    (2*e - d0pr - spr) - n_c, # purple
    (2*n_c + spr) - m_c, # red
    (n_c + spr + l_d - d0pr + 1) - m_c, # green
    (e + (n_c + spr - d0pr + 1)/2) - m_c, # blue
    letter - 3, # types C, D(even), E only
  ],
  temps = {
    l_d_term : l_d,
    o_d : 0
  },
  factor = ({
      G0 : ncflz + m11_h_fl, G1 : -d0pr - ncflz  + m11_h_fl - h_eta,
      # L0: ncflz, L2 : e - m11_h_fl,
    q : (g0 - m11 - n11 + sbar/2 + d0/2 + o1_repl/4)}),
  Ftypes = ["F"],
  subcases = [n11_options, m11_options, normalize_g1]
  )

brown = RingFactorData(
"brown",
ieqs = [
  -n_c, # black
  n11 - (2*e + 1 - o1/2), # yellow
  a2bar_xi1,
  letter - 4, # types D-E only (A-C: spr is too low anyway)
],
temps = {
  vK_xi1f : vK_xi1_strong - o1/4,
  vK_xi1 : vK_xi1_strong,
  l_d_term : l_d
},
factor = ({q : (e - m11 - n11ceil_brown)}),
subcases = [[
  d0_odd.dename() + RingFactorData("F",
    factor = (
      {G0: e, G1: 0}
      ),
    Ftypes = ["F"],
    subcases = [
      [
        RingFactorData("1##", temps = {o1 : 1}),
        RingFactorData("3##", temps = {o1 : 3}),
      ]
    ]),
  RingFactorData("F+Fx",
    ieqs = [m11 - 2*e - 1],
    factor = ({
      G0 : e, G1 : e - d0pr - h_eta,
      # L0: zone_ring(0), L2: zone_ring(0)
    }),
    Ftypes = ["F", "Fx"],
    temps = {o_d : 0, o1 : o1_strong}
    ),
  RingFactorData("F_smeared", # Could be called yellow_end.
    eqns = [m11 - 2*e, n11 - 2*e],
    factor = ({
      G0 : e, G1 : e - d0pr - h_eta,
      # L0: zone_ring(0), L2: zone_ring(0)
    }),
    Ftypes = ["F"],
    temps = {o_d : 0, o1 : 2, letter : 5}
    ),
    
  RingFactorData("xF+xFx",
    factor = ({
      G0 : e, G1 : e - d0pr - h_eta,
      # L0: zone_ring(0), L2: zone_ring(0)
    }),
    Ftypes = ["xF", "xFx"],
    temps = {o_d : 0, o1 : 2 - o1_strong}
    )
  ],
  n11_brown_options
])
yellow_end = RingFactorData(
"yellow_end",
ieqs = [-n_c, # red, green
  n11 - d0, # lemon
  n11 - 1/2, # white
  2*e - 1/2 - n11, # brown
  a2bar_xi1,
  letter - 3, # types C-E only
],
temps = {
  o1 : 2,
  o2 : 0,
  o_d : 0,
  vK_xi1f : vK_xi1_strong - o1/4,
  vK_xi1 : vK_xi1_strong,
  l_d_term : dpf - 1,
},
subcases = [n11_options, [
  RingFactorData(
    "xF",
    factor = ExponentDictionary(
      {
        G0 : n11fl, G1 : -d0pr + n11fl - h_eta,
        # L0 : zone_ring(0), L2: e - n11fl,
      q : (-m11 + d0pr)}),
    temps = {letter : 4}, # type D
    Ftypes = ["xF"],
    ),
  RingFactorData(
    "Ebalx",
    ieqs = [m11 - (2*e + 1)], # m11-smear
    factor = ExponentDictionary(
      {
        G0 : 2*e - n11fl, G1 : -d0pr + n11fl - h_eta,
        # L0 : e - n11fl, L2: zone_ring(0),
      q : (-m11 + d0pr + e - n11fl)}),
    temps = {letter : 4}, # type D
    Ftypes = ["Ebalx"],
    ),
  RingFactorData(
    "F",
    factor = ExponentDictionary(
      {G0 : n11fl, G1 : -d0pr + n11fl - h_eta, q : (-m11 + d0pr)}),
    temps = {o1_strong : 2}, # types C, E
    Ftypes = ["F"],
    ),
  RingFactorData(
    "HF",
    ieqs = [
      m11 - (2*e + 1), # m11-smear
      2*e - 1 - n11
    ],
    factor = ExponentDictionary(
      {G0 : n11fl, G1 : -d0pr + n11fl - h_eta, q : (-m11 + d0pr)}),
    temps = {o1_strong : 2}, # types C, E
    Ftypes = ["eF"],
    ),
  # RingFactorData(
  #   "Fx_yellow_special",
  #   ieqs = [m11 - (2*e + 1)], # m11-smear,
  #   eqns = [n11 - 2*e],
  #   factor = ExponentDictionary(
  #     {G0 : e, G1 : e - d0pr - h_eta,
  #     # L0: zone_ring(0), L2: zone_ring(0)
  #     q : (-m11 + d0pr)}),
  #   Ftypes = ["Fx"]
  #   # If found, limit letter.
  #   )
]]
)
yellow = RingFactorData(
"yellow",
ieqs = [-n_c, # red, green
  n11 - d0 - 1, # lemon
  n11 - 1/2, # beige
  m11 - (d0)/2, # ivory
  2*e - n11, # brown
  a2bar_xi1,
  letter - 3, # types C,D(even),E only
  ],
variables = [k],
temps = {
  o1 : 0,
  o_d : 0,
  vK_xi1f : vK_xi1_strong - o1/4,
  l_d_term : 2*k + h_eta + d0pr,
},
factor = ExponentDictionary(
    {G0 : 2*k + n11ceil, G1: -d0pr + n11ceil - h_eta,
      # L0: k, L2: e - k - n11ceil,
      q : (-m11 + d0pr + k)
    }
  ),
subcases = [n11_options, [
  RingFactorData("+",
    ieqs = [
      k, e - n11ceil - k, # unsmeared sum bounds (+)
      m11 - (2*e - 2*l2 + h_eta), # m11-smear (+)
    ]
  ),
  RingFactorData(
    "-",
    ieqs = [
      k, e - n11ceil - 1 - k, # unsmeared sum bounds (-)
      m11 - (2*e - 2*l2 + 2 + h_eta), # m11-smear (-) -- see p. 322
    ],
    factor = ExponentDictionary({ G0 : 1, G1 : -1 }, -1),
  ),
  RingFactorData(
    "special", # see XIII.33
    ieqs = [
      h_eta - (m11 - n11) - (1 - n11_parity)
    ],
    eqns = [
      k + 1
    ],
    factor = ExponentDictionary({ G0 : 1, G1 : -1 }),
    )
], E_options]
)
beige_lemon = RingFactorData(
  "beige/lemon",
  ieqs = [
    -n_c, # strong
    d0 - o1/2 - n11, # yellow: see XIII.37
    m11 - (d0)/2 - 1/2, # ivory
    a2bar_xi1,
  ],
  variables = [k],
  temps = {
    l : 2*k - ((d0 - o1)/2), # level parity is governed by h_1
    l_d_term : l
  },
  factor = ({G0 : l, G1 : 0, 
    q : (-m11 + (2*l + d0 + o1)/4)
  }),
  subcases = [
    SubcaseSwitch([
      d0_even +
        RingFactorData(eqns = [flavor_type]),
      d0_odd + RingFactorData(
        ieqs = [o3 - o2 + 1, o2 - o3 + 1],
        eqns = [flavor_type - 1]),
    ], priority = 0),
    SubcaseSwitch([
      RingFactorData(
        "bottom", ieqs = [e - l], temps = {vK_xi1f : (2*l + d0 - o1)/4}),
      RingFactorData(
        "top", ieqs = [l - e - 1],
        temps = {vK_xi1f : (2*(2*e - l) + d0 - o1)/4}),
    ], priority = 0),
    E_options,
    SubcaseSwitch([
        RingFactorData("+",
          ieqs = [
            l, l - (n11 - d0pr),
            2*e - l, (2*e - l) - (n11 - d0pr), # unsmeared sum bounds
            m11 - (d0pr + l), # m11-smear
          ]
          ),
        RingFactorData("-",
          ieqs = [
            l, l - (n11 - d0pr),
            2*e - l - 2, (2*e - l) - (n11 - d0pr) - 2, # unsmeared sum bounds (-)
            m11 - (d0pr + l) - 2, # m11-smear (-)
          ],
          factor = ExponentDictionary({ G0 : 1 }, -1),
          ),
        RingFactorData("special",
          eqns = [
            m11 - n11,
            m11 - (d0pr + l) - 1
          ],
          factor = ExponentDictionary({ G0 : 1 }),
          )
    ]),
  ]
)
beige_end = RingFactorData(
"beige_end",
ieqs = [
  -n_c, # strong
  d0/2 - 1/2 - n11, # yellow, green
  m11 - (d0)/2 - 1/2, # ivory
  a2bar_xi1,
],
temps = {
  vK_xi1 : (d0 - 1)/4, # effectively letter A
  vK_xi1f : (d0 - 1 - o1)/4,
  l_d_term : -1/2,
},
factor = {q : (-m11 + (d0 - 1 + o1)/4)}, # sign error fixed: see XIII.43
subcases = [
  [
  RingFactorData("xF",
    ieqs = [letter - 2], # B-E
    factor = ExponentDictionary(
      {G0: 0, G1: 0}),
    Ftypes = ["xF"]
    ),
  RingFactorData("Ebalx",
    ieqs = [letter - 2, # B-E
      m11 - (2*e + d0pr + 1/2)], # m11-smear
    factor = ExponentDictionary(
      {G0: 2*e, G1: 0, q : e}),
    Ftypes = ["Ebalx"]
    ),
  RingFactorData("F",
    temps = {letter : 1},
    factor = ExponentDictionary(
      {G0: 0, G1: 0}),
    Ftypes = ["F"],
    ),
  RingFactorData("eF",
    temps = {letter : 1},
    ieqs = [m11 - (2*e + d0pr + 1/2)], # m11-smear
    Ftypes = ["eF"],
    factor = ExponentDictionary(
      {G0: 0, G1: 0}),
    ),
  ],
  [
    RingFactorData("0##", temps = {o1 : 0}) + d0_odd + d0pr_even,
    RingFactorData("1##", temps = {o1 : 1}) + d0_even + d0pr_odd,
    RingFactorData("2##", temps = {o1 : 2}) + d0_odd + d0pr_odd,
    RingFactorData("3##", temps = {o1 : 3}) + d0_even + d0pr_even,
  ]
])

xi1_options = [
  RingFactorData(
    "strong",
    ieqs = [
      n_c - 1/2, # weak
      2 - flavor_type - o2 + o3, 2 - flavor_type + o2 - o3,
      a2bar_xi1,
    ],
    temps = {
      o_d : mod2(o1 + o_s),
      vK_xi1f : vK_xi1_strong - o1/4,
    },
    variables = [sbarfl],
    subcases = [
      [
        RingFactorData("B-E",
          ieqs = [letter - 2], # B-E,
          temps = {
            spf : o1/2 + 2*sbarfl - dpf - o_d/2 - o_s/2,
          } # sbar mod 2 is governed by o1
          ),
        RingFactorData("A_even",
          temps = {
            letter : 1, o_d : 0,
            dpf : o1/2 + 2*sbarfl - spf - o_d/2 - o_s/2,
          },
          ),
        RingFactorData("A_odd",
          temps = {
            letter : 1, o_d : 1,
            e : o1/2 + 2*sbarfl - spf - o_d/2 - o_s/2,
          },
          ),
      ],
      [black, plum, purple, blue, greenAB, green, red, gray]
    ]
    ),
  RingFactorData(
    "weak",
    ieqs = [
      -n_c, # strong
      m11 - (d0 + 1)/2, # ivory
      a2bar_xi1,
    ],
    temps = {
      # o1_strong : o1,
      o1_strong_parity : o1*(o1 - 2)*(2*o1 - 5)/3, # = o1 % 2. o1 is fixed by now
    },
    subcases = [[brown, yellow_end, yellow]]
      ),
  beige_lemon, # o1 parity not restricted
  beige_end, # o1 parity not restricted
  RingFactorData(
    "white/ivory",
    ieqs = [
      (d0)/2 - m11, # beige/lemon
      vK_xi1f,
      vK_xi1 - (a2bar - a2),
      2*vK_xi1 - m11
    ],
    temps = {l_d_term : e}, # always big enough
    variables = [vK_xi1f],
    factor = ExponentDictionary({G0 : 0, G1 : 0, q : -vK_xi1f}, 1 - 1/q),
    subcases = [
      [
        RingFactorData(),
        RingFactorData("no_minus_1/q",
          # In flavors 002 and 020, because we rule out 1 of (q+1)-many pixels.
          # In flavors 213 and 231, because v^K(xi_1) is auto 0 anyway. 
          eqns = [
            vK_xi1f,
            o2 + o3 - o1 - 2
          ],
          factor = 1/(q - 1)
          )
      ],
      [
      RingFactorData("F",
        Ftypes = ["F"],
        eqns = [flavor_type - o1_strong_parity]
        ),
      RingFactorData("xF",
        Ftypes = ["xF"],
        eqns = [flavor_type + o1_strong_parity - 1]
        ),
    ]]
    )
]
k12pr = k12 + vK_xi1 + o2/4; # k12 is defined by a2bar option
d12 = dpf + (-2 + o3 - o1)*(1 - o_d)/4; # For M12_K. Nonlinear OK since o1, o_d are temp
add_pixel_options = SubcaseSwitch([
          RingFactorData(), # "generic" case
          RingFactorData("add_1/q", # flavor 200 only
            ieqs = [a2bar_type - 2], # a2bar:xi1,d0,a3
            eqns = [o2, k12],
            factor = 1/q,
            ),
          RingFactorData("correction",
            eqns = [o2, k12, a2bar_type - 1], # flavor 200 only. a2bar:xi2
            factor = 1/(q - 1), # undoes 1 - 1/q factor in a2bar:xi2
            )
        ])
xi2_options = [
  RingFactorData(
    "xi2free",
    ieqs = [-m12bar, -m22bar], # TODO: Once everything runs, consolidate the semifree cases.
    factor = {q : (-k12)},
    subcases = [add_pixel_options]
    ),
  RingFactorData(
    "M12_Q",
    ieqs = [
      m12bar - 1/4, -m22bar,
      a2bar - a2 - 1/4, # xi2 is Q-led
    ],
    eqns = [o1 + o2 - o3 - 2], # flavors 020, 200, 132, 312 only
    factor = {G1 : 0},
    
    # m12, d0, k12pr.  M12 => M11 => vK_xi1, k12pr defined
    subcases = [[
      RingFactorData(
        "a2bar:a3",
        temps = {a2bar_type : 4}, # special treatment
        ieqs = [vK_xi1 - m12bar],
        factor = {q : (-k12)},
        subcases = [add_pixel_options]
        ),
      RingFactorData(
        "semifree",
        ieqs = [3 - a2bar_type, # not a2bar:a3
          k12pr - m12, d0pr - m12],
        factor = {q : (-k12)},
        subcases = [add_pixel_options]
        ),
      RingFactorData(
        "M12:d0",
        ieqs = [3 - a2bar_type, # not a2bar:a3
          m12 - d0pr - 1/2, k12pr - d0pr],
          eqns = [flavor_type - (1 - o_d)], # the "d" of XII.312/XIII.53 is in Z
          factor = {q : (-k12 + d0pr - m12)}
        ),
      RingFactorData(
        "M12:k",
        ieqs = [3 - a2bar_type, # not a2bar:a3
          d0pr - k12pr - 1/2,
          m12 - k12pr - 1/2,
          a2bar_type - 2, # must not be constrained by vKxi_2
        ],
        factor = {q : (-k12 + k12pr - m12)}
        ),
      RingFactorData(
        "M12:k=d", # special zone: see XIII.53
        ieqs = [3 - a2bar_type, # not a2bar:a3
          m12 - k12pr - 1/2],
          eqns = [d0pr - k12pr - 1/2],
          temps = {a2bar_type : 1},
          factor = ExponentDictionary({q : (-k12 + k12pr - m12)}, 1/(1 - 1/q))
        )
    ]]
    ),
  RingFactorData(
    "M12_K",
    ieqs = [m12bar - 1/4, -m22bar], # M12 => M11 => vK_xi1, k12pr defined
    eqns = [
      o2 - 2*flavor_type, # flavors 002, 200, 123, 321 only
      a2bar - a2, # xi2 is K-led
    ], 
    # m12 - o1/4, vK_xi1f, d12
    subcases = [
      [
        RingFactorData(
          "semifree",
          ieqs = [vK_xi1f - (m12 - o1/4), d12 - (m12 - o1/4)],
          factor = 1,
          subcases = [add_pixel_options]
          ),
        RingFactorData(
          "M12:d0",
          ieqs = [(m12 - o1/4) - d12 - 1/4, vK_xi1f - d12],
          factor = {q : (d12 - (m12 - o1/4))},
          subcases = [[
            RingFactorData(),
            RingFactorData("correction",
              eqns = [a2bar_type - 1],
              factor = 1/(q - 1),
              # Remove the 1 - 1/q that appears in flavor 200, a2bar:xi2.
          )
          ]]
          ),
      ],
      ]
    ),
  RingFactorData(
    "M22", # See XIII.20-22.
    ieqs = [
      m22bar - 1/4, -m12bar,
      l_d_term + d0pr - m22,
      a2bar_type - 1, # flavors 020, 200, 132, 312 only
      3 - a2bar_type, # types xi1,xi2,d0 only. NOT a3.
    ],
    subcases = [
      [
        RingFactorData("semifree",
          ieqs = [(d0)/4 + l_d_term - m22/2, (a2bar - a2) - m22/2 # - 1/4
            ],
          factor = {q : (-k12)}
          ),
        RingFactorData("linear",
          ieqs = [m22/2 - (a2bar - a2), (d0)/4 + l_d_term - (a2bar - a2),
            a2bar_type - 2, # d0 and xi1 only
          ],
          factor = {q : (-m22ceil)},
          subcases = [
            m22_options
          ]
          ),
      ]
    ]
    )
]

# Unchanging requirements
main_zone = RingFactorData(
  ieqs = [
    # e and t
    t, e - t, e - 1,
    # Ordering of the a_i, a_i^bar, b_i
    a2 - a1, a3 - a2,
    a2bar_min, a2bar_d0, a2bar_a3, 
    b1, b2 - b1,
    # s'
    spr + 1/2,
    b2 - b1 - spr - 1/2,
    # d0
    dpf, e - dpf, o_d, 1 - o_d,
    # Required inactivities
    -m13bar, sbar/2 - n12bar, -n22bar,
    -m33bar
  ],
  variables = [e, t, dpf, o_d, b1f, b2f, spf, o_s,
    a1pf, a2pf, a2pbf, o1, o2, g0, g1],
  factor = ({E_ : e, T_ : t, DPFL_ : dpf,
    B1F_ : b1f, B2F_ : b2f, SPFL_ : spf,
  }),
  subcases = [
    letter_types,
    a_flavors_0,
    a2bar_options,
    xi1_options,
    xi2_options,
    distinctnesses
  ],
  
  # For testing:
  # Ftypes = ["test"]
  );

# Duals

F_shape_dual = (
  {E_ : E_*G0^2*q^2*T_, T_ : q^(-2)/T_, G0 : 1/(q*G0), G1 : G1},
  {xTH : 1/(qTH*xTH), zTH : qTH^-1*xTH^-2*zTH},
  {t : e - t, g0 : 2*e - g0});
bal_shape_dual = ({E_ : E_*q*T_, T_ : q^(-2)/T_}, {}, {t : e - t});
dual_formulas = {
  "F" : ("F", F_shape_dual),
  "xF" : ("Fx", F_shape_dual),
  "Fx" : ("xF", F_shape_dual),
  "xFx" : ("xFx", F_shape_dual),
  "Eside" : ("eF", F_shape_dual),
  "eF" : ("Eside", F_shape_dual),
  "Ebal" : ("Ebal", bal_shape_dual),
  "Ebalx" : ("Ebalx", bal_shape_dual),
}

Ftypes = dual_formulas.keys()
Ftypes_mod_symm = set(Ftypes);
for ft in Ftypes:
  if ft in Ftypes_mod_symm:
    dt = dual_formulas.get(ft)[0];
    if dt != ft:
      Ftypes_mod_symm.remove(dt)

# Spend computer cycles on the main problem.
def toil(letter_ = None, o_d_ = None, verbosity = 0):
  probe = RingFactorData();
  if letter_ is not None:
    probe += RingFactorData(temps = {letter : letter_})
  else:
    for letter_new in [1,2,3,4,5]:
      toil(letter_new, o_d_, verbosity)
    return;
  if o_d_ is not None:
    probe += RingFactorData(temps = {o_d : o_d_});
  else:
    for o_d_new in [0,1]:
      toil(letter_, o_d_new, verbosity);
    return;
      
  print ("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
  ans = run(
    hashing = "full", 
    zone = main_zone, 
    probe = probe,
    dual = True,
    add = True,
    verbosity = verbosity
  )
  return ans[0];

# Use of the resolvent data.

def value(expr):
  return deep_subs(expr, resolvent_data)
def root_dot():
  return value((-2*(e-t) + d0/4, -(e-t) + d0/4 + 1/2, -(e-t)+d0/2));
def make_strong_table():
  load("table_gen_1e21.sage");
  write_tables_e_d0pr_spr(deep_subs([(e,d0pr,spr)], resolvent_data));
def make_weak_table():
  load("table_gen_1e21_light.sage");
  write_tables_e_d0pr_h_1(deep_subs([(e,d0pr,0),(e,d0pr,2)], resolvent_data));
  
def is_asymmetric():
  if resolvent_data.get(t) is None:
    return False;
  return 2*resolvent_data.get(t) != resolvent_data.get(e)
def find_discrep(verbosity = 1):
  run_ans = run(
    hashing = "full", 
    zone = main_zone, 
    probe = RingFactorData(
      variables = [],
      temps = resolvent_data,
    ),
    asymmetric = is_asymmetric(),
    dual = True,
    add = True,
    verbosity = verbosity
  );
  return run_ans[0]; # discrepancy always comes first
def do_all_unramified_cases(verbosity = 0):
  failures = 0;
  for let in [1..5]:
    for o_d_ in [0,1]:
      print("Running", "_ABCDE"[let], o_d_);
      try:
        run_ans = run(
          hashing = "none", 
          zone = main_zone,
          probe = RingFactorData(
            temps = {
              e : 1, t : 0, dpf : 1, o_d : o_d_, letter : let
            }
            ),
          asymmetric = True,
          dual = True,
          add = True,
          verbosity = verbosity
          );
        discrep = run_ans[0];
        if discrep == {}:
          print("Passed", "_ABCDE"[let], o_d_);
        else:
          print("Failed", "_ABCDE"[let], o_d_);
          print("Discrepancy:")
          print(discrep);
          failures += 1;
      except Exception as ee:
        print(ee);
        failures += 1;
  print(failures,"failures")

def do_all_cases(verbosity = 0):
  failures = 0;
  for let in [3..5]: #$[1..5]:
    for o_d_ in [0]: #$ [0,1]:
      print("Running", "_ABCDE"[let], o_d_);
      try:
        run_ans = run(
          hashing = "full", 
          zone = main_zone,
          probe = RingFactorData(
            temps = {
              o_d : o_d_, letter : let
            }
            ),
          asymmetric = False,
          dual = True,
          add = True,
          verbosity = verbosity
          );
        discrep = run_ans[0];
        if discrep == {}:
          print("Passed", "_ABCDE"[let], o_d_);
        else:
          print("Failed", "_ABCDE"[let], o_d_);
          print("Discrepancy:")
          print(discrep);
          failures += 1;
      except Exception as ee:
        print(ee);
        failures += 1;
  print(failures,"failures")
def find_counterexample(verbosity = 1, skip_main = False):
  let = resolvent_data.get(letter);
  if let == 1:
    resolvent_data.update({spf : -1});
    vars_ = [e, t, dpf, o_d, b1f, b2f, o_s];
  elif let == 2:
    vars_ = [e, t, dpf, o_d, spf, b1f, b2f, o_s];
  elif let == 3:
    vars_ = [e, t, dpf, l_m, b1f, b2f, o_s, spf, o_d];
  elif let == 4:
    o_d_ = resolvent_data.get(o_d);
    if o_d_ == 0:
      vars_ = [e, t, dpf, sbarfl, l_m, b1f, b2f, o_s, spf];
    elif o_d_ == 1:
      vars_ = [e, t, dpf, spf, b1f, b2f, o_s];
    else:
      print("Please specify o_d.")
      return;
  elif let == 5:
    vars_ = [e, t, dpf, sbarfl, b1f, b2f, o_s, spf, o_d];
  else:
    print("Please specify letter.")
    return;
  if not skip_main:
    discrep = find_discrep(verbosity=verbosity);
    if discrep == {}:
      print("All good for resolvent", resolvent_data)
      return;
  for var in vars_:
    print("Trying values for", var)
    if resolvent_data.get(var) is None:
      for val in [0..50]: # should be long enough...
        resolvent_data.update({var : val});
        discrep = find_discrep(verbosity=verbosity);
        if discrep != {}:
          # We found a value for the variable that carries on the discrepancy.
          break;
  # All variables are set.
  if is_asymmetric():
    return zone_plot_double();
  else:
    return zone_plot();
def random_counterexample():
  max_iter = 100;
  for i in range(max_iter):
    resolvent_data.clear();
    e_ = randint(1,10)
    resolvent_data.update({e : e_});
    t_ = randint(0, e_)
    resolvent_data.update({t : t_});
    o_d_ = 0; # all odd cases done 
    resolvent_data.update({o_d : o_d_});
    o_s_ = 0;
    resolvent_data.update({o_s : o_s_});
    dpf_ = randint(1, e_);
    resolvent_data.update({dpf : dpf_});
    spf_ = randint(dpf_ - 1, dpf_ + 10);
    resolvent_data.update({spf : spf_});
    if spf_ == dpf_ - 1:
      letter_ = 3;
      h_eta_ = 1;
    elif (spf_ - dpf_) % 2 == 0:
      letter_ = 4;
      h_eta_ = 0;
    else:
      letter_ = 5;
      h_eta_ = 0;
    resolvent_data.update({letter : letter_});
    if letter_ == 5:
      continue;
    l_m_ = randint(0, floor((e_ + 1 - dpf_ + h_eta_)/2));
    resolvent_data.update({l_m : l_m_});
    b1f_ = randint(0,7);
    resolvent_data.update({b1f : b1f_});
    b2f_ = randint(b1f_ - 1, 14)
    resolvent_data.update({b2f : b2f_});
      
    discrep = zone_plot_double(verbosity = -1);
    if discrep != {}:
      zone_plot_double(verbosity = 0);
      return;
  print("Passed", max_iter, "random cases")

pieces_hard = [
  Piece(RingFactorData("{let:X}{o_d_}".format(let=let+9,o_d_=o_d_),
    temps = {o_d : o_d_, letter : let}), symm=True)
    for let in [1..5] for o_d_ in [0..1]
]

pieces_easy = [easy + piece for easy in [
  Piece(RingFactorData("e=1", temps = {e : 1, t : 0})),
  Piece(RingFactorData("red_disc=0", eqns = [red_disc]))
  ] for piece in pieces_hard
]

###########################################################
#                                                         #
#    Code below this point is designed to be modified.    #
#                                                         #
###########################################################

# The current examination case.
resolvent_data = {
  # e: 5, t: 0,
  # dpf : 2,
  o_d : 0,
  # b1f : 1, b2f : 2,
  # spf : 1,
  # o_s : 0,
  letter : 4, # Remember: D when dpf=spf (mod 2)
  # # l_m : 1, # Needed in letters C and D(even) only
}

test_zones = main_zone.fish_all(
  # "D_d0_even,020,a2bar:a3,strong,B-E,gray,>,xi2free,distinct,n11_even,m11_odd"
)
if len(test_zones) != 1:
  print("Choosing 1 of", len(test_zones), "zones");
test_zone = test_zones[0];
test_probe = RingFactorData(
      variables = [],
      temps = resolvent_data,
    )
test_zone_probe = test_zone + test_probe;

def run_it():
  ans = run(
    hashing = "full", 
    zone = test_zone,
    # zone_string = "strong",
    # Ftype_list = ["F"],
    probe = test_probe,
    asymmetric = True,
    dual = True,
    # quick = True,
    add = True,
    verbosity = 6
  )
  return ans[0];

# reset(); attach("quartic_ON_1e21.sage"); resolvent_data = {letter : 3, o_d : 0}; find_counterexample(verbosity = 2);
# reset(); attach("quartic_ON_1e21.sage"); resolvent_data = {letter : 4, o_d : 0}; find_counterexample(verbosity = 2);
# reset(); attach("quartic_ON_1e21.sage"); resolvent_data = {letter : 5, o_d : 0}; find_counterexample(verbosity = 2);
