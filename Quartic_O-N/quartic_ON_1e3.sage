attach("quartic_common.sage")


zone_ring.<e, t, b1f, b2f, o_b, a1pf, a2pf, o1, o2, l_C, k, m11dfl, n11dfl, ntil, m22dfl, efl, l0, l2> = PolynomialRing(QQ)
s = 0;


red_disc = b1f + b2f; # Note red_disc = 0 corrs to O_R
b1 = b1f + o_b/3;
b2 = b2f + 1 - o_b/3

o3 = 3 - o1 - o2;
a1p = a1pf + o1/3;
a2p = a2pf + o2/3;
a1 = a1p - 2*b1;
a2 = a2p - 2*b1;
a3 = b1 + b2 - a1 - a2 - 4*e + 4*t;
a3p = a3 + 2*b1;

h = 3 - 2*o_b;
h1 = -3 + 2*o1;
h2 = -3 + 2*o2;
h3 = -3 + 2*o3;

os = [o1, o2, o_b];
o_denom = 3;
# o_flavors = [ # For examination only
#   {o1 : 0, o2 : 1},
#   {o1 : 0, o2 : 2},
#   {o1 : 1, o2 : 2},
#   {o1 : 1, o2 : 0},
#   {o1 : 2, o2 : 0},
#   {o1 : 2, o2 : 1},
# ]
x_floor = a1pf; x_o = o1; x_exam = a1p;
y_floor = a2pf; y_o = o2; y_exam = a2p;

a = [None, a1, a2, a3];
def mij(i,j):
  return b2 - a[i] - a[j] - (2 if i == j else 3) * (e - t) + s;
def nij(i,j):
  return b1 - a[i] - a[j] - (2 if i == j else 3) * (e - t) + 2*s;

m11 = mij(1,1); m12 = mij(1,2); m22 = mij(2,2); m13 = mij(1,3); m33 = mij(3,3);
n11 = nij(1,1); n12 = nij(1,2); n22 = nij(2,2);

# d for "dot"
m11d = m11 - 2*h1/3;
n11d = n11 + 2*h1/3;
m22d = m22 - 2*h2/3;

l1 = e - l0 - l2;

def a_inact(m,n):
  return -3*mij(m,n);
def a_act(m,n):
  return 3*mij(m,n) - 1;
  
# m11dfl = floor(m11d/2)
m11d_even = RingFactorData("m11d_even", eqns = [m11d - 2*m11dfl], variables = [m11dfl])
m11d_odd = RingFactorData("m11d_odd", eqns = [m11d - 2*m11dfl - 1], variables = [m11dfl])
m11dfl_options = [m11d_even, m11d_odd];
m11dceil = m11d - m11dfl;

# n11dfl = floor(n11d/2)
n11d_even = RingFactorData("n11d_even", eqns = [n11d - 2*n11dfl], variables = [n11dfl])
n11d_odd = RingFactorData("n11d_odd", eqns = [n11d - 2*n11dfl - 1], variables = [n11dfl])
n11dfl_options = [n11d_even, n11d_odd];
n11dceil = n11d - n11dfl;

# ntil = floor((2*e - n11d + 2 + 2h)/4) = ceil((e - n11dceil + h)/2).
# Because only used in strong zones, h = h1.
ntil_near = RingFactorData("ntil_near", eqns = [e - n11dceil + h - 2*ntil], variables = [ntil])
ntil_far = RingFactorData("ntil_far", eqns = [e - n11dceil + h - 2*ntil + 1], variables = [ntil])
ntil_options = [ntil_near, ntil_far]

# m22dfl = floor(m22d/2)
m22d_even = RingFactorData("m22d_even", eqns = [m22d - 2*m22dfl], variables = [m22dfl])
m22d_odd = RingFactorData("m22d_odd", eqns = [m22d - 2*m22dfl - 1], variables = [m22dfl])
m22dfl_options = [m22d_even, m22d_odd]
m22dceil = m22d - m22dfl;

# efl = floor(e/2)
e_even = RingFactorData("e_even", eqns = [e - 2*efl], variables = [efl])
e_odd = RingFactorData("e_odd", eqns = [e - 2*efl - 1], variables = [efl])
efl_options = [e_even, e_odd]
eceil = e - efl;

# E is either F, Ebal, or Eside
E_options = [
  RingFactorData("F", ieqs = [l_C - l0, l2 - l0], Ftypes = ["F"]),
  RingFactorData("Ebal", ieqs = [l_C - l2, l0 - l2 - 1], Ftypes = ["Ebal"]),
  RingFactorData("Eside", ieqs = [l0 - l_C - 1, l2 - l_C - 1], Ftypes = ["Eside"])
]
# For testing purposes only
uniform_E = False;
if uniform_E:
  E_options = [
    RingFactorData("E", Ftypes = ["E"])
  ]


# III. Ring factors.
RINGS_RING.<E_,T_,L_C_,B1F_,B2F_,O_B,q,L0,L2> = PolynomialRing(ZZ, order='invlex');

partial_hash = {
  E_   : Fp.random_element(),
  T_   : Fp.random_element(),
  L_C_ : Fp.random_element(),
  B1F_ : Fp.random_element(),
  B2F_ : Fp.random_element(),
  O_B  : Fp.random_element(),
  q    : Fp.random_element(),
  L0 : L0h, L2 : L2h
};

full_hash = {
  E_   : Fp.random_element(),
  T_   : Fp.random_element(),
  L_C_ : Fp.random_element(),
  B1F_ : Fp.random_element(),
  B2F_ : Fp.random_element(),
  O_B  : Fp.random_element(),
  q    : Fp.random_element(),
  L0   : Fp.random_element(),
  L2   : Fp.random_element(),
};

thought = {
  q : qTH,
  E_ : yTH,
  T_ : 1,
  L_C_ : 1,
  B1F_ : 1,
  B2F_ : 1,
  O_B : 1,
  L0 : xTH/yTH,
  L2 : zTH/yTH
}

# Flavors

rsv1  = RingFactorData("h=1",  temps = {o_b : 1});
rsvm1 = RingFactorData("h=-1", temps = {o_b : 2});
o_b_flavors = [rsv1, rsvm1]

flav012 = RingFactorData("flav1", temps = {o1 : 0, o2 : 1}, factor = ExponentDictionary({q : -4/3}));
flav021 = RingFactorData("flav2", temps = {o1 : 0, o2 : 2}, factor = ExponentDictionary({q : -5/3}));
flav120 = RingFactorData("flav3", temps = {o1 : 1, o2 : 2}, factor = ExponentDictionary({q : -4/3}));
flav102 = RingFactorData("flav4", temps = {o1 : 1, o2 : 0}, factor = ExponentDictionary({q : -5/3}));
flav201 = RingFactorData("flav5", temps = {o1 : 2, o2 : 0}, factor = ExponentDictionary({q : -4/3}));
flav210 = RingFactorData("flav6", temps = {o1 : 2, o2 : 1}, factor = ExponentDictionary({q : -5/3}));
a_flavors = [
  flav012,
  flav021,
  flav120,
  flav102,
  flav201,
  flav210
]

# Distinctness / white-zone totals

distinctnesses = [
  # If all i's are distinct, as they must be.
  RingFactorData(
    ieqs = [3*(a[2] - a[1]) - 1, 3*(a[3] - a[2]) - 1],
    factor = ExponentDictionary({q : 2*a[3] - 2*a[1]})
    ),
]

# Ring factors for \xi_1

v1_options = [
  RingFactorData(
    "black",
    ieqs = [
      n11 - 2*e - 1/3 # blue
    ],
    eqns = [h - h1], # strong
    factor = ExponentDictionary({
      q : 2*e - m11 - n11 + 1,
      L0 : e, L2 : zone_ring(0)}),
    Ftypes = ["F"]
  ),
  RingFactorData(
    "blue",
    ieqs = [
      n11 - (2*e - 1 - 4*l_C + 4*h/3), # green
      2*e - n11, # black
      3*n11 - 2*e, # red
      2*m11 - (2*e + n11) - 2], # gray
    eqns = [h - h1], # strong
    factor = ExponentDictionary(
      {L0 : n11dfl, L2 : ntil, q : -m11d + 1 + e - n11dceil - ntil}),
    Ftypes = ["F"],
    subcases = [n11dfl_options, ntil_options]
  ),
  RingFactorData(
    "green",
    ieqs = [
      n11 - (2*l_C + 2 - 2*h/3), # red
      (2*e - 2 - 4*l_C + 4*h/3) - n11, # blue
      2*l_C + 1/3 - m22 + 2*o2/3 # M22. See remark.
    ], 
    eqns = [h - h1], # strong
    factor = ExponentDictionary(
      {L0: n11dfl, L2: e - l_C - n11dceil + h, q : -m11d + l_C + 1 - h}),
    subcases = [n11dfl_options, [
      RingFactorData(
        "F",
        ieqs = [m11 - (2*l_C + n11 + 1 - 2*h/3) - 1], # gray
        Ftypes = ["F"]
      ),
      RingFactorData(
        "HF",
        ieqs = [m11d - (2*e - 2*l_C)], # m11-smear
        Ftypes = ["HF"]
      )
    ]]
  ),
  RingFactorData(
    "red",
    ieqs = [
      2*e - 3*n11 - 1, # blue
      2*l_C + 1 - 2*h/3 - n11, # green 
      n11 - 1/3, # beige
      m11 - 2*n11 - 1, # gray
      2*k + 1/3 - m22 + 2*o2/3, # M22
    ],
    eqns = [h - h1], # strong
    variables = [k],
    subcases = [n11dfl_options, [
      RingFactorData(
        "+",
        ieqs = [
          k - n11dfl, e - n11d + h - k, # unsmeared sum bounds
          m11d - 2*(n11dceil - h + k), # m11-smear - see p. XII.28
        ],
        factor = ExponentDictionary(
          {L0: k, L2: e - k - n11dceil + h, q: -m11d + k + 1 - h}),
        ),
      RingFactorData(
        "-",
        ieqs = [
          k - n11dfl, e - n11d + h - 1 - k, # unsmeared sum bounds
          m11d - 2*(n11dceil - h + 1 + k), # m11-smear
        ],
        factor = ExponentDictionary(
          {L0: k + 1, L2: e - k - n11dceil + h, q: -m11d + k + 1 - h}, -1)
        )
    ], E_options]
  ),
  RingFactorData(
    "gray",
    ieqs = [
      2*e - n11, # black
      2*e + n11 + 1 - 2*m11, # blue
      2*l_C + n11 + 1 - 2*h/3 - m11, # green 
      2*n11 - m11, # red
      n11 - 1/3, # beige
    ],
    eqns = [h - h1], # strong
    factor = ExponentDictionary(
      {L0: n11dfl, L2: e - m11dfl, q: 1 - m11dceil - n11dceil}
      ),
    Ftypes = ["F"],
    subcases = [
      m11dfl_options,
      n11dfl_options
    ]
  ),
  RingFactorData(
    "beige",
    ieqs = [
      -n11, # white
      m11 - 1/3, # red
      2*k + 1/3 - m22 + 2*o2/3  # M22 restriction
    ],
    variables = [k],
    subcases = [[
      RingFactorData(
        "h1=1",
        eqns = [h1 - 1],
        subcases = [[
          RingFactorData(
            "+",
            ieqs = [
              k, e - k, # unsmeared sum bounds
              m11d - 2*k, # M11 smear
            ],
            factor = ExponentDictionary({L0: k, L2: e - k, q: -m11d + k}),
            ),
          RingFactorData(
            "-",
            ieqs = [
              k, e - 1 - k, # unsmeared sum bounds
              m11d - 2*(k + 1), # M11 smear
            ],
            factor = ExponentDictionary(
              {L0: k + (1+h)/2, L2: e - k + (h-1)/2, q: -m11d + k}, -1),
            ),
        ]]
        ),
      RingFactorData(
        "h1=-1",
        eqns = [h1 + 1],
        subcases = [[
          RingFactorData(
            "+",
            ieqs = [
              k, e - 1 - k, # unsmeared sum bounds
              m11d - (2*k + 1 - h1), # M11 smear
            ],
            factor = ExponentDictionary(
              {L0: k + (1+h)/2, L2: e - k + (h-1)/2, q: -m11d + k + 2}),
            ),
          RingFactorData(
            "-",
            ieqs = [
              k, e - 2 - k, # unsmeared sum bounds
              m11d - (2*(k + 1) + 1 - h1), # M11 smear
            ],
            factor = ExponentDictionary(
              {L0 : k + 1, L2 : e - k - 1, q: -m11d + k + 2}, -1),
            ),
        ]]
        ),
    ], E_options]
    ),
  RingFactorData(
    "beige-special", # XII.35
    ieqs = [-n11], # red
    eqns = [m11 - 1/3],
    factor = ExponentDictionary(
      {L2 : e}
      ),
    Ftypes = ["F"]
    ),
  RingFactorData(
    "white",
    ieqs = [-m11],
    factor = ExponentDictionary({L2 : e}),
    Ftypes = ["F"]
  )
]

# Ring factors for M12 and M22
v2_options = [
  RingFactorData(
    "",
    ieqs = [a_inact(1,2), a_inact(2,2)]
  ),
  RingFactorData(
    "M12_1/3",
    ieqs = [a_act(1,2), a_inact(2,2)],
    eqns = [o3 - 1],
    factor = ExponentDictionary({q : 1/3-m12})
  ),
  RingFactorData(
    "M12_2/3",
    ieqs = [a_act(1,2), a_inact(2,2)],
    eqns = [o3 - 2],
    factor = ExponentDictionary({q : 2/3-m12})
  ),
  RingFactorData(
    ieqs = [a_inact(1,2), a_act(2,2)],
    eqns = [o3 - 0],
    factor = ExponentDictionary({q : (1-h2)/2 - m22dceil}),
    subcases = [m22dfl_options]
  )
]

# Ring factors without resolvent conditions
main_zone = RingFactorData(
  ieqs = [
    # Ordering of the i_n, j_n, t, and e
    a2 - a1, a3 - a2,
    # b1,
    b2 - b1,
    t, e - t, e - 1,
    l_C, e - 2*l_C,
    # Required inactivities
    -m13, -n12, -n22,
    -m33
  ],
  variables = [e, t, b1f, b2f, o_b, l_C, a1pf, a2pf, o1, o2],
  factor = ExponentDictionary(
    {E_ : e, T_ : t, L_C_ : l_C, B1F_ : b1f, B2F_ : b2f, O_B : o_b}),
  subcases = [v1_options, v2_options, o_b_flavors, a_flavors, distinctnesses]
);

# Duals

F_shape_dual = (
  {E_ : E_*q*T_, T_ :  q^(-2)/T_, L0 : L2/q, L2 : q*L0},
  {xTH:zTH/qTH, zTH:xTH*qTH},
  {t : e - t, l0 : l2, l2 : l0});
Ebal_shape_dual = ({E_ : E_*q*T_, T_ : q^(-2)/T_}, {}, {t : e - t});
dual_formulas = {
  "F" : ("F", F_shape_dual),
  "Eside" : ("HF", F_shape_dual),
  "HF" : ("Eside", F_shape_dual),
  "Ebal" : ("Ebal", Ebal_shape_dual)
}

Ftypes = ["F", "HF", "Eside", "Ebal"]
Ftypes_mod_symm = ["F", "HF", "Ebal"]

pieces_easy = [
  Piece(RingFactorData("e=1 (part 1 of 2)", temps = {e : 1, o_b : 1}), symm = True),
  Piece(RingFactorData("e=1 (part 2 of 2)", temps = {e : 1, o_b : 2}), symm = True),
  Piece(RingFactorData("red_disc = 0 (part 1 of 2)", temps = {o_b : 1}, eqns = [red_disc]), symm = True),
  Piece(RingFactorData("red_disc = 0 (part 2 of 2)", temps = {o_b : 2}, eqns = [red_disc]), symm = True),
]

pieces_hard = [
  Piece(RingFactorData("1 of 2", temps = {o_b : 1}), symm = True),
  Piece(RingFactorData("2 of 2", temps = {o_b : 2}), symm = True),
]

###########################################################
#                                                         #
#    Code below this point is designed to be modified.    #
#                                                         #
###########################################################

# The current examination case.
resolvent_data = {
  e : 4, t : 2,
  l_C : 2,
  o_b : 1,
  b1f : 8, b2f : 11
}

def run_it():
  ans = run(
    hashing = "full", 
    zone = main_zone,
    Ftype_list = [],
    # zone_list = ["white,h=1,flav1,distinct"],
    probe = RingFactorData(
      ieqs = [
        # 2*e - 1 - s,
        # - b2
      ],
      eqns =
      dict_to_eqns(resolvent_data) +
      [ 
        # b2 - b1 - s - 1
        # e - 2,
        # t - 1,
        # lCf - 0,
        # s - 0,
        # s_parity - 1
        # b1 - 3, b2 - 7,
      ],
    ),
    dual = True,
    # quick = True,
    add = True,
    verbosity = 2
  )
  return ans[0];
