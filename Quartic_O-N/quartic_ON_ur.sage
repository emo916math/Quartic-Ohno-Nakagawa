attach("quartic_common.sage")


 # II. Zones
                                                        
zone_ring.<e, t, b1, b2, s, a1f, a2f, o1, o2, o3, lCf, k, m11fl, n11fl, ncfl, ntil, m22fl, sfl, l0, l2> = PolynomialRing(QQ);

# Basic indices
a3f = b1 + b2 - a1f - a2f - (o1 + o2 + o3)/2 - 4*e + 4*t + 2*s;
red_disc = b1 + b2; #
a1 = a1f + o1/2;
a2 = a2f + o2/2;
a3 = a3f + o3/2;
a = [None, a1, a2, a3];

# For examination
os = [o1, o2, o3];
x_floor = a1f; x_o = o1; x_exam = a1;
y_floor = a2f; y_o = o2; y_exam = a2;
o_denom = 2;
# o_flavors = [
#   {o1 : 0, o2 : 0, o3 : 0},
#   {o1 : 0, o2 : 1, o3 : 1},
#   {o1 : 1, o2 : 0, o3 : 1},
#   {o1 : 1, o2 : 1, o3 : 0},
# ]

def mij(i,j):
  return b2 - a[i] - a[j] - (2 if i == j else 3) * (e - t) + s;
def nij(i,j):
  return b1 - a[i] - a[j] - (2 if i == j else 3) * (e - t) + 2*s;
  
m11 = mij(1,1); m12 = mij(1,2); m22 = mij(2,2); m13 = mij(1,3); m33 = mij(3,3);
n11 = nij(1,1); n12 = nij(1,2); n22 = nij(2,2);
n_c = n11 - s; # n11 according to the colorful normalization

l1 = e - l0 - l2;

def a_inact(m,n):
  return -2*mij(m,n);
def a_act(m,n):
  return 2*mij(m,n) - 1;

# Parity options.

# sfl = floor(s/2)
s_even = RingFactorData(
  "s_even",
  eqns = [s - 2*sfl],
  ieqs = [lCf, e - 2*lCf]
)
s_odd = RingFactorData(
  "s_odd",
  eqns = [s - 2*sfl - 1, lCf + 1]
)
s_options = [s_even, s_odd]
sceil = s - sfl;
s_parity = s - 2*sfl;

# m11fl = floor(m11/2)
m11_even = RingFactorData("m11_even", eqns = [m11 - 2*m11fl], variables = [m11fl])
m11_odd = RingFactorData("m11_odd", eqns = [m11 - 2*m11fl - 1], variables = [m11fl])
m11_options = [m11_even, m11_odd];
m11ceil = m11 - m11fl;

# n11fl = floor(n11/2)
n11_even = RingFactorData("n11_even", eqns = [n11 - 2*n11fl], variables = [n11fl])
n11_odd = RingFactorData("n11_odd", eqns = [n11 - 2*n11fl - 1], variables = [n11fl])
n11_options = [n11_even, n11_odd];
n11ceil = n11 - n11fl;

# ncfl = floor(n_c/2)
nc_even = RingFactorData("nc_even", eqns = [n_c - 2*ncfl], variables = [ncfl])
nc_odd = RingFactorData("nc_odd", eqns = [n_c - 2*ncfl - 1], variables = [ncfl])
nc_options = [nc_even, nc_odd];
ncceil = n_c - ncfl;
nc_parity = n_c - 2*ncfl;

# m22fl = floor(m22/2)
m22_even = RingFactorData("m22_even", eqns = [mij(2, 2) - 2*m22fl], variables = [m22fl])
m22_odd = RingFactorData("m22_odd", eqns = [mij(2, 2) - 2*m22fl - 1], variables = [m22fl])
m22fl_options = [m22_even, m22_odd]
m22ceil = m22 - m22fl;

# For s even, ntil = floor((2*e - n11 + 2)/4) = ceil((e - n11ceil)/2). Always used with nc_options.
ntil_near = s_even + RingFactorData("ntil_near", eqns = [e - ncceil - sfl - 2*ntil],     variables = [ntil])
ntil_far  = s_even + RingFactorData("ntil_far",  eqns = [e - ncceil - sfl - 2*ntil + 1], variables = [ntil])
ntil_options = [ntil_near, ntil_far]

# E is either F, Ebal, Eside, or exceptionally xF
E_options = [
  RingFactorData("F", ieqs = [lCf - l0, l2 - l0], Ftypes = ["F"]),
  RingFactorData("Ebal", ieqs = [lCf - l2, l0 - l2 - 1], Ftypes = ["Ebal"]),
  RingFactorData("Eside", ieqs = [l0 - lCf - 1, l2 - lCf - 1,
    lCf + l0], # xF
    Ftypes = ["Eside"]),
  RingFactorData("xF", ieqs = [s - 2*l1], eqns = [lCf + 1, l0 - 0],
        Ftypes = ["xF"]),
  RingFactorData("Eside_special",
    ieqs = [2*l1 - s - 1], eqns = [lCf + 1, l0 - 0],
    Ftypes = ["Eside"]),
]
# For testing purposes only
uniform_E = False;
if uniform_E:
  E_options = [
    RingFactorData("E", Ftypes = ["E"])
  ]


# III. Ring factors.
RINGS_RING.<E_,T_,LCF_,SFL_,B1_,B2_,q,L0,L2> = PolynomialRing(ZZ, order='invlex');


PartialHash.<L0h,L2h> = PolynomialRing(Fp);
partial_hash = {
  E_   : Fp.random_element(),
  T_   : Fp.random_element(),
  SFL_ : Fp.random_element(),
  LCF_ : Fp.random_element(),
  B1_  : Fp.random_element(),
  B2_  : Fp.random_element(),
  q    : Fp.random_element(),
  L0 : L0h, L2 : L2h
};

full_hash = {
  E_   : Fp.random_element(),
  T_   : Fp.random_element(),
  SFL_ : Fp.random_element(),
  LCF_ : Fp.random_element(),
  B1_  : Fp.random_element(),
  B2_  : Fp.random_element(),
  q    : Fp.random_element(),
  L0   : Fp.random_element(),
  L2   : Fp.random_element(),
};

# For examination of ring answers in individual zones
# q_ring.<q_> = PolynomialRing(ZZ)
# THOUGHT.<x,y,z> = PolynomialRing(q_ring, order='invlex')

thought = {
  q : qTH,
  E_ : yTH,
  T_ : 1,
  LCF_ : 1,
  SFL_ : 1,
  B1_ : 1,
  B2_ : 1,
  L0 : xTH/yTH,
  L2 : zTH/yTH
}

flav0 = RingFactorData("flav0", temps = {o1 : 0, o2 : 0, o3 : 0});
flav1 = RingFactorData("flav1", temps = {o1 : 0, o2 : 1, o3 : 1}, factor = ExponentDictionary({}, 1/q));
flav2 = RingFactorData("flav2", temps = {o1 : 1, o2 : 0, o3 : 1}, factor = ExponentDictionary({}, 1/q));
flav3 = RingFactorData("flav3", temps = {o1 : 1, o2 : 1, o3 : 0}, factor = ExponentDictionary({}, 1/q));
a_flavors = [flav0, flav1, flav2, flav3]

flav0_temps_only = RingFactorData("", temps = flav0.temps);
flav1_temps_only = RingFactorData("", temps = flav1.temps);
flav2_temps_only = RingFactorData("", temps = flav2.temps);
flav3_temps_only = RingFactorData("", temps = flav3.temps);
off_flavors_temps_only = [flav1_temps_only, flav2_temps_only, flav3_temps_only]

distinctnesses = [
  # If all i's are distinct. 
  RingFactorData(
    "distinct",
    ieqs = [2*(a[2] - a[1]) - 1, 2*(a[3] - a[2]) - 1],
    factor = ExponentDictionary({q : 2*a[3] - 2*a[1]})
    ),
  
  # If a1 == a2 < a3.
  RingFactorData(
    "a1==a2",
    ieqs = [2*(a[3] - a[2]) - 1],
    eqns = [a[2] - a[1]],
    factor = ExponentDictionary({q : 2*a[3] - 2*a[1]}, 1/(1 + 1/q))
    ),
  # If a1 < a2 == a3.
  RingFactorData(
    "a2==a3",
    ieqs = [2*(a[2] - a[1]) - 1],
    eqns = [a[3] - a[2]],
    factor = ExponentDictionary({q : 2*a[3] - 2*a[1]}, 1/(1 + 1/q))
    ),
  # If all i's are equal.
  RingFactorData(
    "a1==a2==a3",
    eqns = [a2 - a1, a3 - a2],
    factor = ExponentDictionary({}, 1/(1 + 1/q)/(1 + 1/q + 1/q^2))
    )
]

# Ring factors for \xi_1
strong_s_options = [s_even + flav0_temps_only, s_odd + flav2_temps_only, s_odd + flav3_temps_only]

v1_options = [
  RingFactorData(
    "black", 
    ieqs = [n_c - 2*e - 1],
    factor = ExponentDictionary({L0 : e, L2 : zone_ring(0), q : 2*e - m11 - n_c - sfl}),
    Ftypes = ["F", "Fx", "Fxx"],
    subcases = [strong_s_options]
    ),
  RingFactorData(
    "purple",
    ieqs = [
      2*e - n_c, # black
      n_c - (2*e - s) - 1, # blue/green/gray
      n_c - 1], # brown
    factor = ExponentDictionary({L0: ncfl, L2 : zone_ring(0), q : e - m11 - ncceil - sfl}),
    Ftypes = ["F", "Fx"],
    subcases = [strong_s_options, nc_options]
    ),
  RingFactorData(
    "blue",
    ieqs = [n_c - (2*e - s - 4*lCf - 2) - 1, # green
      (2*e - s) - n_c, # purple
      3*n_c - (2*e - s), # red
      # 2*lCf + 1 - m22, # M22 is automatic
    m11 - (e + ncceil + sfl) - 1], # gray
    factor = ExponentDictionary({L0 : ncfl, L2 : ntil, q : -m11 + e - sfl - ncceil - ntil}),
    Ftypes = ["F"],
    subcases = [nc_options, ntil_options]
    ) + s_even + flav0_temps_only,
  RingFactorData(
    "green",
    ieqs = [
      (2*e - s) - n_c, # purple
      n_c - (2*lCf + 1) - 1, # red
      (2*e - s - 4*lCf - 2) - n_c, # blue
      n_c - 1, # yellow
      2*lCf + 1 + s_parity - m22, # M22
      m11 - (2*ncceil + 2*lCf + s + 1) - 1 # gray
    ], 
    factor = ExponentDictionary({L0: ncfl, L2: e - lCf - sceil - ncceil, q : -m11 + lCf + s_parity}),
    subcases = [strong_s_options, nc_options, [
      RingFactorData(
        "HF",
        ieqs = [m11 - (2*e - 2*lCf - s_parity)], # HF smears away
        Ftypes = ["HF"]
        ),
      RingFactorData( 
        "F",
        Ftypes = ["F"]
        )
    ]]
    ),
  RingFactorData(
    "gray",
    ieqs = [
      2*e - m11, # purple
      ncceil + e + sfl - m11, # blue
      (2*ncceil + 2*lCf + s + 1) - m11, # green
      2*n_c + s - m11 - 1, # red
      n_c - 1 # yellow
      ], 
    factor = ExponentDictionary({L0: ncfl, L2: e - m11fl, q: -m11ceil - ncceil - sfl}),
    Ftypes = ["F"],
    subcases = [strong_s_options, m11_options, nc_options]
    ),
  RingFactorData(
    "red",
    ieqs = [
      (2*e - s) - 3*n_c - 1, # blue 
      (2*lCf) + 1 - n_c, # green
      n_c - 1, # yellow 
      m11 - 2*n_c - s, # gray
      2*k + 1 - m22, # m22 restriction
    ],
    variables = [k],
    subcases = [n11_options, [
      RingFactorData(
        "+",
        ieqs = [
          k - (n11fl - sfl), e - n11 + sfl - k, # unsmeared sum bounds
          m11 - 2*(n11ceil + k), # m11-smear: k <= floor(m11/2) - ceil(n11/2)
        ],
        factor = ExponentDictionary({L0: k, L2: e - k - n11ceil, q: -m11 + k}),
        ),
      RingFactorData(
        "-",
        ieqs = [
          k - (n11fl - sfl), e - n11 + sfl - 1 - k, # unsmeared sum bounds
          m11 - 2*(n11ceil + 1 + k), # m11-smear: k <= floor(m11/2) - ceil(n11/2) - 1
        ],
        factor = ExponentDictionary({L0: k + 1, L2: e - k - n11ceil, q: -m11 + k}, -1),
        )
    ], E_options]
    ) + s_even + flav0_temps_only,
  RingFactorData(
    "red,special",
    ieqs = [
      (2*e - s) - 3*n_c - 1, # blue 
      (2*lCf) + 1 - n_c, # green
      n_c - 1, # yellow 
      m11 - 2*n_c - s, # gray
    ],
    eqns = [n22], 
    factor = ExponentDictionary({L0: m22fl, L2: e - m11fl, q: -m11 + m22fl - 1}),
    subcases = [E_options] # actually, only F and Eside occur.
    ) + s_even + flav0_temps_only + m22_odd + nc_odd + m11_even,
  RingFactorData(
    "brown",
    ieqs = [n11 - 2*e - 1, # yellow
    s - n11], # purple
    subcases = [n11_options, [
      RingFactorData(
        subcases = [[
          RingFactorData(
            factor = ExponentDictionary({L0: zone_ring(0), L2: zone_ring(0),
            q: e - m11 - n11ceil}), # ceiling option for flavor 0
            Ftypes = ["F", "Fx"]
            ) + s_even,
          RingFactorData(
            factor = ExponentDictionary({L0: zone_ring(0), L2: zone_ring(0),
            q: e - m11 - n11ceil}), # ceiling option for flavor 0
            Ftypes = ["xF", "xFx"]
            ) + s_odd
        ]]
        ) + flav0_temps_only,
      RingFactorData(
        eqns = [o1 - 1],
        subcases = [[
          RingFactorData(
            factor = ExponentDictionary({L0: zone_ring(0), L2: zone_ring(0),
            q: e - m11 - n11fl}), # floor option for flavors 2 and 3
            Ftypes = ["xF", "xFx"]
            ) + s_even,
          RingFactorData(
            factor = ExponentDictionary({L0: zone_ring(0), L2: zone_ring(0),
            q: e - m11 - n11fl}), # floor option for flavors 2 and 3
            Ftypes = ["F", "Fx"]
            ) + s_odd
        ]]
        ),
    ]]
    ),
  RingFactorData(
    "yellow",
    ieqs = [-n_c, # red
      n11 - 1, # beige
      2*e - n11, # brown
    ],
    subcases = [n11_options, s_options, [
      RingFactorData(
        "+",
        ieqs = [
          k, e - n11ceil - k, # unsmeared sum bounds
          m11 - 2*(n11ceil + k), # m11-smear: k <= floor(m11/2) - ceil(n11/2)
          2*k + 1 - m22, # m22 restriction
        ],
        variables = [k],
        factor = ExponentDictionary({L0: k, L2: e - k - n11ceil, q: -m11 + k}),
        ),
      RingFactorData(
        "-",
        ieqs = [
          k, e - n11ceil - k - 1, # unsmeared sum bounds
          m11 - 2*(n11ceil + 1 + k), # m11-smear: k <= floor(m11/2) - ceil(n11/2) - 1
          2*k + 1 - m22, # m22 restriction
        ],
        variables = [k],
        factor = ExponentDictionary({L0: k + 1, L2: e - k - n11ceil, q: -m11 + k}, -1),
        ),
    ], E_options]
    ) + flav0_temps_only,
  RingFactorData(
    "yellow_special",
    ieqs = [-n_c, # red
      n11 - 1, # beige
      2*e - n11, # brown
    ],
    eqns = [m11 - n11],
    factor = ExponentDictionary({L0: zone_ring(0), L2: e - n11fl, q: -n11 - 1}),
    subcases = [[
      s_even + RingFactorData(Ftypes = ["F"]),
      s_odd + RingFactorData(Ftypes = ["xF"]),
    ]]
    ) + n11_odd + flav0_temps_only,
  RingFactorData(
    "yellow_end",
    ieqs = [-n_c, # red
      n11 - 1, # beige
      2*e - n11, # brown
      -m22, # M22 must be off
      o1 - 1, # tiny conic
    ],
    subcases = [n11_options, [
      RingFactorData(
        "xF",
        factor = ExponentDictionary(
          {L0: zone_ring(0), L2: e - n11fl, q: -m11}),
        Ftypes = ["xF"]
        ) + s_even,
      RingFactorData(
        "Ebalx",
        ieqs = [m11 - (2*e + 1)], # m11-smear
        factor = ExponentDictionary(
          {L0: e - n11fl, L2: zone_ring(0), q: -m11 + e - n11fl}),
        Ftypes = ["Ebalx"]
        ) + s_even,
      RingFactorData(
        "F",
        factor = ExponentDictionary(
          {L0: zone_ring(0), L2: e - n11fl, q: -m11}),
        Ftypes = ["F"]
        ) + s_odd,
      RingFactorData(
        "HF",
        ieqs = [
          m11 - (2*e + 1), # m11-smear
          2*e - 1 - n11
        ],
        factor = ExponentDictionary(
          {L0: zone_ring(0), L2: e - n11fl, q: -m11}),
        Ftypes = ["HF"]
        ) + s_odd,
      RingFactorData(
        "Fx_yellow_special",
        ieqs = [m11 - (2*e + 1)], # m11-smear,
        eqns = [n11 - 2*e],
        factor = ExponentDictionary(
          {L0: zone_ring(0), L2: zone_ring(0), q: -m11}),
        Ftypes = ["Fx"]
        ) + s_odd
    ]]
    ),
  RingFactorData(
    "beige",
    ieqs = [-n11, m11 - 1], # beige zone
    variables = [k],
    subcases = [s_options, [
      RingFactorData(
        "+",
        ieqs = [
          k, e - k, m11 - 2*k, # unsmeared sum bounds
          2*k + 1 - m22, # m22 restriction
        ],
        factor = ExponentDictionary({L0: k, L2: e - k, q: -m11 + k}),
        ) + flav0_temps_only,
      RingFactorData(
        "-",
        ieqs = [
          k - 1, e - k - 1, m11 - 2*(k + 1), # unsmeared sum bounds
          2*(k - 1) + 1 - m22, # m22 restriction
        ],
        factor = ExponentDictionary({L0: k, L2: e - k, q: -m11 + k - 1}, -1),
        ) + flav0_temps_only
    ], E_options]
    ),
  RingFactorData(
    "beige_special_1",
    ieqs = [-n11, m11 - 1], # beige zone
    eqns = [e],
    factor = ExponentDictionary({q: -m11 - 1}),
    subcases = [s_options, E_options] # actually, only F (s even) and xF (s odd)
    ) + flav0_temps_only,
  RingFactorData(
    "beige_special_2",
    ieqs = [-n11, # beige zone
      e - 1, # so as not to duplicate beige_special_1
    ],
    eqns = [m11 - 1],
    factor = ExponentDictionary({L2: e}, q^-2),
    subcases = [s_options, E_options] # actually, only F (s even) and xF (s odd)
    ) + flav0_temps_only,
  RingFactorData(
    "beige_special_3",
    ieqs = [-n11, # beige zone
      e - 1, # so as not to duplicate beige_special_1
      m11 - 2, # so as not to duplicate beige_special_2
    ],
    eqns = [a1 - a2],
    factor = ExponentDictionary({L0: m11fl, L2: e - m11fl, q: -m11ceil - 1}),
    subcases = [s_options, m11_options, E_options]
    ) + flav0_temps_only,
  RingFactorData(
    "beige,end",
    ieqs = [-n11, m11 - 1, # beige zone
    ],
    eqns = [o1 - 1],
    subcases = [[
      RingFactorData(
        factor = ExponentDictionary({L0: 0, L2: e, q: -m11}),
        Ftypes = ["xF", "xxF"]
        ) + s_even,
      RingFactorData(
        ieqs = [m11 - (2*e + 1)], # m11-smear
        factor = ExponentDictionary({L0: e, L2: 0, q: -m11 + e}),
        Ftypes = ["Ebalx", "Ebalxx"]
        ) + s_even,
      RingFactorData(
        factor = ExponentDictionary({L0: 0, L2: e, q: -m11}),
        Ftypes = ["F", "xxF"]
        ) + s_odd,
      RingFactorData(
        ieqs = [m11 - (2*e + 1)], # m11-smear
        factor = ExponentDictionary({L0: 0, L2: e, q: -m11}),
        subcases = [[
          RingFactorData(ieqs = [e-1], Ftypes = ["HF"]),
          RingFactorData(eqns = [e], Ftypes = ["Fx"])
        ]]
        ) + s_odd,
      RingFactorData(
        ieqs = [m11 - (2*e + 1)], # m11-smear
        factor = ExponentDictionary({L0: e, L2: 0, q: -m11 + e}),
        Ftypes = ["Ebalxx"]
        ) + s_odd
    ]]
    ),
  RingFactorData(
    "white",
    ieqs = [-m11],
    subcases = [[
      RingFactorData(
        factor = ExponentDictionary({L2 : e}, (1 + 1/q + 1/q^2)),
        Ftypes = ["F"]
        ) + s_even + flav0_temps_only,
      RingFactorData(
        factor = ExponentDictionary({L2 : e}, (1 + 1/q)),
        Ftypes = ["xF", "xxF"],
        subcases = [off_flavors_temps_only]
        ) + s_even,
      RingFactorData(
        factor = ExponentDictionary({L2 : e}, (1 + 1/q + 1/q^2)),
        Ftypes = ["xF"]
        ) + s_odd + flav0_temps_only,
      RingFactorData(
        factor = ExponentDictionary({L2 : e}, (1 + 1/q)),
        Ftypes = ["F", "xxF"],
        subcases = [off_flavors_temps_only]
        ) + s_odd
    ]]
    )
];

# Ring factors for \xi_2: zero flavor, then nonzero flavor.
v2_options = [
  RingFactorData(
    "",
    ieqs = [a_inact(1,2), a_inact(2,2)],
    factor = ExponentDictionary({}, 1 + 1/q)
    ) + flav0_temps_only,
  RingFactorData(
    "M12_on",
    ieqs = [a_act(1,2), a_inact(2,2)],
    factor = ExponentDictionary({q : -m12})
    ) + flav0_temps_only,
  RingFactorData(
    "M22_on",
    ieqs = [a_inact(1,2), a_act(2,2)],
    factor = ExponentDictionary({q : -m22ceil}),
    subcases = [m22fl_options]
    ) + flav0_temps_only,
  RingFactorData(
    "",
    ieqs = [a_inact(1,2), a_inact(2,2)],
    factor = ExponentDictionary({}, 1),
    subcases = [off_flavors_temps_only]
    ),
  RingFactorData(
    "M12_on",
    ieqs = [a_act(1,2), a_inact(2,2)],
    factor = ExponentDictionary({q : 1/2 - m12})
    ) + flav2_temps_only,
]


# There are no v3_options because \xi_3' has soln volume 1 in all cases.

# IV. Fourier transform

F_shape_dual = (
  {E_ : E_*q*T_, T_ :  q^(-2)/T_, L0 : L2/q, L2 : q*L0},
  {xTH:zTH/qTH, zTH:xTH*qTH},
  {t : e - t, l0 : l2, l2 : l0});
Ebal_shape_dual = ({E_ : E_*q*T_, T_ : q^(-2)/T_}, {}, {t : e - t});
dual_formulas = {
  "F" : ("F", F_shape_dual),
  "xF" : ("Fx", F_shape_dual),
  "xxF" : ("Fxx", F_shape_dual),
  "Fx" : ("xF", F_shape_dual),
  "Fxx" : ("xxF", F_shape_dual),
  "xFx" : ("xFx", F_shape_dual),
  "Eside" : ("HF", F_shape_dual),
  "HF" : ("Eside", F_shape_dual),
  "Ebal" : ("Ebal", Ebal_shape_dual),
  "Ebalx" : ("Ebalx", Ebal_shape_dual),
  "Ebalxx" : ("Ebalxx", Ebal_shape_dual),
}

Ftypes = dual_formulas.keys()
Ftypes_mod_symm = set(Ftypes);
for ft in Ftypes:
  if ft in Ftypes_mod_symm:
    dt = dual_formulas.get(ft)[0];
    if dt != ft:
      Ftypes_mod_symm.remove(dt)

# Unchanging requirements
main_zone = RingFactorData(
  ieqs = [
    # Ordering of the a_i, b_i, t, and e
    a2 - a1, a3 - a2,
    b1, b2 - b1,
    t, e - t, # e - 1,
    # s
    s - 0, b2 - b1 - s,
    # Required inactivities
    -m13, s/2 - n12, -n22,
    -m33
  ],
  variables = [e, t, b1, b2, s, sfl, lCf, a1f, a2f, o1, o2, o3],
  factor = ExponentDictionary({E_ : e, T_ : t, LCF_ : lCf, SFL_ : sfl,
    B1_ : b1, B2_ : b2}),
  subcases = [v1_options, v2_options, a_flavors, distinctnesses]
  );

pieces_easy = [
  Piece(RingFactorData("tame", temps = {e : 0}), symm=True),
  Piece(RingFactorData("e=1", temps = {e : 1}), symm=True),
  Piece(RingFactorData("red_disc=0", eqns = [red_disc]), symm=True)
]

pieces_hard = [
  Piece(RingFactorData(), symm=True)
]

###########################################################
#                                                         #
#    Code below this point is designed to be modified.    #
#                                                         #
###########################################################

# The current examination case.
resolvent_data = {
  # e : 8, t : 4,
  # lCf : 0,
  # s : 4,
  # b1 : 4, b2 : 8
}

# The current test case.
def run_it():
  ans = run(
    hashing = "full", 
    zone = main_zone,
    Ftype_list = ["HF"],
    probe = RingFactorData(
      ieqs = [
        # 2*e - 1 - s,
        # - b2
      ],
      eqns =
      dict_to_eqns(resolvent_data) + [ 
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
    asymmetric = False,
    verbosity = 2
  )
  return ans[0];
