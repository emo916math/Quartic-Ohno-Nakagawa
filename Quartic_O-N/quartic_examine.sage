# Items to be defined in the respective splitting types:
# "resolvent_data": A dictionary including values for e, t, and the 
# resolvent parameters in the appropriate splitting type.
# "superzone": A RingFactorData object.

import re

def dict_to_eqns(dic):
  return [var - val for (var,val) in dic.items()]

# Change colors from those that have simple names
# (to be used in the paper) to
# a visually distinct palette.
colors = {
  "black" : "black",
  "plum" : "plum",
  "purple" : "purple",
  "brown" : "tan",
  "blue" : "blue",
  "green" : "green",
  "red" : "red",
  "yellow" : "gold",
  "beige" : "antiquewhite",
  "gray" : "darkgray",
  "white" : "white",
}

def get_color(zone_name):
  blocks = re.split("[^a-z0-9]+", zone_name);
  for block in blocks:
    col = colors.get(block);
    if col is not None:
      return col;
  return "none";

def point_converter(point):
  if len(point) == 2:
    x, y = point;
    x_floor_ = floor(x);
    y_floor_ = floor(y);
    x_o_ = o_denom * (x - x_floor_);
    y_o_ = o_denom * (y - y_floor_);
    return dict_add(resolvent_data,
          {x_floor : x_floor_, y_floor : y_floor_, x_o : x_o_, y_o : y_o_});
  else:
    x, y, z = point;
    x_floor_ = floor(x);
    y_floor_ = floor(y);
    z_floor_ = floor(z);
    x_o_ = o_denom * (x - x_floor_);
    y_o_ = o_denom * (y - y_floor_);
    z_o_ = o_denom * (z - z_floor_);
    return dict_add(resolvent_data,
          {x_floor : x_floor_, y_floor : y_floor_, z_floor : z_floor_,
          x_o : x_o_, y_o : y_o_, z_o : z_o_});

# Converts values to mixed fraction form
def mixed_frac(val, separator=" + "):
  if val in ZZ:
    return str(val);
  elif val in QQ:
    return str(floor(val)) + separator + str(val - floor(val))
  else:
    return str(val);

# arg "point": a pair (x,y), (x,y,z), or (with row) (y,z)
def point_total(point, row = None, added_probe=None, zone_string = "", verbosity=0):
  if not(added_probe):
    added_probe = RingFactorData();
  point_full = (row, ) + tuple(point) if row else point;
  point_str = "(" + ", ".join(mixed_frac(z) for z in point_full) + ")";
  if verbosity > 0:
    print(point_str);
  total, zones = run(
    zone = main_zone,
    hashing = "thought",
    add = True,
    zone_string = zone_string,
    probe = RingFactorData(
      temps = point_converter(point_full)
      ) + added_probe,
    verbosity = verbosity)
  total = {key : val.evaluate() for (key, val) in total.items()}
  if len(zones) > 0 or verbosity > 0:
    if verbosity <= 0:
      print(point_str);
    print(zones);
    print(shorten(total))
    print();
  return total, zones;

def inspect(point, expr=None, verbosity=0, zone_string = ""):
  cvtr = point_converter(point);
  if expr is not None:
    return deep_subs(expr, cvtr)
  print(cvtr);
  for var in [
      "e", "t", "d0pr", "l_C", "l_m",
      "s", "spr", "sbar", "b1", "b2", "o_b", "h",
      "",
      "a1", "a2", "a3", "a2bar", "a3bar",
      "",
      "a1p", "a2p", "a3p", "a2pbf",
      "",
      "h1", "h2",
      "m11", "n11", "m_c", "n_c",
      "m12", "m12bar",
      "m22", "m22bar",
      "",
      "m13", "m13bar",
      "n12", "n12bar",
      "n22", "n22bar",
      "m33", "m33bar",
  ]:
    
    if var == "": # separators
      print();
      continue;
    try:
      val = eval(var);
      val = val.subs(cvtr);
    except NameError:
      # variable is inapplicable to this splitting type
      continue;
    if val in QQ:
      print(var, "=", mixed_frac(val));
    
  point_total(point, verbosity=verbosity,
    zone_string=zone_string);
  


plot_ring.<x_floor_plot,y_floor_plot> = PolynomialRing(QQ)

# Output: A polygon, the convex hull of possibly admissible pairs (a1,a2);
# and the points satisfying a finer version of this.
def pgn_and_pts(row = None):
  assert row is None; # for compatibility across splitting types.
  
  # x_exam, etc. are defined in the respective zones.
  assert x_exam == x_floor + x_o / o_denom;
  assert y_exam == y_floor + y_o / o_denom;
  
  plot_subs = dict_add(resolvent_data,
    {x_floor : x_floor_plot, y_floor : y_floor_plot, x_o : 0, y_o : 0})
  if False:
    print("plot_subs is now", plot_subs)
    for ieq in main_zone.ieqs + main_zone.eqns:
      print("Step 0:", ieq)
      print("Step 1:", ieq.subs(plot_subs))
      print("Step 2:", plot_ring(ieq.subs(plot_subs)));
      print("Step 3:", linear_to_coeffs(plot_ring(ieq.subs(plot_subs))));
  plot_ieqs = [linear_to_coeffs(plot_ring(ieq.subs(plot_subs))) for ieq in main_zone.ieqs]
  plot_eqns = [linear_to_coeffs(plot_ring(eqn.subs(plot_subs))) for eqn in main_zone.eqns]
  pgn_0 = Polyhedron(ieqs = plot_ieqs, eqns = plot_eqns);
  pts = [];
  for a_flav in a_flavors:
    vec = (QQ^2)((a_flav.temps.get(x_o)/o_denom,
                  a_flav.temps.get(y_o)/o_denom))
    pgn = pgn_0.translation(-vec)
    pts += [pt + vec for pt in pgn.integral_points() if
      plot_ring(m12.subs(plot_subs))(*(pt + vec)) <= 0 or
      plot_ring(m22.subs(plot_subs))(*(pt + vec)) <= 0
    ]
  pts.sort();
  return pgn_0, pts;
  
def zone_plot(row = None, verbosity = 0):
  t0 = time.time();
  pgn, pts = pgn_and_pts(row = row);
  print("Plotting resolvent data", resolvent_data);
  print("Drawing outline and", len(pts), "points");
  print();
  
  if pgn.ambient_dim() == 2:
    graphics = pgn.projection().render_outline_2d(
      aspect_ratio = 1, color='silver',
      axes_labels=
        ["$n_{11}$", "$m_{11}$"]
        if x_floor == n11 and y_floor == m11 else
        ["$a_1'$", "$a_2'$"] if row is None else
        ["$" + str(y_floor) + "$",
         "$" + str(z_floor) + "$"],
      axes_labels_size=1,
      title = '$' + ', '.join(
        str(var) + "=" + str(val) for (var,val) in resolvent_data.items()
      ) + (", a_1 = " + str(row) if row is not None else "") + '$'
    );
  else:
    # Occurs in spl. t. 1^2 1.
    graphics = pgn.projection().render_wireframe_3d(
      color='silver',
      axes_labels=[str(x_floor),
                   str(y_floor),
                   str(z_floor)],
      axes_labels_size=1,
      ticks = 1,
      title = '$' + ', '.join(
        str(var) + "=" + str(val) for (var,val) in resolvent_data.items()
      ) + '$'
    );
  
  total = {}; row_total = {}; num_dots = 0;
  old_pt = None; 
  for pt in pts:
    if old_pt is not None and pt[0] != old_pt[0] and row_total != {}:
      print("Total for row (" + mixed_frac(old_pt[0]) + "):")
      print(shorten(row_total));
      print("---------------------------------------");
      row_total = {};
    
    ans, zone_names = point_total(pt, row, verbosity=verbosity);
    try:
      zone_name = zone_names[0];
    except IndexError:
      continue;
    num_dots += 1;
    color = get_color(zone_name)
    for zone_name2 in zone_names:
      color2 = get_color(zone_name2)
      if color2 != color:
        print("Multiple colors for point", pt, ":", color, ",", color2)
        inspect(pt);
        # return;
    if len(pt) == 2:
      graphics += circle(center=pt, radius=1/(2*o_denom),
        fill=True, edgecolor = "black",
        facecolor=color)
    else:
      graphics += sphere(center=pt, size=1/(2*o_denom),
        color=color
        )
    
    dict_add(total, ans, in_place = True);
    dict_add(row_total, ans, in_place = True);
    old_pt = pt;
  if row_total != {}:
    print("Total for row (" + mixed_frac(pts[-1][0]) + "):")
    print(shorten(row_total));
    print("---------------------------------------");
  print();
  print(num_dots, "dots found")
  print("Total time:", time.time() - t0);
  print("Total:")
  print(shorten(total));
  
  if verbosity >= 0:
    show(graphics);
  try:
    is_symmetric = resolvent_data.get(e) == 2*resolvent_data.get(t) and row is None;
  except NameError:
    is_symmetric = False;
  if is_symmetric: 
    print("\nDiscrepancy:")
    return shorten(dict_add(total, thought_dual(total), operation=lambda x,y:x-y))
  else:
    return total;

def row_plot(row = None, verbosity = 0):
  return zone_plot(row = row, verbosity = verbosity)

def zone_plot_double(verbosity = 0):
  if resolvent_data.get(e) == 2*resolvent_data.get(t):
    print("Using zone_plot() instead");
    return zone_plot(verbosity = verbosity);
  
  total1 = zone_plot(verbosity = verbosity);
  resolvent_data.update({t : resolvent_data.get(e) - resolvent_data.get(t)});
  if verbosity > 0:
    print("============================");
  try:
    total2 = zone_plot(verbosity = verbosity);
  except:
    resolvent_data.update({t : resolvent_data.get(e) - resolvent_data.get(t)});
    raise
  resolvent_data.update({t : resolvent_data.get(e) - resolvent_data.get(t)});
  
  discrep = shorten(dict_add(total1, thought_dual(
    deep_subs(total2, {}, -qTH^((2*t-e).subs(resolvent_data))))))
  print("Discrepancy:")
  print(discrep)
  print();
  return discrep;