use wasm_bindgen::prelude::*;
use js_sys::Float32Array;
use js_sys::Array;
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}

macro_rules! console_log {
    // Note that this is using the `log` function imported above during
    // `bare_bones`
    ($($t:tt)*) => (log(&format_args!($($t)*).to_string()));
}
#[derive(Debug)]
#[derive(Copy, Clone)]
pub struct XY_global {
    x: f32,
    y: f32,
}
pub struct NormalizedBox {
    b_L: XY_global,
    b_R: XY_global,
    t_L: XY_global,
    t_R: XY_global,
}
pub struct BoxSides {
    b: (XY_global, XY_global),
    l: (XY_global, XY_global),
    t: (XY_global, XY_global),
    r: (XY_global, XY_global),
}
#[derive(Debug, Clone)]
pub struct Polygon {
    center: XY_global,
    end: XY_global,
    n_sides: u8,
    side_points: Vec<XY_global>,
    side_lines: Vec<(XY_global, XY_global)>,
}

impl Polygon {
    fn new(center_x: f32, center_y: f32, end_x: f32, end_y: f32, n_sides: u8) -> Self {
        let mut side_points: Vec<XY_global> = Vec::new();
        let mut side_lines: Vec<(XY_global, XY_global)> = Vec::new();
        let y = (center_y as f64) - (end_y as f64);
        let x = (center_x as f64) - (end_x as f64);

        // Constants, rotation_in_rads can be reassigned once
        let mut rotation_in_rads: f64 = y.atan2(x);

        let polygon_radius: f32 = get_distance_4p(center_x, center_y, end_x, end_y);
        let rad_interval: f64 = (2.0 * std::f64::consts::PI) / (n_sides as f64);

        if n_sides % 2 != 0 {
            rotation_in_rads += ((std::f64::consts::PI * 2.0) / (n_sides as f64)) * 0.5;
        }

        for i in 1..=n_sides {
            let next_rad: f64 = rad_interval * (i as f64) + rotation_in_rads;
            let out: XY_global = XY_global {
                x: center_x + polygon_radius * (next_rad.cos() as f32),
                y: center_y + polygon_radius * (next_rad.sin() as f32),
            };
            side_points.push(out);
        }

        for i in 0..n_sides {
            if i == n_sides - 1 {
                // Last iteration
                let first_out: XY_global = *side_points.get(i as usize).unwrap();

                let second_out: XY_global = *side_points.get(0).unwrap();

                side_lines.push((first_out, second_out));
            } else {
                let first_out: XY_global = *side_points.get(i as usize).unwrap();

                let second_out: XY_global = *side_points.get((i as usize) + (1 as usize)).unwrap();

                side_lines.push((first_out, second_out));
            }
        }
        Self {
            center: XY_global { x: center_x, y: center_y },
            end: XY_global { x: end_x, y: end_y },
            n_sides: n_sides,
            side_points: side_points,
            side_lines: side_lines,
        }
    }
}

pub struct CircleQuadrants {
    center: XY_global,
    top: XY_global,
    right: XY_global,
    left: XY_global,
    bottom: XY_global,
}

impl BoxSides {
    fn new(s: &NormalizedBox) -> Self {
        Self {
            b: (s.b_L, s.b_R),
            l: (s.b_L, s.t_L),
            r: (s.b_R, s.t_R),
            t: (s.t_L, s.t_R),
        }
    }
}
impl NormalizedBox {
    fn new(sx: f32, sy: f32, ex: f32, ey: f32) -> NormalizedBox<> {
        let mut pos_x: bool = false;
        let mut pos_y: bool = false;

        if sx < ex {
            pos_x = true;
        }
        if sy < ey {
            pos_y = true;
        }

        let b_L: XY_global = XY_global {
            x: if pos_x {
                sx
            } else {
                ex
            },
            y: if pos_y {
                sy
            } else {
                ey
            },
        };
        let b_R: XY_global = XY_global {
            x: if pos_x {
                ex
            } else {
                sx
            },
            y: if pos_y {
                sy
            } else {
                ey
            },
        };
        let t_L: XY_global = XY_global {
            x: if pos_x {
                sx
            } else {
                ex
            },
            y: if pos_y {
                ey
            } else {
                sy
            },
        };
        let t_R: XY_global = XY_global {
            x: if pos_x {
                ex
            } else {
                sx
            },
            y: if pos_y {
                ey
            } else {
                sy
            },
        };

        NormalizedBox { b_L, b_R, t_L, t_R }
    }
}

#[wasm_bindgen]
#[derive(Clone)]
pub struct Geometry {
    sx: f32,
    sy: f32,
    ex: f32,
    ey: f32,
}
#[wasm_bindgen]
pub struct RealGeometry {
    geometry: Vec<Geometry>,
}

#[wasm_bindgen]
pub fn derive_actual_pos(cursor: f32, stage: f32, inverse_z: f32) -> f32 {
    let pos: f32 = (cursor + stage) * inverse_z;
    return (pos * 1000.0).round() / 1000.0;
}
#[wasm_bindgen] //Get distance based on starting x,y and ending x,y
pub fn get_distance_4p(sx: f32, sy: f32, ex: f32, ey: f32) -> f32 {
    let delta_x: f32 = sx - ex;
    let delta_y: f32 = sy - ey;
    let sqrt_x: f32 = delta_x.powf(2.0);
    let sqrt_y: f32 = delta_y.powf(2.0);
    (sqrt_x + sqrt_y).sqrt()
}

#[wasm_bindgen]
pub fn return_jsarr(arr: &Float32Array) -> Float32Array {
    // Best to allocate memory first, then alter it in place
    let r: Vec<f32> = arr.to_vec();
    let out: Vec<f32> = r
        .into_iter()
        .map(|x| x * 2.0)
        .collect::<Vec<_>>();

    Float32Array::from(&out[..])
}

#[wasm_bindgen]
pub fn check_geo_collision(sel_box: &Float32Array, real_geo: &Array) {
    // First lets normalize the selection box
    // Allocating
    let mut sel_arr: [f32; 4] = [0.0; 4];
    let geo: Vec<Geometry> =
        vec![Geometry {sx: 0.0, sy: 0.0, ex: 0.0, ey: 0.0}; real_geo.length() as usize];

    // Copying
    Float32Array::copy_to(sel_box, &mut sel_arr);
    real_geo.for_each(
        &mut (|v, _, _| {
            let f: Option<f64> = v.as_f64();
            if f.is_some() {
            }
        })
    );
}

fn check_line_intersect(f_line: (XY_global, XY_global), s_line: (XY_global, XY_global)) -> bool {
    let (s, e) = f_line;
    let (fs, fe) = s_line;
    let XY_global { x: a, y: b } = s;
    let XY_global { x: c, y: d } = e;
    let XY_global { x: p, y: q } = fs;
    let XY_global { x: r, y: s } = fe;

    let det: f32 = (c - a) * (s - q) - (r - p) * (d - b);
    if det.abs() < f32::EPSILON {
        return false;
    } else {
        let lambda: f32 = ((s - q) * (r - a) + (p - r) * (s - b)) / det;
        let gamma: f32 = ((b - d) * (r - a) + (c - a) * (s - b)) / det;
        if 0.0 < lambda && lambda < 1.0 && 0.0 < gamma && gamma < 1.0 {
            return true;
        }
    }
    return false;
}

#[wasm_bindgen]
pub fn check_line_collision(
    bsx: f32,
    bsy: f32,
    bex: f32,
    bey: f32,
    gsx: f32,
    gsy: f32,
    gex: f32,
    gey: f32
) -> bool {
    let n_box = NormalizedBox::new(bsx, bsy, bex, bey);
    let NormalizedBox { b_L, b_R, t_L, t_R } = n_box;
    /*

if gsx >= b_L.x && gsx <= b_R.x && gsy >= b_L.y && gsy <= t_L.y {
    return true;
}
*/
    if is_between(b_L.x, b_R.x, gsx) && is_between(b_L.y, t_L.y, gsy) {
        return true;
    }

    if is_between(b_L.x, b_R.x, gex) && is_between(b_L.y, t_L.y, gey) {
        return true;
    }

    let box_sides: BoxSides = BoxSides::new(&n_box);
    let b_side_arr: [(XY_global, XY_global); 4] = [
        box_sides.b,
        box_sides.r,
        box_sides.l,
        box_sides.t,
    ];
    let line_to_check: (XY_global, XY_global) = (
        XY_global { x: gsx, y: gsy },
        XY_global { x: gex, y: gey },
    );

    for side in b_side_arr {
        if check_line_intersect(side, line_to_check) {
            return true;
        }
    }

    false
}
#[wasm_bindgen]
pub fn check_rect_collision(
    bsx: f32,
    bsy: f32,
    bex: f32,
    bey: f32,
    gsx: f32,
    gsy: f32,
    gex: f32,
    gey: f32
) -> bool {
    let s_box = NormalizedBox::new(bsx, bsy, bex, bey);
    let i_box = NormalizedBox::new(gsx, gsy, gex, gey);
    let NormalizedBox { b_L, b_R, t_L, t_R } = s_box;
    let NormalizedBox { b_L: i_b_L, b_R: i_b_R, t_L: i_t_L, t_R: i_t_R } = i_box;

    let i_box_arr: [XY_global; 4] = [i_b_L, i_b_R, i_t_L, i_t_R];

    for point in i_box_arr {
        if is_between(b_L.x, b_R.x, point.x) && is_between(b_L.y, t_L.y, point.y) {
            return true;
        }
    }

    let s_box_sides = BoxSides::new(&s_box);
    let i_box_sides = BoxSides::new(&i_box);

    let s_side_arr: [(XY_global, XY_global); 4] = [
        s_box_sides.b,
        s_box_sides.r,
        s_box_sides.l,
        s_box_sides.t,
    ];
    let i_side_arr: [(XY_global, XY_global); 4] = [
        i_box_sides.b,
        i_box_sides.r,
        i_box_sides.l,
        i_box_sides.t,
    ];

    for s_side in s_side_arr {
        for i_side in i_side_arr {
            if check_line_intersect(i_side, s_side) {
                return true;
            }
        }
    }

    return false;
}

#[wasm_bindgen]
pub fn check_circle_collision(
    bsx: f32,
    bsy: f32,
    bex: f32,
    bey: f32,
    gsx: f32,
    gsy: f32,
    gex: f32,
    gey: f32
) -> bool {
    // Can possibly edit to exclude circle centers triggering
    let n_box = NormalizedBox::new(bsx, bsy, bex, bey);
    let NormalizedBox { b_L, b_R, t_L, t_R } = n_box;
    let n_box_arr: [XY_global; 4] = [b_L, b_R, t_L, t_R];
    let circle_radius: f32 = get_distance_4p(gsx, gsy, gex, gey);
    if
        is_between(b_L.x, b_R.x, gsx) &&
        is_between(b_L.y, t_L.y, gsy) &&
        (circle_radius < (b_L.x - b_R.x).abs() || circle_radius < (b_L.y - t_R.y).abs())
    {
        // Possibly check for selection box width / height to be greater than radius to ensure box crosses circle circumference
        return true;
    }

    for box_point in n_box_arr {
        let distance_to_center = get_distance_4p(gsx, gsy, box_point.x, box_point.y);
        if
            distance_to_center < circle_radius &&
            (circle_radius < (b_L.x - b_R.x).abs() || circle_radius < (b_L.y - t_R.y).abs())
        {
            return true;
        }
    }

    let box_sides: BoxSides = BoxSides::new(&n_box);
    let s_side_arr: [(XY_global, XY_global); 4] = [
        box_sides.b,
        box_sides.r,
        box_sides.l,
        box_sides.t,
    ];
    let circle_quads: CircleQuadrants = CircleQuadrants {
        center: XY_global { x: gsx, y: gsy },
        top: XY_global { x: gsx, y: gsy + circle_radius },
        bottom: XY_global { x: gsx, y: gsy - circle_radius },
        right: XY_global { x: gsx + circle_radius, y: gsy },
        left: XY_global { x: gsx - circle_radius, y: gsy },
    };
    let circle_quad_arr: [XY_global; 4] = [
        circle_quads.top,
        circle_quads.bottom,
        circle_quads.left,
        circle_quads.right,
    ];

    for box_side in s_side_arr {
        for quad_line in circle_quad_arr {
            if check_line_intersect((circle_quads.center, quad_line), box_side) {
                return true;
            }
        }
    }

    return false;
}

#[wasm_bindgen]
pub fn check_polygon_collision(
    bsx: f32,
    bsy: f32,
    bex: f32,
    bey: f32,
    gsx: f32,
    gsy: f32,
    gex: f32,
    gey: f32,
    sides: u8
) -> bool {
    let polygon: Polygon = Polygon::new(gsx, gsy, gex, gey, sides);

    let n_box: NormalizedBox = NormalizedBox::new(bsx, bsy, bex, bey);
    let NormalizedBox { b_L, b_R, t_L, t_R } = n_box;

    if
        is_between(b_L.x, b_R.x, gsx) &&
        is_between(b_L.y, t_L.y, gsy) &&
        (get_distance_4p(gsx, gsy, gex, gey) < (b_L.x - b_R.x).abs() ||
            get_distance_4p(gsx, gsy, gex, gey) < (b_L.y - t_R.y).abs())
    {
        return true;
    }

    let box_sides: BoxSides = BoxSides::new(&n_box);

    let s_side_arr: [(XY_global, XY_global); 4] = [
        box_sides.b,
        box_sides.r,
        box_sides.l,
        box_sides.t,
    ];

    for box_side in s_side_arr {
        for poly_side in &polygon.side_lines[..] {
            if check_line_intersect(*poly_side, box_side) {
                return true;
            }
        }
    }
    return false;
}
#[wasm_bindgen]
pub fn find_circle_tan_points(
    circle_center_x: f32,
    circle_center_y: f32,
    circle_radius: f32,
    point_x: f32,
    point_y: f32
) -> Array {
    let out: Array = Array::new();
    let center_point: XY_global = XY_global {
        x: (circle_center_x + point_x) / 2.0,
        y: (circle_center_y + point_y) / 2.0,
    };
    let dx: f32 = center_point.x - circle_center_x;
    let dy: f32 = center_point.y - circle_center_y;
    let distance_between_centers = get_distance_4p(
        circle_center_x,
        circle_center_y,
        center_point.x,
        center_point.y
    );

    let a: f32 =
        circle_radius.powf(2.0) -
        distance_between_centers.powf(2.0) +
        distance_between_centers.powf(2.0);
    let a: f32 = a / (2.0 * distance_between_centers);

    let x2: f32 = circle_center_x + dx * (a / distance_between_centers);
    let y2: f32 = circle_center_y + dy * (a / distance_between_centers);

    let h: f32 = (circle_radius.powf(2.0) - a.powf(2.0)).abs().sqrt();

    let rx: f32 = -dy * (h / distance_between_centers);
    let ry: f32 = dx * (h / distance_between_centers);

    let x0 = (x2 + rx) as f64;
    let y0 = (y2 + ry) as f64;
    let x1 = (x2 - rx) as f64;
    let y1 = (y2 - ry) as f64;
    console_log!("a: {:?}", a);
    console_log!("x2: {:?}", x2);
    console_log!("y2: {:?}", y2);
    console_log!("h: {:?}", h);
    console_log!("rx: {:?}", rx);
    console_log!("ry: {:?}", ry);
    out.push(&JsValue::from_f64(x0));
    out.push(&JsValue::from_f64(y0));
    out.push(&JsValue::from_f64(x1));
    out.push(&JsValue::from_f64(y1));
    return out;
}

fn is_between(n1: f32, n2: f32, between: f32) -> bool {
    if between >= n1 && between <= n2 {
        return true;
    }
    return false;
}
