

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
    ($($t:tt)*) => (log(&format_args!($($t)*).to_string()))
}

#[derive(Copy, Clone)]
pub struct XY_global {
    x:  f32,
    y:  f32
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
    r: (XY_global, XY_global)
}
impl BoxSides {
    fn new (s: &NormalizedBox ) -> Self {
        Self {
            b: (s.b_L, s.b_R),
            l: (s.b_L, s.t_R),
            r: (s.b_R, s.t_R),
            t: (s.t_L, s.t_R)
        }
    }
}
impl NormalizedBox {
    fn new(sx:  f32, sy:  f32, ex:  f32, ey:  f32) -> NormalizedBox<> {
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
            }, y: if pos_y {
                sy
            } else {
                ey
            }
        };
        let b_R: XY_global = XY_global { 
            x: if pos_x {
                ex
            } else { 
                sx 
            }, y: if pos_y {
                sy
            } else {
                ey
            }
        };
        let t_L: XY_global = XY_global { 
            x: if pos_x {
                sx 
            } else { 
                ex
            }, y: if pos_y {
                ey
            } else {
                sy
            }
        };
        let t_R: XY_global = XY_global { 
            x: if pos_x {
                ex
            } else { 
                ex 
            }, y: if pos_y {
                ey
            } else {
                sy
            }
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
    geometry: Vec<Geometry>
}



#[wasm_bindgen]
pub fn derive_actual_pos(cursor: f32, stage: f32, inverse_z: f32 ) -> f32 {
    let pos: f32 = (cursor + stage) * inverse_z;
    return (pos * 1000.0).round() / 1000.0;
}
#[wasm_bindgen] //Get distance based on starting x,y and ending x,y
pub fn get_distance_4p(sx: f32, sy: f32, ex:f32, ey: f32) -> f32 {
    let delta_x: f32 = sx - ex;
    let delta_y: f32 = sy - ey;
    let sqrt_x: f32 = delta_x.powf(2.0);
    let sqrt_y: f32 = delta_y.powf(2.0);
    let _c = Float32Array::new_with_length(10);
    (sqrt_x + sqrt_y).sqrt()
}

#[wasm_bindgen]
pub fn return_jsarr(arr: &Float32Array) -> Float32Array {
    // Best to allocate memory first, then alter it in place
    let r: Vec<f32>= arr.to_vec();
    let out: Vec<f32> = r.into_iter()
        .map(|x| x * 2.0)
        .collect::<Vec<_>>();
   
   Float32Array::from(&out[..])
}

#[wasm_bindgen] 
pub fn check_geo_collision(sel_box: &Float32Array, real_geo: &Array) {
    // First lets normalize the selection box
    // Allocating
    let mut sel_arr: [f32; 4] = [ 0.0; 4 ];
    let geo: Vec<Geometry> = vec![Geometry {sx: 0.0, sy: 0.0, ex: 0.0, ey: 0.0}; real_geo.length() as usize];

    // Copying
    Float32Array::copy_to(sel_box, &mut sel_arr);
    real_geo.for_each(&mut |v, _, _| {
        let f: Option<f64> = v.as_f64();
        if f.is_some() {
            
        }
        
    });
    
}

fn check_line_intersect(f_line: (XY_global, XY_global), s_line: (XY_global, XY_global)) -> bool {

    let (s, e) = f_line;
    let (fs, fe) = s_line;
    let XY_global {x: a, y: b} = s;
    let XY_global {x: c, y: d} = e;
    let XY_global {x: p, y: q} = fs;
    let XY_global {x: r, y: s} = fe;


    let det: f32 = (c - a) * (s - q) - (r - p) * (d - b);
    if det.abs() < f32::EPSILON {
        return false
    } else {
        let lambda: f32 = ((s - q) * (r - a) + (p - r) * (s - b)) / det;
        let gamma: f32 = ((b - d) * (r - a) + (c - a) * (s - b)) / det;
        if 0.0 < lambda && lambda < 1.0 && 0.0 < gamma && gamma < 1.0 { return true; }
    }
    return false
}


#[wasm_bindgen]
pub fn check_line_collision(bsx: f32, bsy: f32, bex: f32, bey: f32, gsx: f32, gsy: f32, gex: f32, gey: f32,) -> bool {
    let n_box = NormalizedBox::new(bsx, bsy, bex, bey);
    let NormalizedBox {
        b_L, 
        b_R, 
        t_L,
        t_R, 
    } = n_box;

    if gsx >= b_L.x && gsx <= b_R.x && gsy >= b_L.y && gsy <= t_L.y { return true; } 
    
    if gex >= b_L.x && gex <= b_R.x && gey >= b_L.y && gey <= t_L.y { return true; } 
    
    let box_sides: BoxSides = BoxSides::new(&n_box);
    let b_side_arr: [(XY_global, XY_global); 4] = [box_sides.b, box_sides.r, box_sides.l, box_sides.t];
    let line_to_check = ( XY_global{x: gsx, y: gsy }, XY_global{ x: gex, y: gey } );
    
    for side in b_side_arr {
        if check_line_intersect(side, line_to_check) { return true; }
    }

    false
}
