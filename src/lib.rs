

use wasm_bindgen::prelude::*;
use js_sys::Float32Array;
use js_sys::Array;


pub struct XY_global {
    x: f32,
    y: f32
}
#[wasm_bindgen]
pub struct NormalizedBox {
    b_L: XY_global,
    b_R: XY_global,
    t_L: XY_global,
    t_R: XY_global
}

impl NormalizedBox {
    fn new() {}
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
pub fn check_geo_collision(sel_box: &Float32Array, real_geo: &Array<>) {
    // First lets normalize the selection box
    // Allocating
    let mut sel_arr: [f32; 4] = [ 0.0; 4 ];
    let geo: Vec<f32> = vec![0.0; real_geo.length() as usize];

    // Copying
    Float32Array::copy_to(sel_box, &mut sel_arr[..]);
}

