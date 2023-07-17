
use wasm_bindgen::prelude::*;


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

    (sqrt_x + sqrt_y).sqrt()
}

pub fn check_collision() {}

