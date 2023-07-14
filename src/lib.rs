
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub fn add(left: usize, right: usize) -> usize {
    left + right
}
#[wasm_bindgen]
pub fn derive_actual_pos(cursor: f32, stage: f32, inverse_z: f32 ) -> f32 {
    (((cursor + stage) * inverse_z)  / 1000.0).round() * 1000.0
}