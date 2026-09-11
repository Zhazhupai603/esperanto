fn main() {
    let dev = candle_core::Device::new_cuda(0).expect("cuda device");
    let a = candle_core::Tensor::randn(0f32, 1f32, (512, 512), &dev).unwrap();
    let b = candle_core::Tensor::randn(0f32, 1f32, (512, 512), &dev).unwrap();
    let c = a.matmul(&b).unwrap();
    let s = c.sum_all().unwrap().to_scalar::<f32>().unwrap();
    println!("gpu matmul ok, sum={s:.3}");
}
