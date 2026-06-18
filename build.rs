fn main() {
    println!("cargo:rerun-if-changed=src/gpu/kernels.cu");

    if std::env::var("CARGO_FEATURE_CUDA").is_ok() {
        cc::Build::new()
            .cuda(true)
            .file("src/gpu/kernels.cu")
            .compile("kira_cuda_kernels");
    }
}
