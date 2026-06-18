pub mod avx2;
pub mod neon;
pub mod scalar;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimdBackend {
    Avx2,
    Neon,
    Scalar,
}

impl SimdBackend {
    pub fn as_str(self) -> &'static str {
        match self {
            SimdBackend::Avx2 => "avx2",
            SimdBackend::Neon => "neon",
            SimdBackend::Scalar => "scalar",
        }
    }
}

pub fn active_backend() -> SimdBackend {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") {
            SimdBackend::Avx2
        } else {
            SimdBackend::Scalar
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        SimdBackend::Neon
    }

    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        SimdBackend::Scalar
    }
}

pub fn hash_bytes32(data: &[u8]) -> u64 {
    match active_backend() {
        SimdBackend::Avx2 => avx2::hash_bytes32(data),
        SimdBackend::Neon => neon::hash_bytes32(data),
        SimdBackend::Scalar => scalar::hash_bytes32(data),
    }
}
