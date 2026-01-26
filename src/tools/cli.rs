//! Global CLI options and configuration.
//!
//! Centralized CLI option handling for openEMS simulations.

use std::path::PathBuf;

/// Global simulation options parsed from CLI or configuration.
#[derive(Debug, Clone)]
pub struct GlobalOptions {
    /// Number of threads (0 = auto)
    pub num_threads: usize,
    /// Verbose level (0-3)
    pub verbose: u8,
    /// Engine type string
    pub engine: String,
    /// Debug mode
    pub debug: bool,
    /// Disable all field dumps
    pub disable_dumps: bool,
    /// Simulation directory
    pub sim_path: Option<PathBuf>,
    /// Output directory
    pub output_path: Option<PathBuf>,
}

impl Default for GlobalOptions {
    fn default() -> Self {
        Self {
            num_threads: 0,
            verbose: 1,
            engine: "parallel".to_string(),
            debug: false,
            disable_dumps: false,
            sim_path: None,
            output_path: None,
        }
    }
}

impl GlobalOptions {
    /// Create new global options.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create options with specified thread count.
    pub fn with_threads(mut self, num_threads: usize) -> Self {
        self.num_threads = num_threads;
        self
    }

    /// Set verbose level.
    pub fn with_verbose(mut self, verbose: u8) -> Self {
        self.verbose = verbose.min(3);
        self
    }

    /// Set engine type.
    pub fn with_engine(mut self, engine: &str) -> Self {
        self.engine = engine.to_string();
        self
    }

    /// Set debug mode.
    pub fn with_debug(mut self, debug: bool) -> Self {
        self.debug = debug;
        self
    }

    /// Set disable dumps flag.
    pub fn with_disable_dumps(mut self, disable: bool) -> Self {
        self.disable_dumps = disable;
        self
    }

    /// Set simulation path.
    pub fn with_sim_path(mut self, path: PathBuf) -> Self {
        self.sim_path = Some(path);
        self
    }

    /// Set output path.
    pub fn with_output_path(mut self, path: PathBuf) -> Self {
        self.output_path = Some(path);
        self
    }

    /// Apply thread configuration.
    pub fn apply_thread_config(&self) {
        if self.num_threads > 0 {
            // Note: This should be called before any parallel work
            let _ = rayon::ThreadPoolBuilder::new()
                .num_threads(self.num_threads)
                .build_global();
        }
    }

    /// Check if verbose output should be shown.
    pub fn is_verbose(&self) -> bool {
        self.verbose > 0
    }

    /// Check if debug output should be shown.
    pub fn is_debug(&self) -> bool {
        self.debug || self.verbose >= 3
    }

    /// Get effective output directory.
    pub fn get_output_dir(&self) -> PathBuf {
        if let Some(ref path) = self.output_path {
            path.clone()
        } else if let Some(ref sim_path) = self.sim_path {
            sim_path.clone()
        } else {
            PathBuf::from(".")
        }
    }
}

/// Argument parsing helper for common options.
pub struct CliParser;

impl CliParser {
    /// Parse engine type from string.
    pub fn parse_engine_type(engine: &str) -> crate::fdtd::EngineType {
        match engine.to_lowercase().as_str() {
            "basic" => crate::fdtd::EngineType::Basic,
            "simd" | "sse" => crate::fdtd::EngineType::Simd,
            "parallel" | "multithread" | "mt" => crate::fdtd::EngineType::Parallel,
            _ => crate::fdtd::EngineType::Parallel,
        }
    }

    /// Parse verbose level from string.
    pub fn parse_verbose(verbose_str: &str) -> u8 {
        match verbose_str.to_lowercase().as_str() {
            "quiet" | "q" | "0" => 0,
            "normal" | "n" | "1" => 1,
            "verbose" | "v" | "2" => 2,
            "debug" | "d" | "3" => 3,
            _ => verbose_str.parse().unwrap_or(1),
        }
    }

    /// Get CPU features string for logging.
    pub fn get_cpu_features() -> String {
        let mut features = Vec::new();

        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("sse") {
                features.push("SSE");
            }
            if is_x86_feature_detected!("sse2") {
                features.push("SSE2");
            }
            if is_x86_feature_detected!("sse4.1") {
                features.push("SSE4.1");
            }
            if is_x86_feature_detected!("avx") {
                features.push("AVX");
            }
            if is_x86_feature_detected!("avx2") {
                features.push("AVX2");
            }
            if is_x86_feature_detected!("avx512f") {
                features.push("AVX-512");
            }
            if is_x86_feature_detected!("fma") {
                features.push("FMA");
            }
        }

        #[cfg(target_arch = "aarch64")]
        {
            features.push("NEON");
        }

        if features.is_empty() {
            "None detected".to_string()
        } else {
            features.join(", ")
        }
    }

    /// Print openEMS banner.
    pub fn print_banner() {
        println!(" ---------------------------------------------------------------------- ");
        println!(
            " | openEMS-Rust {} -- High-Performance FDTD Solver",
            crate::VERSION
        );
        println!(" | (c) 2024 openEMS contributors  GPL license");
        println!(" ---------------------------------------------------------------------- ");
        println!(" | CPU features: {}", Self::get_cpu_features());
        println!(" ---------------------------------------------------------------------- ");
    }
}

/// Environment configuration utilities.
pub struct EnvConfig;

impl EnvConfig {
    /// Get number of threads from environment variable.
    pub fn get_num_threads() -> Option<usize> {
        std::env::var("OPENEMS_NUM_THREADS")
            .ok()
            .and_then(|s| s.parse().ok())
    }

    /// Get verbose level from environment variable.
    pub fn get_verbose() -> Option<u8> {
        std::env::var("OPENEMS_VERBOSE")
            .ok()
            .and_then(|s| s.parse().ok())
    }

    /// Check if running in quiet mode.
    pub fn is_quiet() -> bool {
        std::env::var("OPENEMS_QUIET").is_ok()
    }

    /// Get engine type from environment variable.
    pub fn get_engine() -> Option<String> {
        std::env::var("OPENEMS_ENGINE").ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_global_options_default() {
        let opts = GlobalOptions::default();
        assert_eq!(opts.num_threads, 0);
        assert_eq!(opts.verbose, 1);
        assert_eq!(opts.engine, "parallel");
        assert!(!opts.debug);
    }

    #[test]
    fn test_global_options_builder() {
        let opts = GlobalOptions::new()
            .with_threads(4)
            .with_verbose(2)
            .with_engine("simd")
            .with_debug(true);

        assert_eq!(opts.num_threads, 4);
        assert_eq!(opts.verbose, 2);
        assert_eq!(opts.engine, "simd");
        assert!(opts.debug);
    }

    #[test]
    fn test_parse_engine_type() {
        use crate::fdtd::EngineType;

        assert_eq!(CliParser::parse_engine_type("basic"), EngineType::Basic);
        assert_eq!(CliParser::parse_engine_type("simd"), EngineType::Simd);
        assert_eq!(CliParser::parse_engine_type("sse"), EngineType::Simd);
        assert_eq!(
            CliParser::parse_engine_type("parallel"),
            EngineType::Parallel
        );
        assert_eq!(CliParser::parse_engine_type("mt"), EngineType::Parallel);
        assert_eq!(
            CliParser::parse_engine_type("unknown"),
            EngineType::Parallel
        );
    }

    #[test]
    fn test_parse_verbose() {
        assert_eq!(CliParser::parse_verbose("quiet"), 0);
        assert_eq!(CliParser::parse_verbose("0"), 0);
        assert_eq!(CliParser::parse_verbose("normal"), 1);
        assert_eq!(CliParser::parse_verbose("verbose"), 2);
        assert_eq!(CliParser::parse_verbose("debug"), 3);
        assert_eq!(CliParser::parse_verbose("invalid"), 1);
    }

    #[test]
    fn test_is_verbose() {
        let quiet = GlobalOptions::new().with_verbose(0);
        assert!(!quiet.is_verbose());

        let normal = GlobalOptions::new().with_verbose(1);
        assert!(normal.is_verbose());
    }

    #[test]
    fn test_is_debug() {
        let normal = GlobalOptions::new().with_verbose(1);
        assert!(!normal.is_debug());

        let debug_verbose = GlobalOptions::new().with_verbose(3);
        assert!(debug_verbose.is_debug());

        let debug_flag = GlobalOptions::new().with_debug(true);
        assert!(debug_flag.is_debug());
    }

    #[test]
    fn test_get_output_dir() {
        let opts = GlobalOptions::new();
        assert_eq!(opts.get_output_dir(), PathBuf::from("."));

        let with_sim = GlobalOptions::new().with_sim_path(PathBuf::from("/sim"));
        assert_eq!(with_sim.get_output_dir(), PathBuf::from("/sim"));

        let with_out = GlobalOptions::new()
            .with_sim_path(PathBuf::from("/sim"))
            .with_output_path(PathBuf::from("/out"));
        assert_eq!(with_out.get_output_dir(), PathBuf::from("/out"));
    }

    #[test]
    fn test_cpu_features() {
        let features = CliParser::get_cpu_features();
        // Should return some string, either features or "None detected"
        assert!(!features.is_empty());
    }
}
