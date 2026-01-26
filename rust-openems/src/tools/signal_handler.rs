//! Signal handling for graceful shutdown.
//!
//! Provides cross-platform handling of SIGINT (Ctrl-C) and other signals
//! to allow simulations to terminate gracefully.

use std::sync::atomic::{AtomicBool, AtomicU32, Ordering};
use std::sync::Arc;

/// Global signal state.
static SIGINT_RECEIVED: AtomicBool = AtomicBool::new(false);
static SIGINT_COUNT: AtomicU32 = AtomicU32::new(0);
static HANDLER_INSTALLED: AtomicBool = AtomicBool::new(false);

/// Signal handler for managing simulation interrupts.
pub struct SignalHandler {
    /// Whether to force exit on second SIGINT
    force_exit_on_second: bool,
}

impl SignalHandler {
    /// Create a new signal handler.
    pub fn new() -> Self {
        Self {
            force_exit_on_second: true,
        }
    }

    /// Set whether to force exit on second SIGINT.
    pub fn with_force_exit(mut self, force: bool) -> Self {
        self.force_exit_on_second = force;
        self
    }

    /// Install the signal handler.
    pub fn install(&self) -> Result<(), String> {
        if HANDLER_INSTALLED.swap(true, Ordering::SeqCst) {
            return Ok(()); // Already installed
        }

        #[cfg(unix)]
        {
            self.install_unix()
        }

        #[cfg(windows)]
        {
            self.install_windows()
        }

        #[cfg(not(any(unix, windows)))]
        {
            Ok(()) // No-op on other platforms
        }
    }

    #[cfg(unix)]
    fn install_unix(&self) -> Result<(), String> {
        use std::io::Write;

        // Use signal-hook crate pattern if available, otherwise manual handling
        unsafe {
            let result = libc::signal(libc::SIGINT, sigint_handler as usize);
            if result == libc::SIG_ERR {
                return Err("Failed to install SIGINT handler".to_string());
            }
        }

        Ok(())
    }

    #[cfg(windows)]
    fn install_windows(&self) -> Result<(), String> {
        use std::os::raw::c_int;

        extern "system" {
            fn SetConsoleCtrlHandler(
                handler: extern "system" fn(c_int) -> c_int,
                add: c_int,
            ) -> c_int;
        }

        extern "system" fn console_handler(ctrl_type: c_int) -> c_int {
            const CTRL_C_EVENT: c_int = 0;
            const CTRL_BREAK_EVENT: c_int = 1;

            if ctrl_type == CTRL_C_EVENT || ctrl_type == CTRL_BREAK_EVENT {
                handle_sigint();
                return 1; // Handled
            }
            0 // Not handled
        }

        let result = unsafe { SetConsoleCtrlHandler(console_handler, 1) };
        if result == 0 {
            return Err("Failed to install console control handler".to_string());
        }

        Ok(())
    }

    /// Check if SIGINT was received.
    pub fn received_sigint() -> bool {
        SIGINT_RECEIVED.load(Ordering::SeqCst)
    }

    /// Get the number of times SIGINT was received.
    pub fn sigint_count() -> u32 {
        SIGINT_COUNT.load(Ordering::SeqCst)
    }

    /// Reset the signal state.
    pub fn reset() {
        SIGINT_RECEIVED.store(false, Ordering::SeqCst);
        SIGINT_COUNT.store(0, Ordering::SeqCst);
    }

    /// Create an atomic flag that can be shared with the signal handler.
    pub fn create_stop_flag() -> Arc<AtomicBool> {
        Arc::new(AtomicBool::new(false))
    }
}

impl Default for SignalHandler {
    fn default() -> Self {
        Self::new()
    }
}

/// Handle SIGINT signal.
fn handle_sigint() {
    SIGINT_RECEIVED.store(true, Ordering::SeqCst);
    let count = SIGINT_COUNT.fetch_add(1, Ordering::SeqCst) + 1;

    // Write directly to stderr (signal-safe)
    let msg = if count == 1 {
        "\n[openEMS] Received interrupt signal, stopping gracefully...\n"
    } else {
        "\n[openEMS] Received second interrupt, forcing exit...\n"
    };

    // Use write syscall directly (signal-safe)
    #[cfg(unix)]
    unsafe {
        libc::write(2, msg.as_ptr() as *const libc::c_void, msg.len());
    }

    #[cfg(windows)]
    {
        use std::io::Write;
        let _ = std::io::stderr().write_all(msg.as_bytes());
    }

    // Force exit on second signal
    if count >= 2 {
        std::process::exit(1);
    }
}

#[cfg(unix)]
extern "C" fn sigint_handler(_signum: libc::c_int) {
    handle_sigint();
}

/// Guard that checks for interrupts during simulation.
pub struct InterruptGuard {
    /// Stop flag
    stop_flag: Arc<AtomicBool>,
    /// Check interval (timesteps)
    check_interval: u64,
    /// Last check timestep
    last_check: u64,
}

impl InterruptGuard {
    /// Create a new interrupt guard.
    pub fn new(check_interval: u64) -> Self {
        Self {
            stop_flag: Arc::new(AtomicBool::new(false)),
            check_interval: check_interval.max(1),
            last_check: 0,
        }
    }

    /// Check if simulation should stop.
    pub fn should_stop(&mut self, timestep: u64) -> bool {
        // Only check periodically for performance
        if timestep - self.last_check < self.check_interval {
            return self.stop_flag.load(Ordering::Relaxed);
        }

        self.last_check = timestep;

        // Check global signal state
        if SignalHandler::received_sigint() {
            self.stop_flag.store(true, Ordering::SeqCst);
            return true;
        }

        self.stop_flag.load(Ordering::Relaxed)
    }

    /// Manually request stop.
    pub fn request_stop(&self) {
        self.stop_flag.store(true, Ordering::SeqCst);
    }

    /// Get a clone of the stop flag for sharing.
    pub fn stop_flag(&self) -> Arc<AtomicBool> {
        Arc::clone(&self.stop_flag)
    }

    /// Reset the guard.
    pub fn reset(&mut self) {
        self.stop_flag.store(false, Ordering::SeqCst);
        self.last_check = 0;
    }
}

/// RAII guard for signal handler installation.
pub struct SignalGuard {
    _private: (),
}

impl SignalGuard {
    /// Install signal handlers and return a guard.
    pub fn install() -> Result<Self, String> {
        let handler = SignalHandler::new();
        handler.install()?;
        Ok(Self { _private: () })
    }
}

impl Drop for SignalGuard {
    fn drop(&mut self) {
        // Reset signal state when guard is dropped
        SignalHandler::reset();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_signal_handler_creation() {
        let handler = SignalHandler::new();
        assert!(handler.force_exit_on_second);
    }

    #[test]
    fn test_signal_state() {
        // Reset first
        SignalHandler::reset();

        assert!(!SignalHandler::received_sigint());
        assert_eq!(SignalHandler::sigint_count(), 0);
    }

    #[test]
    fn test_interrupt_guard() {
        let mut guard = InterruptGuard::new(10);

        // Should not stop initially
        assert!(!guard.should_stop(0));
        assert!(!guard.should_stop(5));

        // Request stop
        guard.request_stop();
        assert!(guard.should_stop(20));

        // Reset
        guard.reset();
        assert!(!guard.should_stop(0));
    }

    #[test]
    fn test_stop_flag_sharing() {
        let guard = InterruptGuard::new(1);
        let flag = guard.stop_flag();

        assert!(!flag.load(Ordering::Relaxed));

        guard.request_stop();
        assert!(flag.load(Ordering::Relaxed));
    }
}
