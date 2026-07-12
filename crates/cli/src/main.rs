fn main() {
    // Subcommands (run/resume/estimate/verify/fixtures/inspect) land in P1/P2.
    println!(
        "classdiam {} — see `classdiam --help` once subcommands are implemented",
        env!("CARGO_PKG_VERSION")
    );
}
