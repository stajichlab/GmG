// Embed the git commit in `gmg --version` so a binary on the HPCC can be traced to its source.
use std::process::Command;

fn main() {
    let hash = std::env::var("GMG_GIT_HASH")
        .ok()
        .filter(|h| !h.is_empty())
        .or_else(|| {
            Command::new("git")
                .args(["rev-parse", "--short=9", "HEAD"])
                .output()
                .ok()
                .filter(|o| o.status.success())
                .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        })
        .unwrap_or_else(|| "unknown".to_string());
    println!("cargo:rustc-env=GMG_GIT={hash}");
    println!("cargo:rerun-if-env-changed=GMG_GIT_HASH");
    // HEAD alone is not enough: `git pull` on a checked-out branch moves the ref file, not HEAD, and left a stale hash in --version.
    for p in ["HEAD", "refs", "packed-refs"] {
        println!("cargo:rerun-if-changed=../../.git/{p}");
    }
}
