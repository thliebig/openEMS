#!/bin/bash
# Build openEMS for multiple platforms
# Requires cross-rs: cargo install cross

set -e

VERSION=$(cargo metadata --no-deps --format-version 1 | jq -r '.packages[0].version')
DIST_DIR="dist/v${VERSION}"

echo "Building openEMS v${VERSION} for all platforms..."

mkdir -p "$DIST_DIR"

# Linux x86_64 (static with musl)
echo "Building for Linux x86_64 (static musl)..."
cross build --release --target x86_64-unknown-linux-musl
cp target/x86_64-unknown-linux-musl/release/openems "$DIST_DIR/openems-linux-x86_64"
cp target/x86_64-unknown-linux-musl/release/nf2ff "$DIST_DIR/nf2ff-linux-x86_64"

# macOS x86_64
if [[ "$(uname)" == "Darwin" ]] || command -v cross &> /dev/null; then
    echo "Building for macOS x86_64..."
    cross build --release --target x86_64-apple-darwin 2>/dev/null || true
    if [ -f target/x86_64-apple-darwin/release/openems ]; then
        cp target/x86_64-apple-darwin/release/openems "$DIST_DIR/openems-macos-x86_64"
        cp target/x86_64-apple-darwin/release/nf2ff "$DIST_DIR/nf2ff-macos-x86_64"
    fi
fi

# macOS ARM64 (Apple Silicon)
if [[ "$(uname)" == "Darwin" ]] || command -v cross &> /dev/null; then
    echo "Building for macOS ARM64..."
    cross build --release --target aarch64-apple-darwin 2>/dev/null || true
    if [ -f target/aarch64-apple-darwin/release/openems ]; then
        cp target/aarch64-apple-darwin/release/openems "$DIST_DIR/openems-macos-arm64"
        cp target/aarch64-apple-darwin/release/nf2ff "$DIST_DIR/nf2ff-macos-arm64"
    fi
fi

# Windows x86_64
echo "Building for Windows x86_64..."
cross build --release --target x86_64-pc-windows-gnu 2>/dev/null || true
if [ -f target/x86_64-pc-windows-gnu/release/openems.exe ]; then
    cp target/x86_64-pc-windows-gnu/release/openems.exe "$DIST_DIR/openems-windows-x86_64.exe"
    cp target/x86_64-pc-windows-gnu/release/nf2ff.exe "$DIST_DIR/nf2ff-windows-x86_64.exe"
fi

echo ""
echo "Build complete! Binaries are in $DIST_DIR/"
ls -lh "$DIST_DIR/"
