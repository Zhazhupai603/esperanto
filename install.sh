#!/usr/bin/env sh
# ESPERANTO one-line installer.
#
# Usage (download + install, no sudo):
#   curl -fsSL https://github.com/Zhazhupai603/esperanto/releases/latest/download/install.sh | sh
#
# Options (environment variables):
#   ESPERANTO_PREFIX  install prefix for the binary (default ~/.local)
#   ESPERANTO_TARBALL local tarball path (skips download when set)
set -eu

VERSION="v1.0.5"
ARCH="linux-x86_64"
PKG_BASE="esperanto-1.0.5-${ARCH}"

# CPU build selection: prefer the AVX2-optimized (-v3) build when the local
# CPU supports it, fall back to the portable baseline build.
if grep -qw avx2 /proc/cpuinfo 2>/dev/null; then
    CPUV="v3"
else
    CPUV=""
fi

# Machines with a working NVIDIA driver prefer the GPU-enabled variant of the
# selected CPU build; each missing package falls back to the next candidate.
GPUV=""
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
    GPUV="-gpu"
fi

suffix_for() {
    # $1 = cpu variant (v3|""), $2 = gpu variant (-gpu|"")
    if [ "$1" = "v3" ] && [ "$2" = "-gpu" ]; then
        echo "-v3-gpu"
    elif [ "$1" = "v3" ]; then
        echo "-v3"
    else
        echo "$2"
    fi
}

CANDIDATES="$(suffix_for "$CPUV" "$GPUV") $(suffix_for "$CPUV" "") $(suffix_for "" "$GPUV") $(suffix_for "" "")"

BASE_URL="https://github.com/Zhazhupai603/esperanto/releases/download/${VERSION}"
PREFIX="${ESPERANTO_PREFIX:-$HOME/.local}"
DATA_DIR="${XDG_DATA_HOME:-$HOME/.local/share}"
BIN_DIR="$PREFIX/bin"
BUNDLE_DIR="$DATA_DIR/esperanto/bundle"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

echo "==> ESPERANTO installer (version $VERSION)"

dl() {
    if command -v curl >/dev/null 2>&1; then
        curl -fsSL "$1" -o "$2"
    elif command -v wget >/dev/null 2>&1; then
        wget -q "$1" -O "$2"
    else
        echo "error: need curl or wget" >&2
        exit 1
    fi
}

if [ -n "${ESPERANTO_TARBALL:-}" ]; then
    echo "==> using local tarball: $ESPERANTO_TARBALL"
    cp "$ESPERANTO_TARBALL" "$tmp/pkg.tar.gz"
    # Local tarballs are extracted by their top-level directory name; detect it.
    PKG="$(tar -tzf "$tmp/pkg.tar.gz" | head -1 | cut -d/ -f1)"
else
    PKG=""
    for suf in $CANDIDATES; do
        cand="${PKG_BASE}${suf}"
        echo "==> downloading $BASE_URL/${cand}.tar.gz"
        if dl "$BASE_URL/${cand}.tar.gz" "$tmp/pkg.tar.gz"; then
            PKG="$cand"
            break
        fi
        echo "==> ${cand} unavailable; trying the next build"
    done
    if [ -z "$PKG" ]; then
        echo "error: no downloadable build found for $VERSION" >&2
        exit 1
    fi
fi

echo "==> extracting $PKG"
tar -xzf "$tmp/pkg.tar.gz" -C "$tmp"

src="$tmp/$PKG"
mkdir -p "$BIN_DIR" "$BUNDLE_DIR"

echo "==> installing binary -> $BIN_DIR/esperanto"
cp "$src/bin/esperanto" "$BIN_DIR/esperanto"
chmod +x "$BIN_DIR/esperanto"

echo "==> installing model bundle -> $BUNDLE_DIR"
cp -r "$src/bundle/." "$BUNDLE_DIR/"

echo "==> verifying"
"$BIN_DIR/esperanto" --version

case ":$PATH:" in
    *":$BIN_DIR:"*) ;;
    *)
        echo "warning: $BIN_DIR is not on your PATH." >&2
        echo "         add it with:  export PATH=\"$BIN_DIR:\$PATH\"" >&2
        ;;
esac

echo "==> done. run:  esperanto --help"
