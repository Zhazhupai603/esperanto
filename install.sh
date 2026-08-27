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

VERSION="v1.0.2"
ARCH="linux-x86_64"

# Machines with a working NVIDIA driver get the GPU-enabled build.
VARIANT=""
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
    VARIANT="-gpu"
fi

PKG="esperanto-1.0.2-${ARCH}${VARIANT}"
TARBALL="${PKG}.tar.gz"
BASE_URL="https://github.com/Zhazhupai603/esperanto/releases/download/${VERSION}"

PREFIX="${ESPERANTO_PREFIX:-$HOME/.local}"
DATA_DIR="${XDG_DATA_HOME:-$HOME/.local/share}"
BIN_DIR="$PREFIX/bin"
BUNDLE_DIR="$DATA_DIR/esperanto/bundle"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

echo "==> ESPERANTO installer (version $VERSION)"

if [ -n "${ESPERANTO_TARBALL:-}" ]; then
    echo "==> using local tarball: $ESPERANTO_TARBALL"
    cp "$ESPERANTO_TARBALL" "$tmp/$TARBALL"
else
    echo "==> downloading $BASE_URL/$TARBALL"
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
    if ! dl "$BASE_URL/$TARBALL" "$tmp/$TARBALL"; then
        if [ -n "$VARIANT" ]; then
            echo "==> GPU build unavailable for $VERSION; falling back to the CPU build"
            PKG="esperanto-${VERSION#v}-${ARCH}"
            TARBALL="${PKG}.tar.gz"
            dl "$BASE_URL/$TARBALL" "$tmp/$TARBALL"
        else
            exit 1
        fi
    fi
fi

echo "==> extracting"
tar -xzf "$tmp/$TARBALL" -C "$tmp"

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
