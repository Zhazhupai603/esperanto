//! EA-Myers bit-parallel verifier.
//!
//! Edit distance under the EA substitution rules ([`crate::ea_free`]):
//! identity and the editing pairs (pattern `G` / text `A`, pattern `C` /
//! text `T`) cost 0, everything else costs 1, indels cost 1. The pattern
//! is the read side, the text is the reference (transcript) side; both
//! are matched case-insensitively.
//!
//! * [`infix`] — pattern must match in full, text overhang is free on
//!   both sides (leading column delta `hin = 0`, running column minimum).
//!   Single 128-bit block, `m <= 128`.
//! * [`long`] — the same semantics for `128 < m <= 256` using two chained
//!   blocks with the horizontal carry passed through `hout`.
//! * [`global`] — full-pattern / full-text distance with forced ends
//!   (leading column delta `+1` per text base); single block for
//!   `m <= 128`, two blocks otherwise.

/// 2-bit code with an `N` marker, for the matcher's Peq tables.
#[inline]
fn code5(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// Infix EA distance: full pattern anchored anywhere inside `text`
/// (overhang free on both sides). Single block, pattern length <= 128.
pub fn infix(read: &[u8], text: &[u8]) -> u32 {
    debug_assert!(
        read.len() <= 128,
        "myers::infix is the single-block entry point; use myers::long for m > 128"
    );
    MyersEa::new(read).infix_min(text)
}

/// Global EA distance between the full pattern and the full text
/// (both ends forced). Single block for `m <= 128`, else two blocks.
pub fn global(read: &[u8], text: &[u8]) -> u32 {
    MyersEa::new(read).global(text)
}

/// Two-block variants for `128 < m <= 256`.
pub mod long {
    use super::MyersEa;

    /// Infix EA distance for patterns of length `128 < m <= 256`.
    pub fn infix(read: &[u8], text: &[u8]) -> u32 {
        debug_assert!(
            read.len() <= 256,
            "myers::long::infix supports patterns up to 256 bases"
        );
        MyersEa::new(read).infix_min(text)
    }

    /// Global EA distance for patterns of length `128 < m <= 256`.
    pub fn global(read: &[u8], text: &[u8]) -> u32 {
        MyersEa::new(read).global(text)
    }
}

/// EA-redefined bit-parallel Myers matcher over patterns up to 256 bases.
#[derive(Debug, Clone)]
struct MyersEa {
    m: usize,
    blocks: usize,
    /// Peq[block][code]: pattern positions (within the block) having code.
    peq: [[u128; 5]; 2],
    init_pv: [u128; 2],
    top: [u128; 2],
}

impl MyersEa {
    /// Build the matcher for `pattern` (length 1..=256; asserts the bound).
    fn new(pattern: &[u8]) -> Self {
        assert!(
            !pattern.is_empty() && pattern.len() <= 256,
            "EA-Myers pattern length must be in 1..=256"
        );
        let m = pattern.len();
        let blocks = m.div_ceil(128);
        let mut s = MyersEa {
            m,
            blocks,
            peq: [[0; 5]; 2],
            init_pv: [0; 2],
            top: [0; 2],
        };
        for (i, &b) in pattern.iter().enumerate() {
            let blk = i / 128;
            let bit = i % 128;
            s.peq[blk][code5(b) as usize] |= 1u128 << bit;
        }
        let b0 = m.min(128);
        s.init_pv[0] = if b0 == 128 { !0u128 } else { (1u128 << b0) - 1 };
        s.top[0] = 1u128 << (b0 - 1);
        if blocks == 2 {
            let b1 = m - 128;
            s.init_pv[1] = if b1 == 128 { !0u128 } else { (1u128 << b1) - 1 };
            s.top[1] = 1u128 << (b1 - 1);
        }
        s
    }

    /// EA match vector for one text base inside one block. Text `A` also
    /// pairs with pattern `G`, text `T` with pattern `C`; text `N` never
    /// matches.
    #[inline]
    fn eq_vec(&self, blk: usize, tcode: u8) -> u128 {
        if tcode > 3 {
            return 0;
        }
        let mut eq = self.peq[blk][tcode as usize];
        if tcode == 0 {
            eq |= self.peq[blk][2];
        } else if tcode == 3 {
            eq |= self.peq[blk][1];
        }
        eq
    }

    /// Minimum leading-column-free distance over all text positions.
    fn infix_min(&self, text: &[u8]) -> u32 {
        let mut pv = self.init_pv;
        let mut mv = [0u128; 2];
        let mut score = self.m as i32;
        let mut best = self.m as u32;
        for &b in text {
            let tcode = code5(b);
            let mut hin: i8 = 0;
            let mut cin = false;
            for blk in 0..self.blocks {
                let eq = self.eq_vec(blk, tcode);
                let (h, c) = calculate_block(eq, &mut pv[blk], &mut mv[blk], hin, cin, self.top[blk]);
                hin = h;
                cin = c;
            }
            score += hin as i32;
            if score >= 0 && (score as u32) < best {
                best = score as u32;
            }
        }
        best
    }

    /// Global distance with forced ends.
    fn global(&self, text: &[u8]) -> u32 {
        let mut pv = self.init_pv;
        let mut mv = [0u128; 2];
        let mut score = self.m as i32;
        for &b in text {
            let tcode = code5(b);
            let mut hin: i8 = 1;
            let mut cin = false;
            for blk in 0..self.blocks {
                let eq = self.eq_vec(blk, tcode);
                let (h, c) = calculate_block(eq, &mut pv[blk], &mut mv[blk], hin, cin, self.top[blk]);
                hin = h;
                cin = c;
            }
            score += hin as i32;
        }
        score.max(0) as u32
    }
}

/// One block-column step of the Myers recurrence. `hin` in {-1, 0, 1} is
/// the horizontal delta entering the block's top row; `cin` is the
/// ripple-carry leaving the previous block's match-propagation adder.
/// Returns the horizontal delta leaving the block's bottom row (`hout`)
/// and this block's own adder carry (`cout`) — both must be chained into
/// the next block (dropping `cout` corrupts multi-block columns).
#[inline]
fn calculate_block(eq: u128, pv: &mut u128, mv: &mut u128, hin: i8, cin: bool, top: u128) -> (i8, bool) {
    let xv = eq | *mv;
    let (sum, cout) = (eq & *pv).carrying_add(*pv, cin);
    let xh = (sum ^ *pv) | xv;
    let mut ph = *mv | !(xh | *pv);
    let mut mh = *pv & xh;
    let mut hout: i8 = 0;
    if ph & top != 0 {
        hout = 1;
    } else if mh & top != 0 {
        hout = -1;
    }
    ph <<= 1;
    mh <<= 1;
    if hin == 1 {
        ph |= 1;
    } else if hin == -1 {
        mh |= 1;
    }
    *pv = mh | !(xv | ph);
    *mv = ph & xv;
    (hout, cout)
}
